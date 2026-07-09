/* Copyright 2022-2025 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_TRACKING_PARTICLES_H
#define IMPACTX_TRACKING_PARTICLES_H

#include "elements/All.H"
#include "elements/mixin/accessors.H"
#include "particles/ImpactXParticleContainer.H"
#include "tracking/TrackingDirection.H"
#include "tracking/TrackingState.H"

#include <AMReX_BLProfiler.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>

#include <list>
#include <string>


namespace impactx
{
    /** Traverse the lattice for particle (beam + reference) tracking.
     *
     *  This function performs the nested lattice traversal (looping over periods,
     *  then over the beamline elements in each period, then over the
     *  collective-effect slices within each element). It also maintains the
     *  user-defined hooks, the reference-particle edge update and the global step
     *  counter.
     *
     *  The caller supplies the per-slice physics as a @p collective_kicks K
     *  (space charge, CSR, ISR, ...), an @p element_transport M (the external-field
     *  map of one slice of an element) and a @p slice_diagnostics book-keeping step
     *  (boundary conditions, lost particles, slice diagnostics) run once per slice.
     *
     *  When collective effects are active (@p collective_effects), the kick ``K`` and
     *  transport ``M`` of each finite-length (``Thick``) element are composed with a
     *  second-order, time-symmetric drift-outer (MKM) Strang split per slice:
     *      ``M(ds/2) K(ds) M(ds/2)``
     *  half the slice transport, the full collective kick at the slice midpoint, then
     *  the other half. The half-slice transports are realized by temporarily doubling
     *  the element's ``nslice``; the reference-particle edge is set once per element, so
     *  an s-dependent element map integrates continuously across the two half-maps (the
     *  second half starts where the first ended). This applies the expensive collective
     *  solve once per slice (``nslice`` per element).
     *
     *  Thin (zero-length) elements, and all elements in pure-optics runs (no collective
     *  effects), apply the element transport once per slice with no kick, unchanged.
     *
     *  Hooks, the step counter and the tracking-state element pointer are maintained on
     *  ``Forward`` tracking only. The reference-particle edge is set in both directions,
     *  since s-dependent (soft-edge) element maps integrate their field profile relative to
     *  it. ``Backward`` tracking mirrors the traversal (periods, elements and slices in
     *  reverse) and applies the same MKM split with each element map reversed (``M -> M^-1``)
     *  and each kick inverted (``K -> K^-1``, via the negated ``ds`` that makes ``slice_ds``
     *  negative) — i.e. the exact time-reverse of forward tracking.
     *
     *  The number of periods (``lattice.periods``) and verbosity (``impactx.verbose``)
     *  are read from the inputs. The step counter is advanced from its current value,
     *  so the caller is expected to reset it (to 0) before a forward run.
     *
     *  @tparam HookFn callable ``void (std::string const & name)`` invoking a named hook
     *  @tparam CollectiveKicks callable
     *    ``void (elements::KnownElements & element, amrex::ParticleReal slice_ds)``
     *    applying the collective-effect kick of length @p slice_ds
     *  @tparam ElementTransport callable
     *    ``void (elements::KnownElements & element, int step, int period)``
     *    applying one slice of the element's external-field transport map (no book-keeping)
     *  @tparam SliceDiagnostics callable ``void (int step)`` running the once-per-slice
     *    book-keeping (boundary conditions, lost particles, slice diagnostics)
     *  @param[in,out] lattice the beamline elements to traverse
     *  @param[in,out] pc beam particle container (reference particle / edge update)
     *  @param[in,out] tracking_state tracking state shared with hooks (element, step, period);
     *                 its ``m_direction`` selects the traversal direction
     *  @param[in] collective_effects whether any collective kick is active this run; only then
     *                 are thick elements Strang-split around a mid-slice kick
     *  @param[in] call_hook registry for user-specified callback hooks
     *  @param[in] collective_kicks applied collective effect kicks (SC, CSR, ISR, etc.), per-slice
     *  @param[in] element_transport applies the external-field map of one element slice
     *  @param[in] slice_diagnostics once-per-slice book-keeping and diagnostics
     */
    template <typename HookFn, typename CollectiveKicks, typename ElementTransport, typename SliceDiagnostics>
    void track_lattice_particles (
        std::list<elements::KnownElements> & lattice,
        ImpactXParticleContainer & pc,
        TrackingState & tracking_state,
        bool collective_effects,
        HookFn const & call_hook,
        CollectiveKicks && collective_kicks,
        ElementTransport && element_transport,
        SliceDiagnostics && slice_diagnostics
    )
    {
        bool const forward = tracking_state.forward();

        // periods through the lattice (turns/cycles) and verbosity
        int num_periods = 1;
        amrex::ParmParse("lattice").queryAddWithParser("periods", num_periods);
        int verbose = 1;
        amrex::ParmParse("impactx").queryAddWithParser("verbose", verbose);

        int & step = tracking_state.m_step;

        // visit one element for all of its space-charge slices
        auto visit_element = [&] (
            elements::KnownElements & element_variant,
            int period
        )
        {
            // Set the reference-particle edge at the entrance of this element. Needed in
            // both directions by s-dependent (soft-edge) maps, which integrate their on-axis
            // field profile relative to this edge; under ``Backward`` tracking the reversed
            // map traverses that profile in the opposite z-direction (see the element reverse()).
            pc.SetRefParticleEdge();

            if (forward)
            {
                // book-keeping: currently visited element reference
                tracking_state.m_element = &element_variant;

                // optional, user-defined function call
                call_hook("before_element");
            }
            else
            {
                // back-tracking: invert this element's map in place so the element
                // transport below applies the inverse map M^-1.
                elements::reverse(element_variant);
            }

            // number of slices used for the application of space charge
            int const nslice = elements::nslice(element_variant);
            amrex::ParticleReal const slice_ds = elements::slice_ds(element_variant); // in meters

            // Strang-split a thick element's transport ``M`` around a mid-slice collective
            // kick ``K`` (drift-outer / MKM, #846 / #920) only where it is meaningful: for
            // finite-length elements while collective effects are active. This second-order,
            // time-symmetric split applies the expensive collective solve once per slice.
            // Under ``Backward`` tracking ``slice_ds`` is negative (the map is reversed), so
            // the same composition is the exact inverse, M^-1(ds/2) K^-1(ds) M^-1(ds/2).
            // Thin elements and pure-optics runs apply the transport once, unchanged.
            bool const split = collective_effects
                            && slice_ds != amrex::ParticleReal(0)
                            && elements::is_thick(element_variant);

            // sub-steps for space charge within the element
            for (int s = 0; s < nslice; ++s)
            {
                BL_PROFILE("ImpactX::evolve::slice_step");
                int const slice_step = forward ? s : (nslice - 1 - s);

                if (forward)
                {
                    step++;
                    if (verbose > 0) {
                        amrex::Print() << "\n++++ Starting step=" << step
                                       << " slice_step=" << slice_step;
                    }

                    // optional, user-defined function call
                    call_hook("before_slice");
                }

                if (split)
                {
                    // drift-outer (MKM): M(ds/2) K(ds) M(ds/2). Doubling ``nslice`` makes each
                    // element_transport advance ds/2; the reference-particle edge is fixed for
                    // the element, so an s-dependent map continues across the two half-maps.
                    // (The element is reversed in place for back-tracking, applying M^-1.)
                    elements::ScopedNslice const half_slices(element_variant, 2 * nslice);
                    element_transport(element_variant, step, period);
                    collective_kicks(element_variant, slice_ds);
                    element_transport(element_variant, step, period);
                }
                else
                {
                    // thin element or no collective effects: apply the transport once, no kick
                    element_transport(element_variant, step, period);
                }

                // once-per-slice book-keeping and diagnostics
                slice_diagnostics(step);
            }

            if (forward)
            {
                // optional, user-defined function call
                call_hook("after_element");
            }
            else
            {
                // restore the forward element map (reverse() is an involution)
                elements::reverse(element_variant);
            }
        };

        // loop over channel periods or ring turns
        for (int p = 0; p < num_periods; ++p)
        {
            int const period = forward ? p : (num_periods - 1 - p);

            if (forward)
            {
                // book-keeping update
                tracking_state.m_period = period;
                tracking_state.m_element = &lattice.front();

                // optional, user-defined function call
                call_hook("before_period");
            }

            // loop over lattice elements
            if (forward)
            {
                for (auto & element_variant : lattice)
                {
                    visit_element(element_variant, period);
                }
            }
            else
            {
                for (auto it = lattice.rbegin(); it != lattice.rend(); ++it)
                {
                    visit_element(*it, period);
                }
            }

            // optional, user-defined function call
            if (forward) { call_hook("after_period"); }
        }

        if (forward)
        {
            // avoid dangling references if users manipulate the lattice
            tracking_state.set_no_element();
        }
    }
} // namespace impactx

#endif // IMPACTX_TRACKING_PARTICLES_H
