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
     *  The caller supplies the per-slice physics as three composable operations: the
     *  collective kick @p collective_kicks ``K`` (space charge, CSR, ISR, ...), the
     *  element transport @p element_push ``M`` (the external-field map of one
     *  element slice) and the @p slice_diagnostics book-keeping (boundary conditions,
     *  lost particles, slice diagnostics), which runs once per slice.
     *
     *  With collective effects active (@p collective_effects), each slice of an element
     *  that can be subdivided (@see elements::can_slice) is composed as the second-order,
     *  time-symmetric Strang split
     *      ``M(ds/2) K(ds) M(ds/2)``
     *  half of the slice transport, the collective kick at the slice midpoint, then the
     *  other half. The half transports are realized by temporarily doubling the element's
     *  ``nslice``, which simply halves the ds per slice. The expensive collective solve
     *  still runs only once per slice.
     *  Elements that cannot be subdivided, and all elements once the split is turned off
     *  (@p strang_split), are composed to first order instead: ``K(ds) M(ds)`` forward and
     *  ``M(ds) K(ds)`` backward.
     *  Zero-length elements and runs without collective effects apply ``M`` alone, since a
     *  kick over zero length leaves the momenta unchanged. Typically, thin elements are
     *  embedded in drift elements, which experience the collective effect kicks.
     *
     *  The reference-particle edge is set once per element and in both directions:
     *  s-dependent (soft-edge) maps integrate their on-axis field profile relative to it,
     *  so the second half-map continues where the first one ended.
     *
     *  ``Backward`` tracking mirrors the traversal (periods, elements and slices in
     *  reverse), reverses each element map (``M -> M^-1``) and inverts each kick (the
     *  negated ``ds`` makes ``slice_ds`` negative). This is the exact time-reverse of
     *  forward tracking. Hooks, the step counter and the tracking-state element pointer
     *  are maintained on ``Forward`` tracking only.
     *
     *  The number of periods (``lattice.periods``) and verbosity (``impactx.verbose``)
     *  are read from the inputs. The step counter is advanced from its current value,
     *  so the caller is expected to reset it (to 0) before a forward run.
     *
     *  @tparam HookFn callable ``void (std::string const & name)`` invoking a named hook
     *  @tparam CollectiveKicks callable
     *    ``void (elements::KnownElements & element, amrex::ParticleReal slice_ds)``
     *    applying the collective-effect kick over the length @p slice_ds
     *  @tparam ElementPush callable
     *    ``void (elements::KnownElements & element, int step, int period)``
     *    applying one slice of the element's external-field transport map (no book-keeping)
     *  @tparam SliceDiagnostics callable ``void (int step, int period)`` running the
     *    once-per-slice book-keeping (boundary conditions, lost particles, slice diagnostics)
     *  @param[in,out] lattice the beamline elements to traverse
     *  @param[in,out] pc beam particle container (reference particle / edge update)
     *  @param[in,out] tracking_state tracking state shared with hooks (element, step, period);
     *                 its ``m_direction`` selects the traversal direction
     *  @param[in] collective_effects whether any collective kick is active in this run
     *  @param[in] strang_split whether to Strang-split the collective kicks (second order),
     *                 false composes kick and transport to first order
     *  @param[in] call_hook registry for user-specified callback hooks
     *  @param[in] collective_kicks applies the collective effect kicks (SC, CSR, ISR, ...)
     *  @param[in] element_push applies the external-field map of one element slice
     *  @param[in] slice_diagnostics once-per-slice book-keeping and diagnostics
     */
    template <
        typename HookFn, typename CollectiveKicks,
        typename ElementPush, typename SliceDiagnostics
    >
    void track_lattice_particles (
        std::list<elements::KnownElements> & lattice,
        ImpactXParticleContainer & pc,
        TrackingState & tracking_state,
        bool collective_effects,
        bool strang_split,
        HookFn const & call_hook,
        CollectiveKicks && collective_kicks,
        ElementPush && element_push,
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
            // set the reference-particle edge to the entrance of this element: s-dependent
            // (soft-edge) maps integrate their on-axis field profile relative to it
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
                // backward tracking: invert this element's map in place so the element
                // transport below applies the inverse map M^-1.
                elements::reverse(element_variant);
            }

            // number of slices used for the application of space charge
            int const nslice = elements::nslice(element_variant);
            amrex::ParticleReal const slice_ds = elements::slice_ds(element_variant); // in meters

            // zero-length elements receive no kick: it would leave the momenta unchanged
            bool const kick = collective_effects && slice_ds != amrex::ParticleReal(0);

            // Strang-split the transport ``M`` around the kick ``K`` wherever
            // the element can be subdivided into the two half-maps
            bool const split = kick && strang_split && elements::can_slice(element_variant);

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
                    // M(ds/2) K(ds) M(ds/2): the doubled slice count halves each transport
                    elements::ScopedNslice const half_slices(element_variant, 2 * nslice);
                    element_push(element_variant, step, period);
                    collective_kicks(element_variant, slice_ds);
                    element_push(element_variant, step, period);
                }
                else if (kick)
                {
                    // finite length, but the element cannot be subdivided:
                    // first-order composition
                    if (forward)
                    {
                        collective_kicks(element_variant, slice_ds);
                        element_push(element_variant, step, period);
                    }
                    else
                    {
                        element_push(element_variant, step, period);
                        collective_kicks(element_variant, slice_ds);
                    }
                }
                else
                {
                    // zero-length slice or no collective effects: transport only
                    element_push(element_variant, step, period);
                }

                // once-per-slice book-keeping and diagnostics
                slice_diagnostics(step, period);
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
