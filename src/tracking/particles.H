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
#include "elements/mixin/thick.H"
#include "particles/ImpactXParticleContainer.H"
#include "tracking/TrackingDirection.H"
#include "tracking/TrackingState.H"

#include <AMReX_BLProfiler.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>

#include <list>
#include <string>
#include <type_traits>
#include <variant>


namespace impactx
{
    /** Traverse the lattice for particle (beam + reference) tracking, applying the
     *  per-slice physics as two composable operations.
     *
     *  This function performs the nested lattice traversal (looping over periods,
     *  then over the beamline elements in each period, then over the
     *  collective-effect slices within each element). It also maintains the
     *  user-defined hooks, the reference-particle edge update and the global step
     *  counter.
     *  The caller supplies the per-slice physics as a @p collective_kicks K
     *  (space charge, CSR, ISR, ...), an @p element_transport M (the external-field
     *  map of one slice of an element) and a @p slice_diagnostics book-keeping step
     *  (boundary conditions, lost particles, slice diagnostics) run once per slice.
     *  ``K`` and ``M`` are combined per element with a second-order Strang split, in
     *  one of two equivalent forms:
     *
     *   - kick-outer (KMK), the default and the only correct form for general elements:
     *       ``K(ds/2) [M(ds) K(ds)]^(nslice-1) M(ds) K(ds/2)``
     *     a half-length kick at each element edge, full-length interior kicks between
     *     the per-slice transport maps. Costs ``nslice+1`` collective solves per element.
     *
     *   - drift-outer (MKM), used only where the transport is an exact constant-field
     *     flow (\see elements::mixin::Thick::is_constant_field_flow, so
     *     ``M(ds/2) M(ds/2) == M(ds)`` exactly):
     *       ``M(ds/2) [K(ds) M(ds)]^(nslice-1) K(ds) M(ds/2)``
     *     a half-length transport at each element edge, full-length interior transports
     *     between the per-slice kicks. The two half-transports are realized by temporarily
     *     doubling the element's ``nslice``. Same second order, but it keeps the expensive
     *     collective solve at ``nslice`` (one fewer per element than KMK) by paying an
     *     extra, cheap transport instead.
     *
     *  Both forms are second-order accurate in the slice step size and time-symmetric.
     *
     *  Hooks, the reference-particle edge update, the step counter and the
     *  tracking-state element pointer are maintained on ``Forward`` tracking only.
     *  The ``Backward`` pass is a bare mirrored traversal that applies the same
     *  symmetric composition with the inverted, reversed element maps.
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

        // is this element eligible for the cheaper drift-outer (MKM) Strang split?
        //   (only finite-length ``Thick`` elements carry the is_constant_field_flow property)
        auto mkm_eligible = [] (elements::KnownElements const & ev) -> bool
        {
            return std::visit([] (auto const & element) -> bool {
                using E = std::remove_cv_t<std::remove_reference_t<decltype(element)>>;
                if constexpr (std::is_base_of_v<elements::mixin::Thick, E>) {
                    return E::is_constant_field_flow;
                } else {
                    return false;
                }
            }, ev);
        };

        // set the per-element slice count (only finite-length ``Thick`` elements carry one);
        // used to temporarily halve the transport length for the drift-outer split
        auto set_element_nslice = [] (elements::KnownElements & ev, int n)
        {
            std::visit([n] (auto && element) {
                using E = std::remove_cv_t<std::remove_reference_t<decltype(element)>>;
                if constexpr (std::is_base_of_v<elements::mixin::Thick, E>) {
                    element.m_nslice = n;
                }
            }, ev);
        };

        // visit one element for all of its space-charge slices
        auto visit_element = [&] (
            elements::KnownElements & element_variant,
            int period
        )
        {
            if (forward)
            {
                // update element edge of the reference particle
                pc.SetRefParticleEdge();

                // book-keeping: currently visited element reference
                tracking_state.m_element = &element_variant;

                // optional, user-defined function call
                call_hook("before_element");
            }
            else
            {
                // back-tracking: invert this element's map in place so the element
                // transport below applies the inverse map.
                std::visit([](auto && element){ element.reverse(); }, element_variant);
            }

            // number of slices used for the application of space charge
            int nslice = 1;
            amrex::ParticleReal slice_ds = 0; // in meters
            std::visit([&nslice, &slice_ds](auto && element) {
                nslice = element.nslice();
                slice_ds = element.ds() / nslice;
            }, element_variant);

            // Second-order Strang split of the collective kick ``K`` and the per-slice
            // transport ``M`` (#846 / #920). A kick of zero length (thin elements such as
            // monitors) changes neither momenta nor the gathered positions, so it is
            // skipped exactly (apply_kicks), avoiding redundant field solves.
            bool const apply_kicks = slice_ds > amrex::ParticleReal(0);
            amrex::ParticleReal const half_ds = slice_ds / amrex::ParticleReal(2);

            // Where the transport is an exact constant-field flow (M(ds/2) M(ds/2)==M(ds)),
            // use the drift-outer (MKM) form: it keeps the expensive collective solve at
            // ``nslice`` per element instead of ``nslice+1``. Otherwise the uniform
            // kick-outer (KMK) form, which never sub-divides the element transport.
            bool const drift_outer = apply_kicks && mkm_eligible(element_variant);

            // kick-outer (KMK) entrance half-kick K(ds/2)
            if (apply_kicks && !drift_outer) { collective_kicks(element_variant, half_ds); }

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

                if (drift_outer)
                {
                    // drift-outer (MKM): M(ds/2) K(ds) M(ds/2). Doubling ``nslice`` makes
                    // each element_transport advance ds/2; the half-transports of adjacent
                    // slices merge into M(ds), so K is solved once per slice. The element is
                    // reversed in place for back-tracking, applying M^-1 on ``Backward``.
                    set_element_nslice(element_variant, 2 * nslice);
                    element_transport(element_variant, step, period);
                    collective_kicks(element_variant, slice_ds);
                    element_transport(element_variant, step, period);
                    set_element_nslice(element_variant, nslice);
                }
                else
                {
                    // kick-outer (KMK): interior full kick K(ds) before each transport;
                    // the entrance/exit half-kicks are applied outside the loop.
                    if (apply_kicks && s > 0) { collective_kicks(element_variant, slice_ds); }
                    element_transport(element_variant, step, period);
                }

                // once-per-slice book-keeping and diagnostics
                slice_diagnostics(step);
            }

            // kick-outer (KMK) exit half-kick K(ds/2)
            if (apply_kicks && !drift_outer) { collective_kicks(element_variant, half_ds); }

            if (forward)
            {
                // optional, user-defined function call
                call_hook("after_element");
            }
            else
            {
                // restore the forward element map (reverse() is an involution)
                std::visit([](auto && element){ element.reverse(); }, element_variant);
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
