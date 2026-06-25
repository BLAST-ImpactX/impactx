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
#include "particles/ImpactXParticleContainer.H"
#include "tracking/TrackingDirection.H"
#include "tracking/TrackingState.H"

#include <AMReX_BLProfiler.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>

#include <list>
#include <string>
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
     *  (space charge, CSR, ISR, ...) and an @p element_push M
     *  (the external-field of an element). Both are combined in a second-order
     *  Strang split:
     *      ``K(ds/2) [M(ds) K(ds)]^(nslice-1) M(ds) K(ds/2)``
     *  i.e., a half-length collective kick at each element edge and full-length
     *  interior kicks between the per-slice transport maps. This is second-order
     *  accurate in the slice step size and time-symmetric. It costs one extra
     *  collective kick per element over a first-order ``[K(ds) M(ds)]^nslice``
     *  scheme.
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
     *  @tparam ElementPush callable
     *    ``void (elements::KnownElements & element, int step, int period)``
     *    applying one slice of the element's external-field transport map
     *  @param[in,out] lattice the beamline elements to traverse
     *  @param[in,out] pc beam particle container (reference particle / edge update)
     *  @param[in,out] tracking_state tracking state shared with hooks (element, step, period);
     *                 its ``m_direction`` selects the traversal direction
     *  @param[in] call_hook registry for user-specified callback hooks
     *  @param[in] collective_kicks applied collective effect kicks (SC, CSR, ISR, etc.), per-slice
     *  @param[in] element_push applied the external-field of an element, per-slice
     */
    template <typename HookFn, typename CollectiveKicks, typename ElementPush>
    void track_lattice_particles (
        std::list<elements::KnownElements> & lattice,
        ImpactXParticleContainer & pc,
        TrackingState & tracking_state,
        HookFn const & call_hook,
        CollectiveKicks && collective_kicks,
        ElementPush && element_push
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

            // Second-order Strang splitting of the collective-effect kick ``K``
            // and the per-slice element transport ``M`` over the ``nslice``
            // slices of the element (#846 / #920):
            //     K(ds/2) [M(ds) K(ds)]^(nslice-1) M(ds) K(ds/2)
            // i.e. a half-length kick at each element edge and a full-length
            // interior kick between consecutive transport maps. Each interior
            // ``K(ds)`` is the exact merge of two adjacent half-kicks
            // ``K(ds/2) K(ds/2)``: a collective kick only changes the momenta,
            // not the positions it gathers the self-field from, so successive
            // kicks compose exactly, K(a) K(b) = K(a+b). This second-order scheme
            // therefore costs only one extra field solve per element over the
            // first-order ``[K(ds) M(ds)]^nslice`` scheme.
            amrex::ParticleReal const half_ds = slice_ds / amrex::ParticleReal(2);

            // A kick of zero length (e.g. thin elements such as monitors) changes
            // neither the momenta nor the gathered positions, so it is skipped
            // exactly; this also avoids redundant field solves.
            bool const apply_kicks = slice_ds > amrex::ParticleReal(0);

            // entrance half-kick K(ds/2)
            if (apply_kicks) { collective_kicks(element_variant, half_ds); }

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

                // interior full kick K(ds) between consecutive transport maps
                if (apply_kicks && s > 0) { collective_kicks(element_variant, slice_ds); }

                // per-slice element transport map M(ds)
                // (the element is reversed in place for back-tracking, so this
                //  applies the inverse map on the ``Backward`` pass)
                element_push(element_variant, step, period);
            }

            // exit half-kick K(ds/2)
            if (apply_kicks) { collective_kicks(element_variant, half_ds); }

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
