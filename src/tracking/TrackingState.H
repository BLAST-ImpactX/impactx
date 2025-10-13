/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_TRACKING_STATE_H
#define IMPACTX_TRACKING_STATE_H

#include "elements/All.H"


namespace impactx
{
    /** This class wraps a few states during tracking loops.
     *
     * These states are useful to know in hooks (Python callback functions).
     */
    struct TrackingState
    {
        /** For tracking hooks/callbacks, a global step of the simulation
         *
         * A state of internal simulation steps, increments also for space charge slice steps in elements.
         * We start in "step 0" (initial state).
         */
        int m_step = 0;

        /** For tracking hooks/callbacks, the period in the lattice (e.g., turn or channel period) */
        int m_period = 0;

        /** For tracking hooks/callbacks, the current lattice element */
        elements::KnownElements* m_element = &empty_element;

        void set_no_element () { m_element = &empty_element; }

    private:
        static inline elements::KnownElements empty_element = elements::Empty{};
    };

} // namespace impactx

#endif // IMPACT_X_H
