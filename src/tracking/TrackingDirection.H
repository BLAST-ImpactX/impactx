/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_TRACKING_DIRECTION_H
#define IMPACTX_TRACKING_DIRECTION_H


namespace impactx
{
    /** Direction in which the lattice is traversed.
     *
     * ``Forward`` is ordinary tracking (entrance to exit). ``Backward`` traverses
     * the periods, elements and space-charge slices in reverse, inverts each element's map,
     * and applies the per-slice operations in mirrored order, i.e. time-symmetric back-tracking.
     */
    enum class TrackingDirection { Forward, Backward };

} // namespace impactx

#endif // IMPACTX_TRACKING_DIRECTION_H
