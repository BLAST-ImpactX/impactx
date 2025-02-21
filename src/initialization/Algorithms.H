/* Copyright 2022-2025 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ALGORITHMS_H
#define IMPACTX_ALGORITHMS_H

#include <AMReX_Enum.H>


namespace impactx
{
    AMREX_ENUM(SpaceChargeAlgo,
        False,
        True_3D,
        True_2D
    );

    SpaceChargeAlgo
    get_space_charge_algo ();

} // namespace impactx

#endif // IMPACTX_ALGORITHMS_H
