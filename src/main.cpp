/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell, Ji Qiang
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactX.H"
#include "initialization/InitAMReX.H"


int main(int argc, char* argv[])
{
    using namespace impactx;

    // although ImpactX' init_grids will call this if not done before, we call
    // it here so users can pass command line arguments
    initialization::default_init_AMReX(argc, argv);

    {
        my_run();
    }

}
