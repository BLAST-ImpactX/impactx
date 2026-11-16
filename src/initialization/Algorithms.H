/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
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

#include <string>


namespace impactx
{
    /** Space Charge: Implemented algorithms */
    AMREX_ENUM(SpaceChargeAlgo,
        False,    /**< Disabled: no space charge calculation */
        True_3D,  /**< 3D beam distribution */
        Gauss3D,  /**< Assume a 3D Gaussian beam distribution */
        Gauss2p5D,  /**< Assume a transverse 2D Gaussian beam distribution */
        True_2D,   /**< Averaged 2D transverse beam distribution */
        True_2p5D   /**< Averaged 2D transverse beam distribution with a current along s */
    );

    /** Return the currently active space charge algorithm */
    SpaceChargeAlgo
    get_space_charge_algo ();

    /** A user-friendly string
     *
     * Returns a string that is similar to the spelling in user input options.
     */
    std::string
    to_string (SpaceChargeAlgo sca);

    /** CSR: Implemented algorithms */
    AMREX_ENUM(CSR_Algo,
        False,    /**< Disabled: no csr calculation */
        Ultrarelativistic,  /**< ultrarelativistic steady-state 1D CSR model */
        SteadyState  /**< finite-gamma steady-state 1D CSR model */
    );

    /** Return the currently active csr algorithm */
    CSR_Algo
    get_csr_algo ();

    /** A user-friendly string
     *
     * Returns a string that is similar to the spelling in user input options.
     */
    std::string
    to_string (CSR_Algo csra);

} // namespace impactx

#endif // IMPACTX_ALGORITHMS_H
