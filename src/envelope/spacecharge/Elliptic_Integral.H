/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELLIPTIC_INTEGRAL_H
#define IMPACTX_ELLIPTIC_INTEGRAL_H

#include <AMReX_Extension.H>  // for AMREX_RESTRICT
#include <AMReX_REAL.H>


namespace impactx::envelope::spacecharge
{
    /** Function to return the Carlson elliptic integral denoted
     *  by RD(x,y,z).  This is a reimplementation of
     *  the C++ function written by John Burkardt (GNU LGPL license):
     *
     *  https://people.math.sc.edu/Burkardt/cpp_src
     *        /elliptic_integral/elliptic_integral.html
     *
     *  which in turn is based on:
     *
     * Bille Carlson,
     * Computing Elliptic Integrals by Duplication,
     * Numerische Mathematik,  
     * Volume 33, 1979, pages 1-16.
     *
     * Bille Carlson, Elaine Notis,
     * Algorithm 577, Algorithms for Incomplete Elliptic Integrals,
     * ACM Transactions on Mathematical Software,
     * Volume 7, Number 3, pages 398-403, September 1981.
     *
     * @param[in] x first argument
     * @param[in] y second argument
     * @param[in] z third argument
     * @param[in] errtol error tolerance
     * @returns the elliptic integral RD(x,y,z)
     */
    amrex::ParticleReal
    Elliptic_RD (
        amrex::ParticleReal x, 
        amrex::ParticleReal y,
        amrex::ParticleReal z,
        amrex::ParticleReal errtol
    );

} // namespace impactx

#endif // IMPACTX_ELLIPTIC_INTEGRAL_H
