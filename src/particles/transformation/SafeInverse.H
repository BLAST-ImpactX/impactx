/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_SAFE_INVERSE_H
#define IMPACTX_SAFE_INVERSE_H

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx::transformation
{
    /** Compute inverse denominator with numerical stability for small values.
     *
     * @param denom denominator value
     * @param tol tolerance threshold for edge case handling (~1000x machine epsilon)
     * @return inverse of denominator (1.0 / denom) with edge case handling
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::ParticleReal
    safe_inverse (
        amrex::ParticleReal const denom,
        amrex::ParticleReal const tol =
#ifdef AMREX_SINGLE_PRECISION_PARTICLES
            1.0e-6
#else
            1.0e-12
#endif
    )
    {
        using namespace amrex::literals;
        return std::abs(denom) > tol ?
            1.0_prt / denom :
            // Handle edge case: use linearized approximation for small denominator
            1.0_prt / (denom + (denom >= 0.0_prt ? tol : -tol));
    }

} // namespace impactx::transformation

#endif // IMPACTX_SAFE_INVERSE_H
