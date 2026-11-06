/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_HELPER_EXACTDRIFTSLOPES_H
#define IMPACTX_ELEMENTS_HELPER_EXACTDRIFTSLOPES_H

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx::elements
{
    /** Slopes of the exact drift map, in a region without field
     *
     * @see exact_drift_slopes
     */
    template<typename T_Real>
    struct ExactDriftSlopes
    {
        T_Real inv_pz;     //! 1/pz, so that dx/ds = px/pz and dy/ds = py/pz
        T_Real inv_pz_m1;  //! 1/pz - 1
        T_Real dt_ds;      //! dt/ds = (1/beta - pt)/pz - 1/beta
    };

    /** Slopes of the exact drift map of a particle, in a region without field
     *
     * The longitudinal momentum is
     *   pz^2 = (1/beta - pt)^2 - 1/(beta*gamma)^2 - px^2 - py^2 = 1 - eps,
     *   eps = 2 pt/beta - pt^2 + px^2 + py^2,
     * since 1/beta^2 - 1/(beta*gamma)^2 = 1. The deviations from the reference particle,
     * 1/pz - 1 = eps / (pz (1 + pz)) and dt/ds, are computed from eps, without subtracting
     * terms of order 1/beta from each other: particles that move with the reference
     * particle (px = py = pt = 0) keep x, y and t exactly.
     *
     * @param px particle momentum in x
     * @param py particle momentum in y
     * @param pt particle momentum in t
     * @param ibeta 1/beta of the reference particle
     * @return 1/pz, 1/pz - 1 and dt/ds
     */
    template<typename T_Real>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    ExactDriftSlopes<T_Real>
    exact_drift_slopes (
        T_Real const & px,
        T_Real const & py,
        T_Real const & pt,
        amrex::ParticleReal ibeta
    )
    {
        using namespace amrex::literals; // for _rt and _prt
        using namespace std; // for cmath(float)

        T_Real const eps = 2_prt * ibeta * pt - pt * pt + px * px + py * py;
        T_Real const pz = sqrt(1_prt - eps);
        T_Real const inv_pz = 1_prt / pz;
        T_Real const eps_1pz = eps / (1_prt + pz);

        return {
            inv_pz,
            eps_1pz * inv_pz,
            inv_pz * (ibeta * eps_1pz - pt)
        };
    }

} // namespace impactx::elements

#endif // IMPACTX_ELEMENTS_HELPER_EXACTDRIFTSLOPES_H
