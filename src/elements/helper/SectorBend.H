/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_HELPER_SECTORBEND_H
#define IMPACTX_ELEMENTS_HELPER_SECTORBEND_H

#include "Stumpff.H"
#include "particles/CovarianceMatrix.H"

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>
#include <AMReX_SIMD.H>

#include <cmath>


namespace impactx::elements
{
    /** Geometry of a step along the reference orbit of a sector bend
     *
     * With the curvature h, the bend radius rc = 1/h is not needed: all lengths are
     * finite for any curvature, including h = 0.
     *
     * @see sector_bend
     */
    struct SectorBend
    {
        amrex::ParticleReal h;              //! curvature of the reference orbit, 1/rc [1/m]
        amrex::ParticleReal L;              //! arc length of the reference orbit [m]
        amrex::ParticleReal sin_phi;        //! sin(phi), with the bend angle phi = h L
        amrex::ParticleReal cos_phi;        //! cos(phi)
        amrex::ParticleReal rc_sin_phi;     //! rc sin(phi) [m]
        amrex::ParticleReal rc_versin_phi;  //! rc (1 - cos(phi)) [m]
    };

    /** Geometry of a step along the reference orbit of a sector bend
     *
     * @param h curvature of the reference orbit, 1/rc [1/m]
     * @param L arc length of the reference orbit [m]
     * @return the trigonometric functions of the bend angle phi = h L
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    SectorBend
    sector_bend (
        amrex::ParticleReal h,
        amrex::ParticleReal L
    )
    {
        amrex::ParticleReal const phi = h * L;
        auto const c = stumpff(phi * phi);

        return {h, L, phi * c.c1, c.c0, L * c.c1, L * phi * c.c2};
    }

    /** Exact map of a particle through a step of a sector bend with a uniform field
     *
     * The map of D. L. Bruhwiler et al., Proc. EPAC 98, pp. 1171-1173 (1998), in which
     * the particle moves on a circle of radius rc p_perp in the bend plane:
     *
     *   pxf = px cos(phi) + (pzi - 1 - h x) sin(phi),
     *   xf  = -rc + (rc + x) cos(phi) + rc (pzf + px sin(phi) - pzi cos(phi)),
     *   yf  = y + rc theta py,
     *   tf  = t - rc theta (pt - 1/beta) - L/beta,
     *
     * with theta = phi + asin(px/p_perp) - asin(pxf/p_perp) the angle of the particle's
     * arc, p_perp^2 = 1 - 2 pt/beta + pt^2 - py^2 and pz^2 = p_perp^2 - px^2. The terms of
     * order rc and L/beta are eliminated analytically:
     *
     *   rc (px - pxf)  = px rc (1 - cos(phi)) - (pzi - 1) rc sin(phi) + x sin(phi),
     *   rc (pzf - pzi) = rc (px - pxf) (px + pxf) / (pzi + pzf),
     *   p_perp^2 sin(theta - phi) = h rc (px - pxf) (px (px + pxf)/(pzi + pzf) + pzi),
     *   p_perp^2 cos(theta - phi) = pzi pzf + px pxf,
     *
     * and rc (theta - phi) follows from the arcsine, or atan2 beyond 45 degrees.
     *
     * The update of each coordinate is then accurate to rounding relative to its change,
     * also for |h x| << 1, and particles that move with the reference particle keep their
     * coordinates exactly.
     *
     * @param[in,out] x particle position in x
     * @param[in,out] px particle momentum in x
     * @param[in,out] y particle position in y
     * @param[in,out] t particle position in t
     * @param py particle momentum in y
     * @param pt particle momentum in t
     * @param bend geometry of the step, with h != 0
     * @param ibeta 1/beta of the reference particle
     * @return mask of the particles outside of the domain of the map, |px| or |pxf| >= p_perp
     */
    template<typename T_Real>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    auto
    exact_sector_bend (
        T_Real & AMREX_RESTRICT x,
        T_Real & AMREX_RESTRICT px,
        T_Real & AMREX_RESTRICT y,
        T_Real & AMREX_RESTRICT t,
        T_Real const & AMREX_RESTRICT py,
        T_Real const & AMREX_RESTRICT pt,
        SectorBend const & bend,
        amrex::ParticleReal ibeta
    )
    {
        using namespace amrex::literals; // for _rt and _prt
        using namespace std; // for cmath(float)

        // p_perp^2 - 1, and pzi - 1 = (p_perp^2 - 1 - px^2) / (1 + pzi)
        T_Real const pperp2_m1 = pt * pt - 2_prt * ibeta * pt - py * py;
        T_Real const pperp2 = 1_prt + pperp2_m1;
        T_Real const pzi = sqrt(pperp2 - px * px);
        T_Real const pzi_m1 = (pperp2_m1 - px * px) / (1_prt + pzi);

        // momentum at the end of the step
        T_Real const pxf = px * bend.cos_phi + (pzi_m1 - bend.h * x) * bend.sin_phi;
        T_Real const pzf = sqrt(pperp2 - pxf * pxf);

        // rc (px - pxf), and (pzf - pzi) / (px - pxf)
        T_Real const rc_dpx = px * bend.rc_versin_phi - pzi_m1 * bend.rc_sin_phi + x * bend.sin_phi;
        T_Real const dpz_dpx = (px + pxf) / (pzi + pzf);

        // theta - phi, with theta the angle of the particle's arc, from p_perp^2 times its
        // sine and cosine: the arcsine is accurate up to 45 degrees, steeper particles
        // (rare) take atan2
        T_Real const sin_pp2 = bend.h * rc_dpx * (px * dpz_dpx + pzi);
        T_Real const cos_pp2 = pzi * pzf + px * pxf;
        auto const steep = cos_pp2 < abs(sin_pp2);
#ifdef AMREX_USE_SIMD
        T_Real sin_dtheta = sin_pp2 / pperp2;
        amrex::simd::stdx::where(steep, sin_dtheta) = 0_prt;
        T_Real dtheta = asin(sin_dtheta);
        if (amrex::simd::stdx::any_of(steep)) {
            amrex::simd::stdx::where(steep, dtheta) = atan2(sin_pp2, cos_pp2);
        }
#else
        T_Real const dtheta = steep ? atan2(sin_pp2, cos_pp2) : asin(sin_pp2 / pperp2);
#endif
        T_Real const rc_dtheta = dtheta / bend.h;

        auto const outside = (pperp2 <= px * px) || (pperp2 <= pxf * pxf);

        // advance position and momentum
        x = x * bend.cos_phi + rc_dpx * dpz_dpx + pzi_m1 * bend.rc_versin_phi
            + px * bend.rc_sin_phi;
        px = pxf;
        y = y + (bend.L + rc_dtheta) * py;
        t = t - bend.L * pt - rc_dtheta * (pt - ibeta);

        return outside;
    }

    /** Linear transport map of a sector bend with a quadrupole gradient
     *
     * The horizontal motion is a harmonic oscillator of strength gx = h^2 + k and the
     * vertical one of strength gy = -k, of any sign. With the Stumpff functions c_j of
     * z = g ds^2, the map entries have no cancellation, also for weak bends, and are
     * finite at h = 0 and g = 0.
     *
     * @param ds arc length of the reference orbit [m]
     * @param h curvature of the reference orbit, 1/rc [1/m]
     * @param k quadrupole strength [1/m^2], horizontally focusing for k > 0
     * @param beta relativistic beta of the reference particle
     * @param ibetagam2 1/(beta gamma)^2 of the reference particle
     * @return 6x6 transport matrix
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Map6x6
    sector_bend_transport_map (
        amrex::ParticleReal ds,
        amrex::ParticleReal h,
        amrex::ParticleReal k,
        amrex::ParticleReal beta,
        amrex::ParticleReal ibetagam2
    )
    {
        amrex::ParticleReal const gx = k + h * h;
        amrex::ParticleReal const gy = -k;
        auto const cx = stumpff(gx * ds * ds);
        auto const cy = stumpff(gy * ds * ds);

        // 1 / (beta rc)
        amrex::ParticleReal const hb = h / beta;

        Map6x6 R = Map6x6::Identity();

        R(1,1) = cx.c0;
        R(1,2) = ds * cx.c1;
        R(1,6) = -ds * ds * cx.c2 * hb;
        R(2,1) = -gx * ds * cx.c1;
        R(2,2) = cx.c0;
        R(2,6) = -ds * cx.c1 * hb;
        R(3,3) = cy.c0;
        R(3,4) = ds * cy.c1;
        R(4,3) = -gy * ds * cy.c1;
        R(4,4) = cy.c0;
        R(5,1) = ds * cx.c1 * hb;
        R(5,2) = ds * ds * cx.c2 * hb;
        R(5,6) = ds * ibetagam2 - ds * ds * ds * cx.c3 * hb * hb;

        return R;
    }

} // namespace impactx::elements

#endif // IMPACTX_ELEMENTS_HELPER_SECTORBEND_H
