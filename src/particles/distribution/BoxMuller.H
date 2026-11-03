/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DISTRIBUTION_BOXMULLER_H
#define IMPACTX_DISTRIBUTION_BOXMULLER_H

#include <ablastr/constant.H>

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Random.H>
#include <AMReX_REAL.H>

#include <cmath>


/** Normal random variables from the Box-Muller transform
 *
 * Two uniform draws give a radius and an angle. The draw for the radius, u1, is in (0,1], so that
 * its logarithm is finite.
 */
namespace impactx::distribution::box_muller
{
    /** Sample a pair of independent standard normal random variables
     *
     * @param[out] g1,g2 the normal random variables
     * @param[in] engine a random number engine (with associated state)
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void pair (
        amrex::ParticleReal & AMREX_RESTRICT g1,
        amrex::ParticleReal & AMREX_RESTRICT g2,
        amrex::RandomEngine const & engine
    )
    {
        using namespace amrex::literals;
        using ablastr::constant::math::pi;

        amrex::ParticleReal const u1 = amrex::RandomPositive(engine);
        amrex::ParticleReal const u2 = amrex::Random(engine);
        amrex::ParticleReal const ln1 = std::sqrt(-2_prt*std::log(u1));
        g1 = ln1 * std::cos(2_prt*pi*u2);
        g2 = ln1 * std::sin(2_prt*pi*u2);
    }

    /** Sample one standard normal random variable
     *
     * This is one variable of a pair, for when the other one is not needed.
     *
     * @param[in] engine a random number engine (with associated state)
     * @return the normal random variable
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::ParticleReal single (amrex::RandomEngine const & engine)
    {
        using namespace amrex::literals;
        using ablastr::constant::math::pi;

        amrex::ParticleReal const u1 = amrex::RandomPositive(engine);
        amrex::ParticleReal const u2 = amrex::Random(engine);
        amrex::ParticleReal const ln1 = std::sqrt(-2_prt*std::log(u1));
        return ln1 * std::cos(2_prt*pi*u2);
    }

    /** Sample a pair of independent normal random variables, truncated in their radius
     *
     * The radius is r^2 = -2 ln(w), with w = u1 + (1 - u1)*outside_cut uniform in (outside_cut, 1].
     * w >= u1 > 0 holds also if outside_cut rounds to 0.
     *
     * @param[out] g1,g2 the normal random variables
     * @param[in] outside_cut fraction of the untruncated distribution outside the cut;
     *                        0 means no truncation
     * @param[in] r2scale variance of the truncated distribution, relative to the untruncated one,
     *                    to scale the variables to unit variance
     * @param[in] engine a random number engine (with associated state)
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void truncated_pair (
        amrex::ParticleReal & AMREX_RESTRICT g1,
        amrex::ParticleReal & AMREX_RESTRICT g2,
        amrex::ParticleReal outside_cut,
        amrex::ParticleReal r2scale,
        amrex::RandomEngine const & engine
    )
    {
        using namespace amrex::literals;
        using ablastr::constant::math::pi;

        amrex::ParticleReal const u1 = amrex::RandomPositive(engine);
        amrex::ParticleReal const u2 = amrex::Random(engine);
        amrex::ParticleReal ln1;
        if (outside_cut > 0.0_prt) {
            amrex::ParticleReal const r2 = -2_prt * std::log(u1 + (1_prt - u1)*outside_cut);
            ln1 = std::sqrt(r2 / r2scale);
        } else {
            ln1 = std::sqrt(-2_prt*std::log(u1));
        }
        g1 = ln1 * std::cos(2_prt*pi*u2);
        g2 = ln1 * std::sin(2_prt*pi*u2);
    }

} // namespace impactx::distribution::box_muller

#endif // IMPACTX_DISTRIBUTION_BOXMULLER_H
