/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Ji Qiang, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DISTRIBUTION_TRIANGLE
#define IMPACTX_DISTRIBUTION_TRIANGLE

#include "particles/ReferenceParticle.H"
#include "particles/distribution/mixin/beamellipse.H"

#include <ablastr/constant.H>

#include <AMReX_Random.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx::distribution
{
    struct Triangle
    : public mixin::BeamEllipse
    {
        /** A Triangle distribution for LPA applications.
         *
         * Return sampling from a ramped, triangular current profile
         * with a Gaussian energy spread (possibly correlated).  The
         * transverse distribution is a 4D waterbag.
         *
         * @param lambdax,lambday,lambdat for zero correlation, these are the related
         *                       RMS sizes (in meters)
         * @param lambdapx,lambdapy,lambdapt RMS momentum
         * @param muxpx,muypy,mutpt correlation length-momentum
         * @param meanx,meany,meant offsets of the mean (centroid) positions from those of the reference particle
         * @param meanpx,meanpy,meanpt offsets of the mean (centroid) momenta from those of the reference particle
         * @param dispx,disppx,dispy,disppy dispersion and its derivative in horizontal and vertical directions
         */
        Triangle (
                amrex::ParticleReal const lambdax, amrex::ParticleReal const lambday,
                amrex::ParticleReal const lambdat, amrex::ParticleReal const lambdapx,
                amrex::ParticleReal const lambdapy, amrex::ParticleReal const lambdapt,
                amrex::ParticleReal const muxpx=0.0, amrex::ParticleReal const muypy=0.0,
                amrex::ParticleReal const mutpt=0.0, amrex::ParticleReal const meanx=0.0,
                amrex::ParticleReal const meany=0.0, amrex::ParticleReal const meant=0.0,
                amrex::ParticleReal const meanpx=0.0, amrex::ParticleReal const meanpy=0.0,
                amrex::ParticleReal const meanpt=0.0, amrex::ParticleReal const dispx=0.0,
                amrex::ParticleReal const disppx=0.0, amrex::ParticleReal const dispy=0.0,
                amrex::ParticleReal const disppy=0.0
        )
          : BeamEllipse(lambdax, lambday, lambdat, lambdapx, lambdapy, lambdapt, muxpx, muypy, mutpt,
                        meanx, meany, meant, meanpx, meanpy, meanpt, dispx, disppx, dispy, disppy)
        {
        }

        /** Initialize the distribution.
         *
         * Nothing to do here.
         *
         * @param bunch_charge charge of the beam in C
         * @param ref the reference particle
         */
        void initialize ([[maybe_unused]] amrex::ParticleReal bunch_charge, [[maybe_unused]] RefPart const & ref)
        {
        }

        /** Close and deallocate all data and handles.
         *
         * Nothing to do here.
         */
        void
        finalize ()
        {
        }

        /** Return 1 6D particle coordinate
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param engine a random number engine (with associated state)
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            amrex::ParticleReal & x,
            amrex::ParticleReal & y,
            amrex::ParticleReal & t,
            amrex::ParticleReal & px,
            amrex::ParticleReal & py,
            amrex::ParticleReal & pt,
            amrex::RandomEngine const& engine
        ) const
        {

            using namespace amrex::literals;
            using ablastr::constant::math::pi;

            amrex::ParticleReal ln1, norm, u0, u1, u2;
            amrex::ParticleReal g1, g2, g3, g4, g5;
            amrex::ParticleReal d;

            // Sample the t coordinate for a ramped triangular profile (unit
            // variance):
            u0 = amrex::Random(engine);
            t = std::sqrt(2_prt)*(2_prt-3_prt*std::sqrt(u0));

            // Generate five standard normal random variables using Box-Muller:
            u1 = amrex::Random(engine);
            u2 = amrex::Random(engine);
            ln1 = std::sqrt(-2_prt*std::log(u1));
            g1 = ln1 * std::cos(2_prt*pi*u2);
            g2 = ln1 * std::sin(2_prt*pi*u2);
            u1 = amrex::Random(engine);
            u2 = amrex::Random(engine);
            ln1 = std::sqrt(-2_prt*std::log(u1));
            g3 = ln1 * std::cos(2_prt*pi*u2);
            g4 = ln1 * std::sin(2_prt*pi*u2);
            u1 = amrex::Random(engine);
            u2 = amrex::Random(engine);
            ln1 = std::sqrt(-2_prt*std::log(u1));
            g5 = ln1 * std::cos(2_prt*pi*u2);

            // Use one of these normal random variables for pt:
            pt = g5;

            // Normalize the rest to produce uniform samples on the unit sphere:
            norm = std::sqrt(g1*g1+g2*g2+g3*g3+g4*g4);
            g1 /= norm;
            g2 /= norm;
            g3 /= norm;
            g4 /= norm;

            // Scale to produce uniform samples in a 4D ball (unit variance):
            d = 4_prt;  // unit ball dimension
            u1 = amrex::Random(engine);   // uniform sample
            u2 = std::sqrt(d+2_prt) * std::pow(u1,1_prt/d);
            x = g1*u2;
            y = g2*u2;
            px = g3*u2;
            py = g4*u2;

            // Transform to the beam ellipse:
            to_beam_ellipse(x, y, t, px, py, pt);
        }
    };

} // namespace impactx::distribution

#endif // IMPACTX_DISTRIBUTION_TRIANGLE
