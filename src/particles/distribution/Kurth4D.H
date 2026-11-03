/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DISTRIBUTION_KURTH4D
#define IMPACTX_DISTRIBUTION_KURTH4D

#include "particles/ReferenceParticle.H"
#include "particles/distribution/mixin/beamellipse.H"

#include <ablastr/constant.H>

#include <AMReX_Math.H>
#include <AMReX_Random.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx::distribution
{
    struct Kurth4D
    : public mixin::BeamEllipse
    {
        /** A 4D Kurth distribution transversely + a uniform distribution
         *  in t + a Gaussian distribution in pt
         *
         * Return sampling from a 4D Kurth + uniform distribution.
         *
         * @param lambdax,lambday,lambdat for zero correlation, these are the related
         *                                   RMS sizes (in meters)
         * @param lambdapx,lambdapy,lambdapt RMS momentum
         * @param muxpx,muypy,mutpt correlation length-momentum
         * @param meanx,meany,meant offsets of the mean (centroid) positions from those of the reference particle
         * @param meanpx,meanpy,meanpt offsets of the mean (centroid) momenta from those of the reference particle
         * @param dispx,disppx,dispy,disppy dispersion and its derivative in horizontal and vertical directions
         */
        Kurth4D (
            amrex::ParticleReal lambdax,
            amrex::ParticleReal lambday,
            amrex::ParticleReal lambdat,
            amrex::ParticleReal lambdapx,
            amrex::ParticleReal lambdapy,
            amrex::ParticleReal lambdapt,
            amrex::ParticleReal muxpx=0.0,
            amrex::ParticleReal muypy=0.0,
            amrex::ParticleReal mutpt=0.0,
            amrex::ParticleReal meanx=0.0,
            amrex::ParticleReal meany=0.0,
            amrex::ParticleReal meant=0.0,
            amrex::ParticleReal meanpx=0.0,
            amrex::ParticleReal meanpy=0.0,
            amrex::ParticleReal meanpt=0.0,
            amrex::ParticleReal dispx=0.0,
            amrex::ParticleReal disppx=0.0,
            amrex::ParticleReal dispy=0.0,
            amrex::ParticleReal disppy=0.0
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
            amrex::ParticleReal & AMREX_RESTRICT x,
            amrex::ParticleReal & AMREX_RESTRICT y,
            amrex::ParticleReal & AMREX_RESTRICT t,
            amrex::ParticleReal & AMREX_RESTRICT px,
            amrex::ParticleReal & AMREX_RESTRICT py,
            amrex::ParticleReal & AMREX_RESTRICT pt,
            amrex::RandomEngine const & engine
        ) const
        {
            using namespace amrex::literals;
            using amrex::Math::powi;
            using ablastr::constant::math::pi;

            amrex::ParticleReal v,phi,r,u1,u2,ln1;
            amrex::ParticleReal alpha,u,Lz,pmax,pr,pphi;

            // Sample and transform to define (x,y):
            v = amrex::Random(engine);
            phi = amrex::Random(engine);
            phi = 2_prt*pi*phi;
            r = std::sqrt(v);
            x = r * std::cos(phi);
            y = r * std::sin(phi);

            // Random samples used to define Lz:
            u = amrex::Random(engine);
            Lz = r*(2.0_prt*u-1.0_prt);

            // Random samples used to define pr:
            alpha = amrex::Random(engine);
            alpha = pi*alpha;
            pmax = 1.0_prt - powi<2>((Lz/r)) - powi<2>(r) + powi<2>(Lz);
            pmax = std::sqrt(pmax);
            pr = pmax * std::cos(alpha);
            pphi = Lz/r;

            // Transformations used to obtain (px,py):
            px = pr * std::cos(phi)-pphi * std::sin(phi);
            py = pr * std::sin(phi)+pphi * std::cos(phi);

            // Sample and transform to define (t,pt):
            t = amrex::Random(engine);
            t = 2.0_prt*(t-0.5_prt);
            u1 = amrex::Random(engine);
            u2 = amrex::Random(engine);
            ln1 = std::sqrt(-2_prt*std::log(u1));
            pt = ln1 * std::cos(2_prt*pi*u2);

            // Scale to produce the identity covariance matrix:
            amrex::ParticleReal const c = std::sqrt(3.0_prt);
            x = 2_prt*x;
            y = 2_prt*y;
            t = c*t;
            px = 2_prt*px;
            py = 2_prt*py;
            // pt = pt;

            // Transform to the beam ellipse:
            to_beam_ellipse(x, y, t, px, py, pt);
        }
    };

} // namespace impactx::distribution

#endif // IMPACTX_DISTRIBUTION_KURTH4D
