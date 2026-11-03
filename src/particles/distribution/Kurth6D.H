/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DISTRIBUTION_KURTH6D
#define IMPACTX_DISTRIBUTION_KURTH6D

#include "particles/ReferenceParticle.H"
#include "particles/distribution/mixin/beamellipse.H"

#include <ablastr/constant.H>

#include <AMReX_Math.H>
#include <AMReX_Random.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx::distribution
{
    struct Kurth6D
    : public mixin::BeamEllipse
    {
        /** A 6D Kurth distribution
         *
         * Return sampling from a 6D Kurth distribution.
         *
         * R. Kurth, Quarterly of Applied Mathematics vol. 32, pp. 325-329 (1978)
         * C. Mitchell, K. Hwang and R. D. Ryne, IPAC2021, WEPAB248 (2021)
         *
         * @param lambdax,lambday,lambdat for zero correlation, these are the related
         *                                   RMS sizes (in meters)
         * @param lambdapx,lambdapy,lambdapt RMS momentum
         * @param muxpx,muypy,mutpt correlation length-momentum
         * @param meanx,meany,meant offsets of the mean (centroid) positions from those of the reference particle
         * @param meanpx,meanpy,meanpt offsets of the mean (centroid) momenta from those of the reference particle
         * @param dispx,disppx,dispy,disppy dispersion and its derivative in horizontal and vertical directions
         */
        Kurth6D (
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

            amrex::ParticleReal v,costheta,sintheta,phi,r;
            amrex::ParticleReal L,alpha,pmax,pr,beta,p1,p2;

            // Random samples used to define (x,y,z):
            v = amrex::Random(engine);
            costheta = amrex::Random(engine);
            costheta = 2_prt*(costheta-0.5_prt);
            sintheta = std::sqrt(1_prt-powi<2>(costheta));
            phi = amrex::Random(engine);
            phi = 2_prt*pi*phi;

            // Transformations for (x,y,t):
            r = std::pow(v,1_prt/3_prt);
            x = r*sintheta * std::cos(phi);
            y = r*sintheta * std::sin(phi);
            t = r*costheta;

            // Random samples used to define L:
            L = amrex::Random(engine);
            L = r*std::sqrt(L);

            // Random samples used to define pr:
            alpha = amrex::Random(engine);
            alpha = pi*alpha;
            pmax = 1_prt - powi<2>(L/r) - powi<2>(r) + powi<2>(L);
            pmax = std::sqrt(pmax);
            pr = pmax * std::cos(alpha);

            // Random samples used to define ptangent:
            beta = amrex::Random(engine);
            beta = 2_prt*pi*beta;
            p1 = L/r * std::cos(beta);  // This is phi component
            p2 = L/r * std::sin(beta);  // This is theta component

            // Transformation from spherical to Cartesian coord.:
            px = pr*sintheta * std::cos(phi) + p2*costheta * std::cos(phi) - p1 * std::sin(phi);
            py = pr*sintheta * std::sin(phi) + p2*costheta * std::sin(phi) + p1 * std::cos(phi);
            pt = pr*costheta - p2*sintheta;

            // Scale to produce the identity covariance matrix:
            amrex::ParticleReal const c = std::sqrt(5.0_prt);
            x = c*x;
            y = c*y;
            t = c*t;
            px = c*px;
            py = c*py;
            pt = c*pt;

            // Transform to the beam ellipse:
            to_beam_ellipse(x, y, t, px, py, pt);
        }
    };

} // namespace impactx::distribution

#endif // IMPACTX_DISTRIBUTION_KURTH6D
