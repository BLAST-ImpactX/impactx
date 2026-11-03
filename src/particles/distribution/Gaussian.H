/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DISTRIBUTION_GAUSSIAN
#define IMPACTX_DISTRIBUTION_GAUSSIAN

#include "particles/ReferenceParticle.H"
#include "particles/distribution/mixin/beamellipse.H"

#include <ablastr/constant.H>

#include <AMReX_Random.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx::distribution
{
    struct Gaussian
    : public mixin::BeamEllipse
    {
        /** A Gaussian distribution
         *
         * Return sampling from a 6D Gaussian distribution, optionally with
         * truncation in each phase plane.
         *
         * @param lambdax,lambday,lambdat phase space axis intercepts (position);
         *                                for zero correlation, these are the related RMS sizes (in meters)
         * @param lambdapx,lambdapy,lambdapt phase space axis intercepts (momentum)
         * @param muxpx,muypy,mutpt correlation length-momentum
         * @param meanx,meany,meant offsets of the mean (centroid) positions from those of the reference particle
         * @param meanpx,meanpy,meanpt offsets of the mean (centroid) momenta from those of the reference particle
         * @param dispx,disppx,dispy,disppy dispersion and its derivative in horizontal and vertical directions
         * @param cutx,cuty,cutt dimensionless cut parameter (number of sigma) in each phase plane;
         *                       a value of 0 means no truncation (default)
         */
        Gaussian (
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
            amrex::ParticleReal disppy=0.0,
            amrex::ParticleReal cutx=0.0,
            amrex::ParticleReal cuty=0.0,
            amrex::ParticleReal cutt=0.0
        )
        : BeamEllipse(lambdax, lambday, lambdat, lambdapx, lambdapy, lambdapt, muxpx, muypy, mutpt,
                      meanx, meany, meant, meanpx, meanpy, meanpt, dispx, disppx, dispy, disppy),
          m_cutX(cutx), m_cutY(cuty), m_cutT(cutt)
        {
            truncation_constants(m_cutX, m_inside_cutX, m_r2scaleX);
            truncation_constants(m_cutY, m_inside_cutY, m_r2scaleY);
            truncation_constants(m_cutT, m_inside_cutT, m_r2scaleT);
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

        /** Constants of the radius in a phase plane, for a Gaussian truncated at a cut
         *
         * @param[in] cut dimensionless cut parameter (number of sigma); 0 means no truncation
         * @param[out] inside_cut fraction of the untruncated distribution inside the cut
         * @param[out] r2scale variance of the truncated distribution, relative to the untruncated one
         */
        static void truncation_constants (
            amrex::ParticleReal cut,
            amrex::ParticleReal & inside_cut,
            amrex::ParticleReal & r2scale
        )
        {
            using namespace amrex::literals;

            if (cut > 0.0_prt) {
                amrex::ParticleReal f = cut*cut/2.0_prt;
                amrex::ParticleReal ef = std::exp(-f);
                inside_cut = 1_prt - ef;
                r2scale = 1_prt - f*ef/(1_prt - ef);
            } else {
                // not used
                inside_cut = 1_prt;
                r2scale = 1_prt;
            }
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
            using ablastr::constant::math::pi;

            amrex::ParticleReal ln1,u1,u2;

            // Generate three pairs of (optionally truncated) normal random variables using Box-Muller:

            u1 = amrex::Random(engine);
            u2 = amrex::Random(engine);
            if (m_cutX > 0.0_prt) {
                amrex::ParticleReal r2 = -2_prt * std::log(1_prt - u1*m_inside_cutX);
                ln1 = std::sqrt(r2 / m_r2scaleX);
            } else {
                ln1 = std::sqrt(-2_prt*std::log(u1));
            }
            x = ln1 * std::cos(2_prt*pi*u2);
            px = ln1 * std::sin(2_prt*pi*u2);

            u1 = amrex::Random(engine);
            u2 = amrex::Random(engine);
            if (m_cutY > 0.0_prt) {
                amrex::ParticleReal r2 = -2_prt * std::log(1_prt - u1*m_inside_cutY);
                ln1 = std::sqrt(r2 / m_r2scaleY);
            } else {
                ln1 = std::sqrt(-2_prt*std::log(u1));
            }
            y = ln1 * std::cos(2_prt*pi*u2);
            py = ln1 * std::sin(2_prt*pi*u2);

            u1 = amrex::Random(engine);
            u2 = amrex::Random(engine);
            if (m_cutT > 0.0_prt) {
                amrex::ParticleReal r2 = -2_prt * std::log(1_prt - u1*m_inside_cutT);
                ln1 = std::sqrt(r2 / m_r2scaleT);
            } else {
                ln1 = std::sqrt(-2_prt*std::log(u1));
            }
            t = ln1 * std::cos(2_prt*pi*u2);
            pt = ln1 * std::sin(2_prt*pi*u2);

            // Transform to the beam ellipse:
            to_beam_ellipse(x, y, t, px, py, pt);
        }

        amrex::ParticleReal m_cutX, m_cutY, m_cutT;  //! dimensionless cutoff for each phase plane; 0 means no cut

    private:
        // constants that are independent of the individually sampled particle,
        // computed from the parameters above in the constructor
        amrex::ParticleReal m_inside_cutX, m_inside_cutY, m_inside_cutT;  //! fraction of the untruncated distribution inside the cut
        amrex::ParticleReal m_r2scaleX, m_r2scaleY, m_r2scaleT;  //! variance of the truncated distribution, relative to the untruncated one
    };

} // namespace impactx::distribution

#endif // IMPACTX_DISTRIBUTION_GAUSSIAN
