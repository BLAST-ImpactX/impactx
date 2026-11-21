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
#include "particles/distribution/BoxMuller.H"
#include "particles/distribution/mixin/beamellipse.H"

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
            double muxpx=0.0,
            double muypy=0.0,
            double mutpt=0.0,
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
                      meanx, meany, meant, meanpx, meanpy, meanpt, dispx, disppx, dispy, disppy)
        {
            truncation_constants(cutx, m_outside_cutX, m_r2scaleX);
            truncation_constants(cuty, m_outside_cutY, m_r2scaleY);
            truncation_constants(cutt, m_outside_cutT, m_r2scaleT);
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
         * @param[out] outside_cut fraction of the untruncated distribution outside the cut;
         *                         0 means no truncation
         * @param[out] r2scale variance of the truncated distribution, relative to the untruncated one
         */
        static void truncation_constants (
            amrex::ParticleReal cut,
            amrex::ParticleReal & outside_cut,
            amrex::ParticleReal & r2scale
        )
        {
            using namespace amrex::literals;

            if (cut > 0.0_prt) {
                amrex::ParticleReal f = cut*cut/2.0_prt;
                amrex::ParticleReal ef = std::exp(-f);
                outside_cut = ef;
                r2scale = 1_prt - f*ef/(1_prt - ef);
            } else {
                outside_cut = 0_prt;
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
            // Generate three pairs of (optionally truncated) normal random variables using Box-Muller:
            box_muller::truncated_pair(x, px, m_outside_cutX, m_r2scaleX, engine);
            box_muller::truncated_pair(y, py, m_outside_cutY, m_r2scaleY, engine);
            box_muller::truncated_pair(t, pt, m_outside_cutT, m_r2scaleT, engine);

            // Transform to the beam ellipse:
            to_beam_ellipse(x, y, t, px, py, pt);
        }

    private:
        // constants that are independent of the individually sampled particle,
        // computed from the cut parameters in the constructor
        //! fraction of the untruncated distribution outside the cut; 0 means no truncation
        amrex::ParticleReal m_outside_cutX, m_outside_cutY, m_outside_cutT;
        amrex::ParticleReal m_r2scaleX, m_r2scaleY, m_r2scaleT;  //! variance of the truncated distribution, relative to the untruncated one
    };

} // namespace impactx::distribution

#endif // IMPACTX_DISTRIBUTION_GAUSSIAN
