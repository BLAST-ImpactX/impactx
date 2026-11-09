/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DISTRIBUTION_MIXIN_BEAMELLIPSE_H
#define IMPACTX_DISTRIBUTION_MIXIN_BEAMELLIPSE_H

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx::distribution::mixin
{
    /** The phase space ellipse of a beam: its second moments, dispersion and centroid.
     *
     * Distributions sample particles with zero mean and an identity covariance
     * matrix and then transform them to this ellipse.
     */
    struct BeamEllipse
    {
        /** A beam ellipse
         *
         * @param lambdax,lambday,lambdat phase space axis intercepts (position);
         *                                for zero correlation, these are the related RMS sizes (in meters)
         * @param lambdapx,lambdapy,lambdapt phase space axis intercepts (momentum)
         * @param muxpx,muypy,mutpt correlation length-momentum
         * @param meanx,meany,meant offsets of the mean (centroid) positions from those of the reference particle
         * @param meanpx,meanpy,meanpt offsets of the mean (centroid) momenta from those of the reference particle
         * @param dispx,disppx,dispy,disppy dispersion and its derivative in horizontal and vertical directions
         */
        BeamEllipse (
            amrex::ParticleReal lambdax,
            amrex::ParticleReal lambday,
            amrex::ParticleReal lambdat,
            amrex::ParticleReal lambdapx,
            amrex::ParticleReal lambdapy,
            amrex::ParticleReal lambdapt,
            amrex::ParticleReal muxpx,
            amrex::ParticleReal muypy,
            amrex::ParticleReal mutpt,
            amrex::ParticleReal meanx,
            amrex::ParticleReal meany,
            amrex::ParticleReal meant,
            amrex::ParticleReal meanpx,
            amrex::ParticleReal meanpy,
            amrex::ParticleReal meanpt,
            amrex::ParticleReal dispx,
            amrex::ParticleReal disppx,
            amrex::ParticleReal dispy,
            amrex::ParticleReal disppy
        )
        : m_lambdaX(lambdax), m_lambdaY(lambday), m_lambdaT(lambdat), m_lambdaPx(lambdapx), m_lambdaPy(lambdapy),
          m_lambdaPt(lambdapt), m_muxpx(muxpx), m_muypy(muypy), m_mutpt(mutpt), m_meanx(meanx), m_meany(meany),
          m_meant(meant), m_meanpx(meanpx), m_meanpy(meanpy), m_meanpt(meanpt), m_dispx(dispx), m_disppx(disppx),
          m_dispy(dispy), m_disppy(disppy)
        {
            using namespace amrex::literals;

            m_rootX = std::sqrt(1.0_prt-m_muxpx*m_muxpx);
            m_rootY = std::sqrt(1.0_prt-m_muypy*m_muypy);
            m_rootT = std::sqrt(1.0_prt-m_mutpt*m_mutpt);
        }

        /** Transform a particle from zero mean and an identity covariance matrix to the beam ellipse
         *
         * @param[inout] x particle position in x
         * @param[inout] y particle position in y
         * @param[inout] t particle position in t
         * @param[inout] px particle momentum in x
         * @param[inout] py particle momentum in y
         * @param[inout] pt particle momentum in t
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void to_beam_ellipse (
            amrex::ParticleReal & AMREX_RESTRICT x,
            amrex::ParticleReal & AMREX_RESTRICT y,
            amrex::ParticleReal & AMREX_RESTRICT t,
            amrex::ParticleReal & AMREX_RESTRICT px,
            amrex::ParticleReal & AMREX_RESTRICT py,
            amrex::ParticleReal & AMREX_RESTRICT pt
        ) const
        {
            amrex::ParticleReal a1, a2;

            // Transform to produce the desired second moments/correlations:
            a1 = m_lambdaX * x / m_rootX;
            a2 = m_lambdaPx * (-m_muxpx * x / m_rootX + px);
            x = a1;
            px = a2;
            a1 = m_lambdaY * y / m_rootY;
            a2 = m_lambdaPy * (-m_muypy * y / m_rootY + py);
            y = a1;
            py = a2;
            a1 = m_lambdaT * t / m_rootT;
            a2 = m_lambdaPt * (-m_mutpt * t / m_rootT + pt);
            t = a1;
            pt = a2;

            // Apply dispersive correlations
            x = x - m_dispx * pt;
            px = px - m_disppx * pt;
            y = y - m_dispy * pt;
            py = py - m_disppy * pt;

            // Apply centroid translation
            x = x + m_meanx;
            px = px + m_meanpx;
            y = y + m_meany;
            py = py + m_meanpy;
            t = t + m_meant;
            pt = pt + m_meanpt;
        }

        amrex::ParticleReal m_lambdaX, m_lambdaY, m_lambdaT;  //! related position axis intercepts (length) of the phase space ellipse
        amrex::ParticleReal m_lambdaPx, m_lambdaPy, m_lambdaPt;  //! related momentum axis intercepts of the phase space ellipse
        amrex::ParticleReal m_muxpx, m_muypy, m_mutpt;  //! correlation length-momentum
        amrex::ParticleReal m_meanx, m_meany, m_meant;  //! spatial coordinates of centroid offset
        amrex::ParticleReal m_meanpx, m_meanpy, m_meanpt;  //! momentum coordinates of centroid offset
        amrex::ParticleReal m_dispx, m_disppx, m_dispy, m_disppy;  //! dispersion and its derivative

    private:
        // constants that are independent of the individually sampled particle,
        // computed from the parameters above in the constructor
        amrex::ParticleReal m_rootX, m_rootY, m_rootT;  //! sqrt(1 - mu^2) of each phase plane
    };

} // namespace impactx::distribution::mixin

#endif // IMPACTX_DISTRIBUTION_MIXIN_BEAMELLIPSE_H
