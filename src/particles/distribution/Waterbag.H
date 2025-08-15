/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Ji Qiang, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DISTRIBUTION_WATERBAG
#define IMPACTX_DISTRIBUTION_WATERBAG

#include "particles/ReferenceParticle.H"

#include <ablastr/constant.H>

#include <AMReX_Random.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx::distribution
{
    struct Waterbag
    {
        /** A Waterbag distribution
         *
         * Return sampling from a 6D Waterbag distribution.
         *
         * @param lambda_x,lambda_y,lambda_t for zero correlation, these are the related
         *                                   RMS sizes (in meters)
         * @param lambda_px,lambda_py,lambda_pt RMS momentum
         * @param mu_x_px,mu_y_py,mu_t_pt correlation length-momentum
         * @param mean_x,mean_y,mean_t offsets of the mean (centroid) positions from those of the reference particle
         * @param mean_px,mean_py,mean_pt offsets of the mean (centroid) momenta from those of the reference particle
         * @param dispersion_x,dispersion_px,dispersion_y,dispersion_py dispersion and its derivative in horizontal and vertical directions
         */
        Waterbag (
            amrex::ParticleReal lambda_x,
            amrex::ParticleReal lambda_y,
            amrex::ParticleReal lambda_t,
            amrex::ParticleReal lambda_px,
            amrex::ParticleReal lambda_py,
            amrex::ParticleReal lambda_pt,
            amrex::ParticleReal mu_x_px=0.0,
            amrex::ParticleReal mu_y_py=0.0,
            amrex::ParticleReal mu_t_pt=0.0,
            amrex::ParticleReal mean_x=0.0,
            amrex::ParticleReal mean_y=0.0,
            amrex::ParticleReal mean_t=0.0,
            amrex::ParticleReal mean_px=0.0,
            amrex::ParticleReal mean_py=0.0,
            amrex::ParticleReal mean_pt=0.0,
            amrex::ParticleReal dispersion_x=0.0,
            amrex::ParticleReal dispersion_px=0.0,
            amrex::ParticleReal dispersion_y=0.0,
            amrex::ParticleReal dispersion_py=0.0
        )
        : m_lambda_x(lambda_x), m_lambda_y(lambda_y), m_lambda_t(lambda_t), m_lambda_px(lambda_px), m_lambda_py(lambda_py),
          m_lambda_pt(lambda_pt), m_mu_x_px(mu_x_px), m_mu_y_py(mu_y_py), m_mu_t_pt(mu_t_pt), m_mean_x(mean_x), m_mean_y(mean_y),
          m_mean_t(mean_t), m_mean_px(mean_px), m_mean_py(mean_py), m_mean_pt(mean_pt), m_dispersion_x(dispersion_x), m_dispersion_px(dispersion_px),
          m_dispersion_y(dispersion_y), m_dispersion_py(dispersion_py)
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
            using ablastr::constant::math::pi;

            amrex::ParticleReal ln1,norm,u1,u2;
            amrex::ParticleReal g1,g2,g3,g4,g5,g6;
            amrex::ParticleReal root,a1,a2;

            // Generate six standard normal random variables using Box-Muller:
            u1 = amrex::Random(engine);
            u2 = amrex::Random(engine);
            ln1 = std::sqrt(-2_prt*std::log(u1));
            g1 = ln1*std::cos(2_prt*pi*u2);
            g2 = ln1*std::sin(2_prt*pi*u2);
            u1 = amrex::Random(engine);
            u2 = amrex::Random(engine);
            ln1 = std::sqrt(-2_prt*std::log(u1));
            g3 = ln1*std::cos(2_prt*pi*u2);
            g4 = ln1*std::sin(2_prt*pi*u2);
            u1 = amrex::Random(engine);
            u2 = amrex::Random(engine);
            ln1 = std::sqrt(-2_prt*std::log(u1));
            g5 = ln1*std::cos(2_prt*pi*u2);
            g6 = ln1*std::sin(2_prt*pi*u2);

            // Normalize to produce uniform samples on the unit sphere:
            norm = std::sqrt(g1*g1+g2*g2+g3*g3+g4*g4+g5*g5+g6*g6);
            g1 /= norm;
            g2 /= norm;
            g3 /= norm;
            g4 /= norm;
            g5 /= norm;
            g6 /= norm;

            // Scale to produce uniform samples in a ball (unit variance):
            u1 = amrex::Random(engine);
            u2 = std::sqrt(8_prt)*std::pow(u1,1_prt/6_prt);
            x = g1*u2;
            y = g2*u2;
            t = g3*u2;
            px = g4*u2;
            py = g5*u2;
            pt = g6*u2;

            // Transform to produce the desired second moments/correlations:
            root = std::sqrt(1.0_prt-m_mu_x_px*m_mu_x_px);
            a1 = m_lambda_x * x / root;
            a2 = m_lambda_px * (-m_mu_x_px * x / root + px);
            x = a1;
            px = a2;
            root = std::sqrt(1.0_prt-m_mu_y_py*m_mu_y_py);
            a1 = m_lambda_y * y / root;
            a2 = m_lambda_py * (-m_mu_y_py * y / root + py);
            y = a1;
            py = a2;
            root = std::sqrt(1.0_prt-m_mu_t_pt*m_mu_t_pt);
            a1 = m_lambda_t * t / root;
            a2 = m_lambda_pt * (-m_mu_t_pt * t / root + pt);
            t = a1;
            pt = a2;

            // Apply dispersive correlations
            x = x - m_dispersion_x * pt;
            px = px - m_dispersion_px * pt;
            y = y - m_dispersion_y * pt;
            py = py - m_dispersion_py * pt;

            // Apply centroid translation
            x = x + m_mean_x;
            px = px + m_mean_px;
            y = y + m_mean_y;
            py = py + m_mean_py;
            t = t + m_mean_t;
            pt = pt + m_mean_pt;
        }

        amrex::ParticleReal m_lambda_x, m_lambda_y, m_lambda_t;  //! related position axis intercepts (length) of the phase space ellipse
        amrex::ParticleReal m_lambda_px, m_lambda_py, m_lambda_pt;  //! related momentum axis intercepts of the phase space ellipse
        amrex::ParticleReal m_mu_x_px, m_mu_y_py, m_mu_t_pt;  //! correlation length-momentum
        amrex::ParticleReal m_mean_x, m_mean_y, m_mean_t;  //! spatial coordinates of centroid offset
        amrex::ParticleReal m_mean_px, m_mean_py, m_mean_pt;  //! momentum coordinates of centroid offset
        amrex::ParticleReal m_dispersion_x, m_dispersion_px, m_dispersion_y, m_dispersion_py;  //! dispersion and its derivative
    };

} // namespace impactx::distribution

#endif // IMPACTX_DISTRIBUTION_WATERBAG
