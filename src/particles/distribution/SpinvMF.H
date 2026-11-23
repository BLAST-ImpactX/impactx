/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DISTRIBUTION_SPINVMF
#define IMPACTX_DISTRIBUTION_SPINVMF

#include <ablastr/constant.H>

#include <AMReX_Math.H>
#include <AMReX_Random.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx::distribution
{
    struct SpinvMF
    {
        /** A von Mises-Fisher (vMF) distribution on the unit 2-sphere.
         *
         * This is used for initializing particle spin.  There is a natural
         * bijective correspondence between vMF distributions and mean
         * (polarization) vectors.
         *
         * The algorithm used here is a simplification of the algorithm described in:
         * C. Pinzon and K. Jung, "Fast Python sampler of the von Mises Fisher distribution", in the special
         * case of the 2-sphere.  Additional references used include:
         *
         * K. V. Mardia and P. E. Jupp, Directional Statistics, Wiley, 1999;
         * S. Kang and H-S. Oh, "Novel sampling method for the von Mises-Fisher distribution", Stat.
         * and Comput. 34, 106 (2024), https://doi.org/10.1007/s11222-024-10419-3.
         *
         * Return sampling from a vMF distribution.
         *
         * @param mux,muy,muz components of the polarization vector: its direction is the mean
         *                    spin direction and its magnitude, less than 1, is the degree of
         *                    polarization. A zero vector describes an unpolarized beam.
         */
        SpinvMF (
            amrex::ParticleReal mux,
            amrex::ParticleReal muy,
            amrex::ParticleReal muz
        )
        : m_muX(mux), m_muY(muy), m_muZ(muz)
        {
            using namespace amrex::literals;  // for _rt and _prt

            // magnitude of the polarization vector (= polarization magnitude)
            amrex::ParticleReal const pmag = std::sqrt(
                mux * mux +
                muy * muy +
                muz * muz
            );

            m_kappa = inverse_Langevin(pmag);

            // Mean direction: without polarization, the distribution is uniform
            // on the sphere (kappa = 0), so any direction is correct.
            if (pmag == 0_prt) {
                m_dirX = 0_prt;
                m_dirY = 0_prt;
                m_dirZ = 1_prt;
            } else {
                m_dirX = mux / pmag;
                m_dirY = muy / pmag;
                m_dirZ = muz / pmag;
            }

            // Transverse part of the mean direction, for the local frame of the sampling:
            m_dirPerp = std::sqrt(m_dirX * m_dirX + m_dirY * m_dirY);
            m_a = (m_dirPerp == 0_prt) ? 0_prt : m_dirX / m_dirPerp;
            m_b = (m_dirPerp == 0_prt) ? 1_prt : m_dirY / m_dirPerp;
        }

        /** Return 1 3D spin (unit) vector
         *
         * @param sx particle spin component in x
         * @param sy particle spin component in y
         * @param sz particle spin component in z
         * @param engine a random number engine (with associated state)
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            amrex::ParticleReal & AMREX_RESTRICT sx,
            amrex::ParticleReal & AMREX_RESTRICT sy,
            amrex::ParticleReal & AMREX_RESTRICT sz,
            amrex::RandomEngine const & engine
        ) const
        {
            using namespace amrex::literals;
            using amrex::Math::powi;
            using ablastr::constant::math::pi;

            // Sample two iid random variables:
            amrex::ParticleReal x1 = amrex::Random(engine);
            amrex::ParticleReal x2 = amrex::Random(engine);

            // Basic transformation of the random variables:
            amrex::ParticleReal c = std::cos(2_prt*pi*x1);
            amrex::ParticleReal s = std::sin(2_prt*pi*x1);
            amrex::ParticleReal t;
            // The form of the function below is obtained from eq (30) of:
            // D. Frisch and U. D. Hanebeck, "Deterministic Von Mises-Fisher Sampling on the Sphere Using
            // Fibonacci Lattices", IEEE, 2023, DOI: 10.1109/SDF-MFI59545.2023.10361396
            if(m_kappa > 0) {
               t = 1_prt + std::log1p(x2 * std::expm1(-2_prt * m_kappa)) / m_kappa;
            } else {
               t = 1_prt - 2_prt * x2;
            }

            // Evaluation of the spin components
            amrex::ParticleReal u = std::sqrt(1_prt-powi<2>(t));
            sx = u * (m_b*c + m_a*m_dirZ*s) + t*m_dirX;
            sy = u * (-m_a*c + m_b*m_dirZ*s) + t*m_dirY;
            sz = u * (-m_dirPerp*s) + t*m_dirZ;

        }


        /** This function evaluates the inverse Langevin function, in order to return
         *  the value of concentration (kappa) required to produce a given polarization magnitude.
         *  It is used to prepare for spin sampling using the vMF distribution.
         *
         *  Here, we use an approximate expression due to R. Petrosyan.  See:
         *  R. Petrosyan, Rheologica Acta 56, 21-26, 2016; doi:10.1007/s00397-016-0977-9
         *  R. Jedynak, Mathematics and Mechanics of Solids 24, 1-25, 2018; doi:10.1177/1081286518811395
         *
         * @param pmag polarization magnitude
         * @returns required value of concentration kappa
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        static amrex::ParticleReal
        inverse_Langevin (amrex::ParticleReal pmag)
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;

            AMREX_ASSERT_WITH_MESSAGE(pmag >= 0_prt, "Polarization magnitude must be greater or equal 0.");
            AMREX_ASSERT_WITH_MESSAGE(pmag < 1_prt, "Polarization magnitude must be less than 1.");

            // Evaluation of the inverse Langevin function.  Here, we use the approximate expression
            // due to R. Petrosyan.  See:
            // R. Petrosyan, Rheologica Acta 56, 21-26, 2016; doi:10.1007/s00397-016-0977-9
            // R. Jedynak, Mathematics and Mechanics of Solids 24, 1-25, 2018; doi:10.1177/1081286518811395

            amrex::ParticleReal concentration =
                3_prt * pmag +
                powi<2>(pmag) / 5_prt * std::sin(7_prt * pmag * 0.5_prt) +
                powi<3>(pmag) / (1_prt - pmag);

            return concentration;
        }

        amrex::ParticleReal m_muX, m_muY, m_muZ;  //! components of the polarization vector
        amrex::ParticleReal m_kappa;  //! concentration parameter

    private:
        amrex::ParticleReal m_dirX, m_dirY, m_dirZ;  //! unit vector of the mean spin direction
        amrex::ParticleReal m_dirPerp;  //! magnitude of the transverse (x,y) part of the mean direction
        amrex::ParticleReal m_a, m_b;  //! transverse unit direction (x,y) of the mean direction
    };

} // namespace impactx::distribution

#endif // IMPACTX_DISTRIBUTION_SPINVMF
