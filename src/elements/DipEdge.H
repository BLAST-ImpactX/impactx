/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DIPEDGE_H
#define IMPACTX_DIPEDGE_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/alignment.H"
#include "mixin/beamoptic.H"
#include "mixin/thin.H"
#include "mixin/lineartransport.H"
#include "mixin/named.H"
#include "mixin/nofinalize.H"

#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_SIMD.H>

#include <cmath>
#include <stdexcept>


namespace impactx::elements
{
    struct DipEdge
    : public mixin::Named,
      public mixin::BeamOptic<DipEdge>,
      public mixin::LinearTransport<DipEdge>,
      public mixin::Thin,
      public mixin::Alignment,
      public mixin::NoFinalize
      // At least on Intel AVX512, there is a small overhead to vectorize this element, see
      // https://github.com/BLAST-ImpactX/impactx/pull/1002
      // public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "DipEdge";
        using PType = ImpactXParticleContainer::ParticleType;

        /** Edge focusing associated with bend entry or exit
         *
         * This model assumes a first-order effect of nonzero gap.
         * Here we use the linear fringe field map, given to first order in g/rc (gap / radius of curvature).
         *
         * References:
         *   K. L. Brown, SLAC Report No. 75 (1982).
         *   K. Hwang and S. Y. Lee, PRAB 18, 122401 (2015).
         *
         * @param psi Pole face entry/exit angle in rad.
         * @param rc Radius of curvature in m.
         * @param g Gap parameter in m.
         * @param R Length scale parameter in m.
         * @param K0 Fringe field integeral (unitless).
         * @param K1 Fringe field integeral (unitless).
         * @param K2 Fringe field integral (unitless).
         * @param K3 Fringe field integral (unitless).
         * @param K4 Fringe field integral (unitless).
         * @param K5 Fringe field integral (unitless).
         * @param K6 Fringe field integral (unitless).
         * @param model Model (0 linear, 1 nonlinear).
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param name a user defined and not necessarily unique name of the element
         */
        DipEdge (
            amrex::ParticleReal psi,
            amrex::ParticleReal rc,
            amrex::ParticleReal g,
            amrex::ParticleReal R,
            amrex::ParticleReal K0,
            amrex::ParticleReal K1,
            amrex::ParticleReal K2,
            amrex::ParticleReal K3,
            amrex::ParticleReal K4,
            amrex::ParticleReal K5,
            amrex::ParticleReal K6,
            int model,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0,
            std::optional<std::string> name = std::nullopt
        )
        : Named(std::move(name)),
          Alignment(dx, dy, rotation_degree),
          m_psi(psi), m_rc(rc), m_g(g), m_R(R), m_K0(K0), m_K1(K1), m_K2(K2), m_K3(K3), m_K4(K4), m_K5(K5), m_K6(K6), m_model(model)
        {
        }

        /** Push all particles */
        using BeamOptic::operator();

        /** Compute and cache the constants for the push.
         *
         * In particular, used to pre-compute and cache variables that are
         * independent of the individually tracked particle.
         *
         * @param refpart reference particle (unused)
         */
        void compute_constants ([[maybe_unused]] RefPart const & refpart)
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;

            Alignment::compute_constants(refpart);

            // edge focusing matrix elements (zero gap)
            m_R21 = std::tan(m_psi) / m_rc;

            // first-order effect of nonzero gap
            auto const [sin_psi, cos_psi] = amrex::Math::sincos(m_psi);
            amrex::ParticleReal const vf = (1.0_prt + powi<2>(sin_psi))
                / powi<3>(cos_psi)
                * m_g * m_K2 / powi<2>(m_rc);

            m_R43 = vf - m_R21;
        }

        /** This is a dipedge functor, so that a variable of this type can be used like a
         *  dipedge function.
         *
         * The @see compute_constants method must be called before pushing particles through this operator.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t (unchanged)
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t (unchanged)
         * @param idcpu particle global index (unused)
         * @param refpart reference particle (unused)
         */
        template<typename T_Real=amrex::ParticleReal, typename T_IdCpu=uint64_t>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            T_Real & AMREX_RESTRICT x,
            T_Real & AMREX_RESTRICT y,
            [[maybe_unused]] T_Real & AMREX_RESTRICT t,
            T_Real & AMREX_RESTRICT px,
            T_Real & AMREX_RESTRICT py,
            [[maybe_unused]] T_Real & AMREX_RESTRICT pt,
            [[maybe_unused]] T_IdCpu & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
        ) const
        {

            using namespace amrex::literals; // for _rt and _prt

            // shift due to alignment errors of the element
            shift_in(x, y, px, py);

            // apply edge focusing
            px = px + m_R21 * x;
            py = py + m_R43 * y;

            // undo shift due to alignment errors of the element
            shift_out(x, y, px, py);
        }

        /** This pushes the reference particle. */
        using Thin::operator();

        /** This pushes the covariance matrix. */
        using LinearTransport::operator();

        /** This function returns the linear transport map.
         *
         * @param[in] refpart reference particle
         * @returns 6x6 transport matrix
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        Map6x6
        transport_map ([[maybe_unused]] RefPart const & AMREX_RESTRICT refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;

            // initialize linear map matrix elements
            Map6x6 R = Map6x6::Identity();

            // edge focusing matrix elements (zero gap)
            amrex::ParticleReal const R21 = std::tan(m_psi) / m_rc;

            // first-order effect of nonzero gap
            auto const [sin_psi, cos_psi] = amrex::Math::sincos(m_psi);
            amrex::ParticleReal const vf = (1.0_prt + powi<2>(sin_psi))
                / powi<3>(cos_psi)
                * m_g * m_K2 / powi<2>(m_rc);

            amrex::ParticleReal const R43 = vf - R21;

            // set non-identity matrix elements
            R(2,1) = R21;
            R(4,3) = R43;

            return R;
        }

        amrex::ParticleReal m_psi; //! pole face angle in rad
        amrex::ParticleReal m_rc; //! bend radius in m
        amrex::ParticleReal m_g; //! gap parameter in m
        amrex::ParticleReal m_R; //! scale length in m
        amrex::ParticleReal m_K0; //! fringe field integral
        amrex::ParticleReal m_K1; //! fringe field integral
        amrex::ParticleReal m_K2; //! fringe field integral
        amrex::ParticleReal m_K3; //! fringe field integral
        amrex::ParticleReal m_K4; //! fringe field integral
        amrex::ParticleReal m_K5; //! fringe field integral
        amrex::ParticleReal m_K6; //! fringe field integral
        int m_model; //! model specification


      private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_R21, m_R43;
    };

} // namespace impactx

#endif // IMPACTX_DIPEDGE_H
