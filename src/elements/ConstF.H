/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_CONSTF_H
#define IMPACTX_CONSTF_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/alignment.H"
#include "mixin/pipeaperture.H"
#include "mixin/beamoptic.H"
#include "mixin/thick.H"
#include "mixin/lineartransport.H"
#include "mixin/named.H"
#include "mixin/nofinalize.H"

#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <stdexcept>

namespace impactx::elements
{
    struct ConstF
    : public mixin::Named,
      public mixin::BeamOptic<ConstF>,
      public mixin::LinearTransport<ConstF>,
      public mixin::Thick,
      public mixin::Alignment,
      public mixin::PipeAperture,
      public mixin::NoFinalize
    {
        static constexpr auto type = "ConstF";
        using PType = ImpactXParticleContainer::ParticleType;

        /** A linear Constant Focusing element
         *
         * @param ds Segment length in m.
         * @param kx Focusing strength for x in 1/m.
         * @param ky Focusing strength for y in 1/m.
         * @param kt Focusing strength for t in 1/m.
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         * @param nslice number of slices used for the application of space charge
         * @param name a user defined and not necessarily unique name of the element
         */
        ConstF (
            amrex::ParticleReal ds,
            amrex::ParticleReal kx,
            amrex::ParticleReal ky,
            amrex::ParticleReal kt,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0,
            amrex::ParticleReal aperture_x = 0,
            amrex::ParticleReal aperture_y = 0,
            int nslice = 1,
            std::optional<std::string> name = std::nullopt
        )
        : Named(std::move(name)),
          Thick(ds, nslice),
          Alignment(dx, dy, rotation_degree),
          PipeAperture(aperture_x, aperture_y),
          m_kx(kx), m_ky(ky), m_kt(kt)
        {
        }

        /** Push all particles */
        using BeamOptic::operator();

        /** Compute and cache the constants for the push.
         *
         * In particular, used to pre-compute and cache variables that are
         * independent of the individually tracked particle.
         *
         * @param refpart reference particle
         */
        void compute_constants (RefPart const & refpart)
        {
            using namespace amrex::literals; // for _rt and _prt

            Alignment::compute_constants(refpart);

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // find beta*gamma^2
            amrex::ParticleReal const betgam2 = amrex::Math::powi<2>(refpart.pt) - 1.0_prt;

            // trigo
	    // Note: The convention here is m_kx[1/m] = sign(k)*sqrt(abs(k))
	    // where k[1/m^2] is in the usual MADX convention
	    if (m_kx > 0)
	    {
	      std::cout << "kx+" << std::endl;
	      auto const [sin_kxds, cos_kxds] = amrex::Math::sincos(m_kx * slice_ds);
	      m_cos_kxds = cos_kxds;
	      m_const_x  = -m_kx * sin_kxds;
	      m_sincx    = sin_kxds / m_kx;
	    }
	    else if (m_kx < 0)
	    {
	      auto const sinh_kxds = std::sinh(-m_kx*slice_ds);
	      auto const cosh_kxds = std::cosh(-m_kx*slice_ds);
	      m_cos_kxds = cosh_kxds;
	      m_const_x  = -m_kx * sinh_kxds;
	      m_sincx    = -sinh_kxds / m_kx;
	    }
	    else //m_kx == 0
	    {
	      m_cos_kxds = 1.0;
	      m_const_x  = 0.0;
	      m_sincx    = slice_ds;
	    }

	    if (m_ky > 0)
	    {
	      auto const [sin_kyds, cos_kyds] = amrex::Math::sincos(m_ky * slice_ds);
	      m_cos_kyds = cos_kyds;
	      m_const_y  = -m_ky * sin_kyds;
	      m_sincy    = sin_kyds / m_ky;
	    }
	    else if (m_ky < 0)
	    {
	      auto const sinh_kyds = std::sinh(-m_ky*slice_ds);
	      auto const cosh_kyds = std::cosh(-m_ky*slice_ds);
	      m_cos_kyds = cosh_kyds;
	      m_const_y  = -m_ky * sinh_kyds;
	      m_sincy    = -sinh_kyds / m_ky;
	    }
	    else //m_ky == 0
	    {
	      m_cos_kyds = 1.0;
	      m_const_y  = 0.0;
	      m_sincy    = slice_ds;
	    }

	    if (m_kt < 0)
	    {
	      throw std::runtime_error(std::string(type) + ": must have kt >= 0!");
	    }
            auto const [sin_ktds, cos_ktds] = amrex::Math::sincos(m_kt * slice_ds);
            m_cos_ktds = cos_ktds;
            m_const_t = -m_kt * betgam2 * sin_ktds;
            // intermediate quantities - to avoid division by zero
            amrex::ParticleReal const sinct = m_kt > 0 ? std::sin(m_kt * slice_ds) / m_kt : slice_ds;
            m_const_pt = sinct / betgam2;
        }

        /** This pushes a single particle, relative to the reference particle
         *
         * The @see compute_constants method must be called before pushing particles through this operator.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param idcpu particle global index
         * @param refpart reference particle (unused)
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            amrex::ParticleReal & AMREX_RESTRICT x,
            amrex::ParticleReal & AMREX_RESTRICT y,
            amrex::ParticleReal & AMREX_RESTRICT t,
            amrex::ParticleReal & AMREX_RESTRICT px,
            amrex::ParticleReal & AMREX_RESTRICT py,
            amrex::ParticleReal & AMREX_RESTRICT pt,
            uint64_t & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // shift due to alignment errors of the element
            shift_in(x, y, px, py);

            // initialize output values
            amrex::ParticleReal xout = x;
            amrex::ParticleReal yout = y;
            amrex::ParticleReal tout = t;
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal ptout = pt;

            // advance position and momentum
            xout = m_cos_kxds * x + m_sincx    * px;
            pxout = m_const_x * x + m_cos_kxds * px;

            yout = m_cos_kyds * y + m_sincy    * py;
            pyout = m_const_y * y + m_cos_kyds * py;

            tout = m_cos_ktds * t + m_const_pt * pt;
            ptout = m_const_t * t + m_cos_ktds * pt;

            // assign updated values
            x = xout;
            y = yout;
            t = tout;
            px = pxout;
            py = pyout;
            pt = ptout;

            // apply transverse aperture
            apply_aperture(x, y, idcpu);

            // undo shift due to alignment errors of the element
            shift_out(x, y, px, py);
        }

        /** This pushes the reference particle.
         *
         * @param[in,out] refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (RefPart & AMREX_RESTRICT refpart) const {

            using namespace amrex::literals; // for _rt and _prt

            // assign input reference particle values
            amrex::ParticleReal const x = refpart.x;
            amrex::ParticleReal const px = refpart.px;
            amrex::ParticleReal const y = refpart.y;
            amrex::ParticleReal const py = refpart.py;
            amrex::ParticleReal const z = refpart.z;
            amrex::ParticleReal const pz = refpart.pz;
            amrex::ParticleReal const t = refpart.t;
            amrex::ParticleReal const pt = refpart.pt;
            amrex::ParticleReal const s = refpart.s;

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // assign intermediate parameter
            amrex::ParticleReal const step = slice_ds / std::sqrt(amrex::Math::powi<2>(pt) - 1.0_prt);

            // advance position and momentum (straight element)
            refpart.x = x + step*px;
            refpart.y = y + step*py;
            refpart.z = z + step*pz;
            refpart.t = t - step*pt;

            // advance integrated path length
            refpart.s = s + slice_ds;

        }

        /** This pushes the covariance matrix. */
        using LinearTransport::operator();

        /** This function returns the linear transport map.
         *
         * @returns 6x6 transport matrix
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        Map6x6
        transport_map ([[maybe_unused]] RefPart const & AMREX_RESTRICT refpart) const
        {

            using namespace amrex::literals; // for _rt and _prt

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // find beta*gamma^2
            amrex::ParticleReal const betgam2 = amrex::Math::powi<2>(refpart.pt) - 1.0_prt;

            // trigo
	    // Note: The convention here is m_kx[1/m] = sign(k)*sqrt(abs(k))
	    // where k[1/m^2] is in the usual MADX convention

	    amrex::ParticleReal sin_kxds, cos_kxds, const_x, sincx;
	    if (m_kx > 0) {
	      std::tie(sin_kxds, cos_kxds) = amrex::Math::sincos(m_kx * slice_ds);
	      const_x = -m_kx * sin_kxds;
	      sincx   = sin_kxds / m_kx;
	    }
	    else if (m_kx < 0) {
	      sin_kxds = std::sinh(-m_kx*slice_ds);
	      cos_kxds = std::cosh(-m_kx*slice_ds);
	      const_x  = -m_kx * sin_kxds;
	      sincx    = -sin_kxds / m_kx;
	    }
	    else { //m_kx == 0
	      cos_kxds = 1.0;
	      const_x  = 0.0;
	      sincx    = slice_ds;
	    }

	    amrex::ParticleReal sin_kyds, cos_kyds, const_y, sincy;
	    if (m_ky > 0) {
	      std::tie(sin_kyds, cos_kyds) = amrex::Math::sincos(m_ky * slice_ds);
	      const_y = -m_ky * sin_kyds;
	      sincy   = sin_kyds / m_ky;
	    }
	    else if (m_ky < 0) {
	      sin_kyds = std::sinh(-m_ky*slice_ds);
	      cos_kyds = std::cosh(-m_ky*slice_ds);
	      const_y  = -m_ky * sin_kyds;
	      sincy    = -sin_kyds / m_ky;
	    }
	    else { //m_ky == 0
	      cos_kyds = 1.0;
	      const_y  = 0.0;
	      sincy    = slice_ds;
	    }

	    if (m_kt < 0) {
	      throw std::runtime_error(std::string(type) + ": must have kt >= 0!");
	    }
            auto const [sin_ktds, cos_ktds] = amrex::Math::sincos(m_kt * slice_ds);
            amrex::ParticleReal const_t = -m_kt * betgam2 * sin_ktds;
            // intermediate quantities - to avoid division by zero
            amrex::ParticleReal const sinct = m_kt > 0 ? std::sin(m_kt * slice_ds) / m_kt : slice_ds;
            amrex::ParticleReal const_pt = sinct / betgam2;

            // assign linear map matrix elements
            Map6x6 R = Map6x6::Identity();

            R(1,1) = cos_kxds;
            R(1,2) = sincx;
            R(2,1) = const_x;
            R(2,2) = cos_kxds;

            R(3,3) = cos_kyds;
            R(3,4) = sincy;
            R(4,3) = const_y;
            R(4,4) = cos_kyds;

            R(5,5) = cos_ktds;
            R(5,6) = const_pt;
            R(6,5) = const_t;
            R(6,6) = cos_ktds;

            return R;
        }

        amrex::ParticleReal m_kx; //! focusing x strength in 1/m
        amrex::ParticleReal m_ky; //! focusing y strength in 1/m
        amrex::ParticleReal m_kt; //! focusing t strength in 1/m

      private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_sincx, m_sincy;
        amrex::ParticleReal m_const_x, m_const_y, m_const_t, m_const_pt;
        amrex::ParticleReal m_cos_kxds, m_cos_kyds, m_cos_ktds;
    };

} // namespace impactx

#endif // IMPACTX_CONSTF_H
