/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_CFBEND_H
#define IMPACTX_CFBEND_H

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
    struct CFbend
    : public mixin::Named,
      public mixin::BeamOptic<CFbend>,
      public mixin::LinearTransport<CFbend>,
      public mixin::Thick,
      public mixin::Alignment,
      public mixin::PipeAperture,
      public mixin::NoFinalize
    {
        static constexpr auto type = "CFbend";
        using PType = ImpactXParticleContainer::ParticleType;

        /** An combined-function bend, consisting of an ideal sector bend with
         *  an upright quadrupole focusing component.
         *
         * @param ds Segment length in m.
         * @param rc Radius of curvature in m.
         * @param k  Quadrupole strength in m^(-2) (MADX convention)
         *           = (gradient in T/m) / (rigidity in T-m)
         *           k > 0 horizontal focusing
         *           k < 0 horizontal defocusing
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         * @param nslice number of slices used for the application of space charge
         * @param name a user defined and not necessarily unique name of the element
         */
        CFbend (
            amrex::ParticleReal ds,
            amrex::ParticleReal rc,
            amrex::ParticleReal k,
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
          m_rc(rc), m_k(k)
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

            // find beta*gamma^2, beta
            amrex::ParticleReal const betgam2 = amrex::Math::powi<2>(refpart.pt) - 1_prt;
            amrex::ParticleReal const bet = refpart.beta();
            amrex::ParticleReal const ibetgam2 = 1_prt / betgam2;
            amrex::ParticleReal const b2rc2 = amrex::Math::powi<2>(bet) * amrex::Math::powi<2>(m_rc);

            // update horizontal and longitudinal phase space variables
            amrex::ParticleReal const gx = m_k + amrex::Math::powi<-2>(m_rc);
            amrex::ParticleReal const omega_x = std::sqrt(std::abs(gx));

            // update vertical phase space variables
            amrex::ParticleReal const gy = -m_k;
            amrex::ParticleReal const omega_y = std::sqrt(std::abs(gy));

            // trigonometry
            auto const [sinx, cosx] = amrex::Math::sincos(omega_x * slice_ds);
            amrex::ParticleReal const sinhx = std::sinh(omega_x * slice_ds);
            amrex::ParticleReal const coshx = std::cosh(omega_x * slice_ds);
            auto const [siny, cosy] = amrex::Math::sincos(omega_y * slice_ds);
            amrex::ParticleReal const sinhy = std::sinh(omega_y * slice_ds);
            amrex::ParticleReal const coshy = std::cosh(omega_y * slice_ds);

            amrex::ParticleReal const igbrc = 1_prt / (gx * bet * m_rc);
            amrex::ParticleReal const iobrc = 1_prt / (omega_x * bet * m_rc);
            amrex::ParticleReal const igobr = 1_prt / (gx * omega_x * b2rc2);

            m_R11 = gx > 0_prt ? cosx : coshx;
            m_R12 = gx > 0_prt ? sinx / omega_x : sinhx / omega_x;
            m_R16 = gx > 0_prt ? -(1_prt - cosx) * igbrc : -(1_prt - coshx) * igbrc;
            m_R21 = gx > 0_prt ? -omega_x * sinx : omega_x * sinhx;
            m_R22 = gx > 0_prt ? cosx : coshx;
            m_R26 = gx > 0_prt ? -sinx * iobrc : -sinhx * iobrc;
            m_R33 = gy > 0_prt ? cosy : coshy;
            m_R34 = gy > 0_prt ? siny / omega_y : sinhy / omega_y;
            m_R43 = gy > 0_prt ? -omega_y * siny : omega_y * sinhy;
            m_R44 = gy > 0_prt ? cosy : coshy;
            m_R51 = gx > 0_prt ? sinx * iobrc : sinhx * iobrc;
            m_R52 = gx > 0_prt ? (1_prt - cosx) * igbrc : (1_prt - coshx) * igbrc;
            m_R56 = gx > 0_prt ?
                slice_ds * ibetgam2 + (sinx  - omega_x * slice_ds) * igobr :
                slice_ds * ibetgam2 + (sinhx - omega_x * slice_ds) * igobr;
        }

        /** This is a cfbend functor, so that a variable of this type can be used like a cfbend function.
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

            // initialize output values of momenta
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal ptout = pt;

            // update horizontal and longitudinal phase space variables
            //   advance position and momentum: (de)focusing
            x =     m_R11 * xout + m_R12 * px + m_R16 * pt;
            pxout = m_R21 * xout + m_R22 * px + m_R26 * pt;
            t =     m_R51 * xout + m_R52 * px + m_R56 * pt + tout;
            ptout = pt;

            // update vertical phase space variables
            //   advance position and momentum (de)focusing
            y =     m_R33 * yout + m_R34 * py;
            pyout = m_R43 * yout + m_R44 * py;

            // assign updated momenta
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
        void operator() (RefPart & AMREX_RESTRICT refpart) const
        {
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
            amrex::ParticleReal const theta = slice_ds/m_rc;
            amrex::ParticleReal const B = std::sqrt(amrex::Math::powi<2>(pt)-1.0_prt)/m_rc;

            // calculate expensive terms once
            auto const [sin_theta, cos_theta] = amrex::Math::sincos(theta);

            // advance position and momentum (bend)
            refpart.px = px*cos_theta - pz*sin_theta;
            refpart.py = py;
            refpart.pz = pz*cos_theta + px*sin_theta;
            refpart.pt = pt;

            refpart.x = x + (refpart.pz - pz)/B;
            refpart.y = y + (theta/B)*py;
            refpart.z = z - (refpart.px - px)/B;
            refpart.t = t - (theta/B)*pt;

            // advance integrated path length
            refpart.s = s + slice_ds;
        }

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

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // find beta*gamma^2, beta
            amrex::ParticleReal const betgam2 = amrex::Math::powi<2>(refpart.pt) - 1_prt;
            amrex::ParticleReal const bet = refpart.beta();
            amrex::ParticleReal const ibetgam2 = 1_prt / betgam2;
            amrex::ParticleReal const b2rc2 = amrex::Math::powi<2>(bet) * amrex::Math::powi<2>(m_rc);

            // update horizontal and longitudinal phase space variables
            amrex::ParticleReal const gx = m_k + amrex::Math::powi<-2>(m_rc);
            amrex::ParticleReal const omega_x = std::sqrt(std::abs(gx));

            // update vertical phase space variables
            amrex::ParticleReal const gy = -m_k;
            amrex::ParticleReal const omega_y = std::sqrt(std::abs(gy));

            // trigonometry
            auto const [sinx, cosx] = amrex::Math::sincos(omega_x * slice_ds);
            amrex::ParticleReal const sinhx = std::sinh(omega_x * slice_ds);
            amrex::ParticleReal const coshx = std::cosh(omega_x * slice_ds);
            auto const [siny, cosy] = amrex::Math::sincos(omega_y * slice_ds);
            amrex::ParticleReal const sinhy = std::sinh(omega_y * slice_ds);
            amrex::ParticleReal const coshy = std::cosh(omega_y * slice_ds);

            amrex::ParticleReal const igbrc = 1_prt / (gx * bet * m_rc);
            amrex::ParticleReal const iobrc = 1_prt / (omega_x * bet * m_rc);
            amrex::ParticleReal const igobr = 1_prt / (gx * omega_x * b2rc2);

            // initialize linear map matrix elements
            Map6x6 R = Map6x6::Identity();

            R(1,1) = gx > 0_prt ? cosx : coshx;
            R(1,2) = gx > 0_prt ? sinx / omega_x : sinhx / omega_x;
            R(1,6) = gx > 0_prt ? -(1_prt - cosx) * igbrc : -(1_prt - coshx) * igbrc;
            R(2,1) = gx > 0_prt ? -omega_x * sinx : omega_x * sinhx;
            R(2,2) = gx > 0_prt ? cosx : coshx;
            R(2,6) = gx > 0_prt ? -sinx * iobrc : -sinhx * iobrc;
            R(3,3) = gy > 0_prt ? cosy : coshy;
            R(3,4) = gy > 0_prt ? siny / omega_y : sinhy / omega_y;
            R(4,3) = gy > 0_prt ? -omega_y * siny : omega_y * sinhy;
            R(4,4) = gy > 0_prt ? cosy : coshy;
            R(5,1) = gx > 0_prt ? sinx * iobrc : sinhx * iobrc;
            R(5,2) = gx > 0_prt ? (1_prt - cosx) * igbrc : (1_prt - coshx) * igbrc;
            R(5,6) = gx > 0_prt ?
                slice_ds * ibetgam2 + (sinx  - omega_x * slice_ds) * igobr :
                slice_ds * ibetgam2 + (sinhx - omega_x * slice_ds) * igobr;

            return R;
        }

        amrex::ParticleReal m_rc; //! bend radius in m
        amrex::ParticleReal m_k;  //! quadrupole strength in m^(-2)

    private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_R11, m_R12, m_R16, m_R21, m_R22, m_R26, m_R33, m_R34, m_R43, m_R44, m_R51, m_R52, m_R56;
    };

} // namespace impactx

#endif // IMPACTX_CFBEND_H
