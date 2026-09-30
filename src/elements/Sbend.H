/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_SBEND_H
#define IMPACTX_SBEND_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/gpuelement.H"
#include "mixin/alignment.H"
#include "mixin/pipeaperture.H"
#include "mixin/beamoptic.H"
#include "mixin/bendgeometry.H"
#include "mixin/nofinalize.H"
#include "mixin/lineartransport.H"
#include "mixin/spintransport.H"
#include "mixin/synchrotronradiation.H"

#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_SIMD.H>

#include <cmath>
#include <optional>

namespace impactx::elements
{
    struct SbendGPU
    : public mixin::BeamOptic<SbendGPU>,
      public mixin::LinearTransport<SbendGPU>,
      public mixin::ThickBend,
      public mixin::Alignment,
      public mixin::PipeAperture,
      public mixin::SpinTransport,
      public mixin::SynchrotronRadiation,
      public mixin::NoFinalize
      // At least on Intel AVX512, there is a small overhead to vectorize this element in DP, and no benefit in SP, see
      // https://github.com/BLAST-ImpactX/impactx/pull/1002
      // public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "Sbend";
        using PType = ImpactXParticleContainer::ParticleType;

        /** An ideal sector bend
         *
         * The bend is specified by exactly one of rc, phi or B, or by phi together with B,
         * @see mixin::ThickBend
         *
         * @param ds Segment length in m.
         * @param rc Radius of curvature in m.
         * @param phi Bend angle in degrees.
         * @param B Magnetic field in T.
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         * @param nslice number of slices used for the application of space charge
         */
        SbendGPU (
            amrex::ParticleReal ds,
            std::optional<amrex::ParticleReal> rc,
            std::optional<amrex::ParticleReal> phi = std::nullopt,
            std::optional<amrex::ParticleReal> B = std::nullopt,
            amrex::ParticleReal dx = DEFAULT_dx,
            amrex::ParticleReal dy = DEFAULT_dy,
            amrex::ParticleReal rotation_degree = DEFAULT_rotation_degree,
            amrex::ParticleReal aperture_x = DEFAULT_aperture_x,
            amrex::ParticleReal aperture_y = DEFAULT_aperture_y,
            int nslice = DEFAULT_nslice
        )
        : ThickBend(type, ds, nslice, rc, phi, B),
          Alignment(dx, dy, rotation_degree),
          PipeAperture(aperture_x, aperture_y)
        {
        }

        /** Reverse the element (negate ds and the bend angle) */
        void reverse () { ThickBend::reverse(); }

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
            using amrex::Math::powi;

            Alignment::compute_constants(refpart);

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // access reference particle values
            amrex::ParticleReal const ibet = 1.0_prt / refpart.beta();
            amrex::ParticleReal const ibetagam2 = 1_prt / powi<2>(refpart.beta_gamma());

            // compute intermediate constants (curvature, bend angle, trigonometry)
            bool const straight = is_straight();
            amrex::ParticleReal const rc = signed_rc(refpart);
            amrex::ParticleReal const h = straight ? 0_prt : 1_prt/rc;
            amrex::ParticleReal const theta = slice_angle(refpart);
            auto const [sin_theta, cos_theta] = amrex::Math::sincos(theta);

            // nonvanishing linear matrix elements (a straight bend is a drift)
            m_R11 = straight ? 1_prt : cos_theta;
            m_R12 = straight ? slice_ds : rc * sin_theta;
            m_R16 = straight ? 0_prt : -rc * ibet * (1_prt - cos_theta);
            m_R21 = straight ? 0_prt :-sin_theta * h;
            // m_R22 = m_R11
            m_R26 = straight ? 0_prt : -sin_theta * ibet;
            m_R34 = straight ? slice_ds : rc * theta;
            // m_R51 = -m_R26
            // m_R52 = -m_R16
            m_R56 = straight ? slice_ds * ibetagam2 : rc * (-theta + sin_theta * ibet * ibet);

            // access reference particle values for spin calculation
            amrex::ParticleReal G = refpart.gyromagnetic_anomaly;
            amrex::ParticleReal const gamma = refpart.gamma();
            amrex::ParticleReal const beta = refpart.beta();
            amrex::ParticleReal const gyro_const = 1_prt + G * gamma;

            // trigonometry for spin calculation
            amrex::ParticleReal const spin_ref_angle = G * gamma * theta;
            auto const [sin_Gtheta, cos_Gtheta] = amrex::Math::sincos(spin_ref_angle);

            // nonvanishing components of spin rotation at the design point
            m_lambday = -spin_ref_angle;

            // elements of the spin-orbit coupling matrix
            m_A14 = (gamma - 1_prt)/gamma * (1_prt - cos_Gtheta);
            m_A21 = -gyro_const * h * sin_theta;
            m_A22 = -gyro_const * (1_prt - cos_theta);
            m_A26 = -gyro_const * ibet * sin_theta + G * beta * gamma * theta;
            m_A34 = (gamma - 1_prt)/gamma * sin_Gtheta;

        }

        /** This is a sbend functor, so that a variable of this type can be used like a sbend function.
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
        template<typename T_Real=amrex::ParticleReal, typename T_IdCpu=uint64_t>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            T_Real & AMREX_RESTRICT x,
            T_Real & AMREX_RESTRICT y,
            T_Real & AMREX_RESTRICT t,
            T_Real & AMREX_RESTRICT px,
            T_Real const & AMREX_RESTRICT py,
            T_Real const & AMREX_RESTRICT pt,
            [[maybe_unused]] T_IdCpu const & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // initialize output values
            T_Real xout = x;
            T_Real yout = y;
            T_Real tout = t;

            // initialize output values of momenta
            T_Real pxout = px;
            // T_Real const pyout = py;
            // T_Real const ptout = pt;

            // advance position and momentum (sector bend)
            T_Real const R22 =  m_R11;
            T_Real const R51 = -m_R26;
            T_Real const R52 = -m_R16;

            xout = m_R11 * x
                 + m_R12 * px
                 + m_R16 * pt;

            pxout = m_R21 * x
                  +   R22 * px
                  + m_R26 * pt;

            yout =         y
                 + m_R34 * py;

            // pyout = py;

            tout =   R51 * x
                 +   R52 * px
                 +         t
                 + m_R56 * pt;

            // ptout = pt;

            // assign updated values
            x = xout;
            y = yout;
            t = tout;
            px = pxout;
            // py = pyout;
            // pt = ptout;
        }

        /** This pushes the reference particle. */
        using ThickBend::operator();

        /** This pushes the covariance matrix. */
        using LinearTransport::operator();

        /** This function returns the linear transport map.
         *
         * @returns 6x6 transport matrix
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        Map6x6
        transport_map (RefPart const & AMREX_RESTRICT refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // access reference particle values to find beta*gamma^2
            amrex::ParticleReal const pt_ref = refpart.pt;
            amrex::ParticleReal const betgam2 = powi<2>(pt_ref) - 1.0_prt;
            amrex::ParticleReal const bet = std::sqrt(betgam2/(1.0_prt + betgam2));

            // initialize linear map matrix elements
            Map6x6 R = Map6x6::Identity();

            // treat the special case of a straight bend (drift)
            if (is_straight()) {
               R(1,2) = slice_ds;
               R(3,4) = slice_ds;
               R(5,6) = slice_ds / betgam2;

            } else {

               // calculate expensive terms once
               amrex::ParticleReal const rc = signed_rc(refpart);
               amrex::ParticleReal const theta = slice_angle(refpart);
               auto const [sin_theta, cos_theta] = amrex::Math::sincos(theta);

               // assign linear map matrix elements
               R(1,1) = cos_theta;
               R(1,2) = rc*sin_theta;
               R(1,6) = - (rc/bet)*(1.0_prt - cos_theta);
               R(2,1) = -sin_theta/rc;
               R(2,2) = cos_theta;
               R(2,6) = - sin_theta/bet;
               R(3,4) = rc*theta;
               R(5,1) = sin_theta/bet;
               R(5,2) = rc/bet*(1.0_prt - cos_theta);
               R(5,6) = rc*(-theta+sin_theta/(bet*bet));

            }

            // apply the transverse rotation (roll) alignment error
            return rotate_aligned_map(R);
        }

        /** This operator pushes the particle spin.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param idcpu particle global index
         * @param refpart reference particle (unused)
         * @param sx particle spin x-component
         * @param sy particle spin y-component
         * @param sz particle spin z-component
         */
        template<typename T_Real=amrex::ParticleReal, typename T_IdCpu=uint64_t>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void spin_and_phasespace_push (
            T_Real & AMREX_RESTRICT x,
            T_Real & AMREX_RESTRICT y,
            T_Real & AMREX_RESTRICT t,
            T_Real & AMREX_RESTRICT px,
            T_Real const & AMREX_RESTRICT py,
            T_Real const & AMREX_RESTRICT pt,
            T_Real & AMREX_RESTRICT sx,
            T_Real & AMREX_RESTRICT sy,
            T_Real & AMREX_RESTRICT sz,
            T_IdCpu const & AMREX_RESTRICT idcpu,
            RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // initialize the three components of the axis-angle vector
            T_Real lambdax = 0_prt;
            T_Real lambday = 0_prt;
            T_Real lambdaz = 0_prt;

            // set the angle-axis generator based on the phase space variables
            lambdax = m_A14 * py;
            lambday = m_A21 * x + m_A22 * px + m_A26 * pt;
            lambdaz = m_A34 * py;

            // push the spin vector using the generator just determined
            rotate_spin(lambdax,lambday,lambdaz,sx,sy,sz);

            // axis-angle vector components generating the reference spin map
            lambdax = 0.0_prt;
            lambday = m_lambday;
            lambdaz = 0.0_prt;

            // push the spin vector using the generator just determined
            rotate_spin(lambdax,lambday,lambdaz,sx,sy,sz);

            // phase space push
            (*this)(x, y, t, px, py, pt, idcpu, refpart);
        }

      private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_R11, m_R12, m_R16, m_R21, m_R26, m_R34, m_R56;
        amrex::ParticleReal m_lambday, m_A14, m_A21, m_A22, m_A26, m_A34;
    };

    /** Sbend: GPU state plus common host-only metadata */
    using Sbend = Element<SbendGPU>;

} // namespace impactx

IMPACTX_PUSH_EXTERN_TEMPLATE(impactx::elements::SbendGPU)

#endif // IMPACTX_SBEND_H
