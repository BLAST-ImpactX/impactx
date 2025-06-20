/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_EXACTSBEND_H
#define IMPACTX_EXACTSBEND_H

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
#include <AMReX_SIMD.H>

#include <cmath>
#include <stdexcept>


namespace impactx::elements
{
    struct ExactSbend
    : public mixin::Named,
      public mixin::BeamOptic<ExactSbend>,
      public mixin::LinearTransport<ExactSbend>,
      public mixin::Thick,
      public mixin::Alignment,
      public mixin::PipeAperture,
      public mixin::NoFinalize,
      public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "ExactSbend";
        static constexpr amrex::ParticleReal degree2rad = ablastr::constant::math::pi / 180.0;
        using PType = ImpactXParticleContainer::ParticleType;

        /** The body of an ideal sector bend, using the exact nonlinear transfer map.
         *  The map corresponds to the map described in:
         *  D. L. Bruhwiler et al, in Proc. of EPAC 98, pp. 1171-1173 (1998).
         *
         *  In the ultrarelativistic limit, this is equivalent to the map described in:
         *  E. Forest et al, Part. Accel. 45, pp. 65-94 (1994).
         *
         *  The model consists of a uniform bending field B_y with a hard edge.
         *  Pole faces are normal to the entry and exit velocity of the reference particle.
         *
         * @param ds Segment length in m.
         * @param phi Segment bend angle in degrees.
         * @param B Magnetic field strength in T.
         *  When B = 0 (default), the reference bending radius is defined by r0 = length / (angle in rad), corresponding to
         *  a magnetic field of B = rigidity / r0; otherwise the reference bending radius is defined by r0 = rigidity / B.
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         * @param nslice number of slices used for the application of space charge
         * @param name a user defined and not necessarily unique name of the element
         */
        ExactSbend (
            amrex::ParticleReal ds,
            amrex::ParticleReal phi,
            amrex::ParticleReal B,
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
          m_phi(phi * degree2rad), m_B(B)
        {
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal
        rc (RefPart const & refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt

            return m_B != 0_prt ? refpart.rigidity_Tm() / m_B : m_ds / m_phi;
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
            using amrex::Math::powi;

            Alignment::compute_constants(refpart);

            // reference particle's orbital radius
            m_rc = this->rc(refpart);

            // arc length and angle of the current slice
            m_slice_ds = m_ds / nslice();
            amrex::ParticleReal const bend_sign = m_ds * m_rc / std::abs(m_ds * m_rc);
            m_slice_phi = bend_sign * std::abs(m_phi) / nslice();  // the sign is determined by rc/ds

            // trigonometric evaluations
            auto const [sin_phi, cos_phi] = amrex::Math::sincos(m_slice_phi);
            m_sin_phi = sin_phi;
            m_cos_phi = cos_phi;

            // access reference particle values to find relativistic factors
            m_ibeta = 1_prt / refpart.beta();
            m_ibetagam2 = 1_prt / powi<2>(refpart.beta_gamma());

        }

        /** This is an ExactSbend functor, so that a variable of this type can be used like an ExactSbend function.
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
            T_Real & AMREX_RESTRICT py,
            T_Real & AMREX_RESTRICT pt,
            T_IdCpu & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;
            using namespace std; // for cmath(float)
            namespace stdx = amrex::simd::stdx;

            // shift due to alignment errors of the element
            shift_in(x, y, px, py);

            // access position data
            T_Real const xout = x;
            T_Real const yout = y;
            T_Real const tout = t;

            // initialize output values of momenta
            T_Real pxout = px;
            T_Real pyout = py;
            T_Real ptout = pt;

            // treat the special case of zero field
            if (m_phi==0_prt && m_B==0_prt) {
               // compute the radical in the denominator (= pz):
               T_Real const inv_pzden = 1_prt / sqrt(
                   powi<2>(pt - m_ibeta) -
                   m_ibetagam2 -
                   powi<2>(px) -
                   powi<2>(py)
               );

               // advance position and momentum (exact drift)
               x = xout + m_slice_ds * px * inv_pzden;
               // pxout = px;
               y = yout + m_slice_ds * py * inv_pzden;
               // pyout = py;
               t = tout - m_slice_ds * (m_ibeta + (pt - m_ibeta) * inv_pzden);
               // ptout = pt;

               // assign updated momenta
               px = pxout;
               py = pyout;
               pt = ptout;

            } else {

               // assign intermediate quantities
               T_Real const pperp2 = powi<2>(pt)-2.0_prt * m_ibeta * pt - powi<2>(py)+1.0_prt;
               T_Real const px2 = powi<2>(px);
               // determine if particle lies within the domain of map definition
               auto const mask = pperp2 <= px2;
               if (stdx::any_of(mask) == true)
               {
                   amrex::ParticleIDWrapper<T_IdCpu>{idcpu}.make_invalid(mask);
               } else
               {
                   T_Real const pperp = sqrt(pperp2);
                   T_Real const pzi = sqrt(powi<2>(pperp) - powi<2>(px));
                   T_Real const rho = m_rc + xout;

                   // update momenta
                   pxout = px * m_cos_phi + (pzi - rho / m_rc) * m_sin_phi;
                   pyout = py;
                   ptout = pt;

                   // angle of momentum rotation
                   T_Real const px2f = powi<2>(pxout);
                   // determine if particle lies within domain of map definition
                   auto const mask2 = pperp2 <= px2f;
                   if (stdx::any_of(mask2) == true)
                   {
                       amrex::ParticleIDWrapper<T_IdCpu>{idcpu}.make_invalid(mask2);
                   } else
                   {
                       T_Real const pzf = sqrt(powi<2>(pperp)-powi<2>(pxout));
                       T_Real const theta = m_slice_phi + asin(px/pperp) - asin(pxout/pperp);

                       // update position coordinates
                       x = -m_rc + rho*m_cos_phi + m_rc*(pzf + px*m_sin_phi - pzi*m_cos_phi);
                       y = yout + theta * m_rc * py;
                       t = tout - theta * m_rc * (pt - 1.0_prt * m_ibeta) - m_slice_phi * m_rc * m_ibeta;

                       // assign updated momenta
                       px = pxout;
                       py = pyout;
                       pt = ptout;
                   }
               }
            }

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
            using amrex::Math::powi;

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

            // special case of zero field = an exact drift
            if (m_phi==0_prt && m_B==0_prt) {
               // advance position and momentum (drift)
               amrex::ParticleReal const step = slice_ds /std::sqrt(powi<2>(pt)-1.0_prt);
               refpart.x = x + step*px;
               refpart.y = y + step*py;
               refpart.z = z + step*pz;
               refpart.t = t - step*pt;

            } else {

               // assign intermediate parameters
               amrex::ParticleReal const rc = (m_B != 0_prt) ? refpart.rigidity_Tm() / m_B : m_ds / m_phi;
               amrex::ParticleReal const B = refpart.beta_gamma() /rc;
               amrex::ParticleReal const bend_sign = m_ds * rc / std::abs(m_ds * rc);
               amrex::ParticleReal const theta = bend_sign * std::abs(m_phi) / nslice(); // sign is determined by rc*ds

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

            }

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
            throw std::runtime_error(std::string(type) + ": Envelope tracking is not yet implemented!");
            return Map6x6::Identity();
        }

        amrex::ParticleReal m_phi; //! bend angle in radians
        amrex::ParticleReal m_B;  //! magnetic field in T

      private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_slice_ds, m_slice_phi, m_ibeta, m_ibetagam2, m_rc, m_sin_phi, m_cos_phi;
    };

} // namespace impactx

#endif // IMPACTX_EXACTSBEND_H
