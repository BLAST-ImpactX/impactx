/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
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
#include "mixin/alignment.H"
#include "mixin/pipeaperture.H"
#include "mixin/beamoptic.H"
#include "mixin/thick.H"
#include "mixin/named.H"
#include "mixin/nofinalize.H"
#include "mixin/lineartransport.H"

#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx::elements
{
    struct Sbend
    : public mixin::Named,
      public mixin::BeamOptic<Sbend>,
      public mixin::LinearTransport<Sbend>,
      public mixin::Thick,
      public mixin::Alignment,
      public mixin::PipeAperture,
      public mixin::NoFinalize
    {
        static constexpr auto type = "Sbend";
        using PType = ImpactXParticleContainer::ParticleType;

        /** An ideal sector bend
         *
         * @param ds Segment length in m.
         * @param rc Radius of curvature in m.
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         * @param nslice number of slices used for the application of space charge
         * @param name a user defined and not necessarily unique name of the element
         */
        Sbend (
            amrex::ParticleReal ds,
            amrex::ParticleReal rc,
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
          m_rc(rc)
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

            // access reference particle values
            amrex::ParticleReal const ibet = 1.0_prt / refpart.beta();
            amrex::ParticleReal const ibetagam2 = 1_prt / std::pow(refpart.beta_gamma(), 2);

            // treat the special case of zero field (equivalent to a drift)
            if (m_rc==0_prt) {
               m_R11 = 1.0_prt;
               m_R12 = slice_ds;
               m_R16 = 0.0_prt;
               m_R21 = 0.0_prt;
               m_R26 = 0.0_prt;
               m_R34 = slice_ds;
               m_R56 = slice_ds * ibetagam2;

            } else {

               // calculate expensive terms once
               amrex::ParticleReal const theta = slice_ds / m_rc;
               auto const [sin_theta, cos_theta] = amrex::Math::sincos(theta);

               m_R11 = cos_theta;
               m_R12 = m_rc * sin_theta;
               m_R16 = -m_rc * ibet * (1.0_prt - cos_theta);
               m_R21 = -sin_theta / m_rc;
               // m_R22 = m_R11
               m_R26 = -sin_theta * ibet;
               m_R34 = m_rc * theta;
               // m_R51 = -m_R26
               // m_R52 = -m_R16
               m_R56 = m_rc * (-theta + sin_theta * ibet * ibet);
           }
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

            // intialize output values
            amrex::ParticleReal xout = x;
            amrex::ParticleReal yout = y;
            amrex::ParticleReal tout = t;

            // initialize output values of momenta
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal const pyout = py;
            amrex::ParticleReal const ptout = pt;

            // advance position and momentum (sector bend)
            amrex::ParticleReal const R22 =  m_R11;
            amrex::ParticleReal const R51 = -m_R26;
            amrex::ParticleReal const R52 = -m_R16;

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

            // treat the special case of zero field (drift)
            if (m_rc == 0.0_prt) {
               // advance position and momentum (drift)
               amrex::ParticleReal const step = slice_ds /std::sqrt(std::pow(pt,2)-1.0_prt);           
               refpart.x = x + step*px;
               refpart.y = y + step*py;
               refpart.z = z + step*pz;
               refpart.t = t - step*pt;

            } else {

               // assign intermediate parameter
               amrex::ParticleReal const theta = slice_ds/m_rc;
               amrex::ParticleReal const B = std::sqrt(std::pow(pt,2)-1.0_prt)/m_rc;

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
         * @returns 6x6 transport matrix
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        Map6x6
        transport_map (RefPart const & AMREX_RESTRICT refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // access reference particle values to find beta*gamma^2
            amrex::ParticleReal const pt_ref = refpart.pt;
            amrex::ParticleReal const betgam2 = std::pow(pt_ref, 2) - 1.0_prt;
            amrex::ParticleReal const bet = std::sqrt(betgam2/(1.0_prt + betgam2));

            // calculate expensive terms once
            amrex::ParticleReal const theta = slice_ds/m_rc;
            auto const [sin_theta, cos_theta] = amrex::Math::sincos(theta);

            // initialize linear map matrix elements
            Map6x6 R = Map6x6::Identity();

            // assign linear map matrix elements
            R(1,1) = cos_theta;
            R(1,2) = m_rc*sin_theta;
            R(1,6) = - (m_rc/bet)*(1.0_prt - cos_theta);
            R(2,1) = -sin_theta/m_rc;
            R(2,2) = cos_theta;
            R(2,6) = - sin_theta/bet;
            R(3,4) = m_rc*theta;
            R(5,1) = sin_theta/bet;
            R(5,2) = m_rc/bet*(1.0_prt - cos_theta);
            R(5,6) = m_rc*(-theta+sin_theta/(bet*bet));

            return R;
        }

        amrex::ParticleReal m_rc; //! bend radius in m

      private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_R11, m_R12, m_R16, m_R21, m_R26, m_R34, m_R56;
    };

} // namespace impactx

#endif // IMPACTX_SBEND_H
