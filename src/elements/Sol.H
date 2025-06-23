/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_SOL_H
#define IMPACTX_SOL_H

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
    struct Sol
    : public mixin::Named,
      public mixin::BeamOptic<Sol>,
      public mixin::LinearTransport<Sol>,
      public mixin::Thick,
      public mixin::Alignment,
      public mixin::PipeAperture,
      public mixin::NoFinalize
    {
        static constexpr auto type = "Sol";
        using PType = ImpactXParticleContainer::ParticleType;

        /** An ideal hard-edge Solenoid magnet
         *
         * @param ds Segment length in m.
         * @param ks Solenoid strength in m^(-1) (MADX convention)
         *           = (magnetic field Bz in T) / (rigidity in T-m)
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         * @param nslice number of slices used for the application of space charge
         * @param name a user defined and not necessarily unique name of the element
         */
        Sol (
            amrex::ParticleReal ds,
            amrex::ParticleReal ks,
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
          m_ks(ks)
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

            // access reference particle values to find beta*gamma^2
            amrex::ParticleReal const betgam2 = amrex::Math::powi<2>(refpart.pt) - 1.0_prt;

            // compute phase advance per unit length (in rad/m) and
            // rotation angle (in rad)
            amrex::ParticleReal const alpha = m_ks / 2.0_prt;
            amrex::ParticleReal const theta = alpha*slice_ds;
            auto const [sin_theta, cos_theta] = amrex::Math::sincos(theta);

            m_sin_theta = sin_theta;
            m_cos_theta = cos_theta;
            m_sin_theta_ialpha = sin_theta / alpha;
            m_neg_alpha_sin_theta = -alpha * sin_theta;
            m_ds_bg2 = slice_ds / betgam2;
        }

        /** This is a sol functor, so that a variable of this type can be used like a sol function.
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
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal ptout = pt;

            // advance positions and momenta using map for focusing
            xout = m_cos_theta        * x
                 + m_sin_theta_ialpha * px;
            pxout = m_neg_alpha_sin_theta * x
                  + m_cos_theta           * px;

            yout = m_cos_theta        * y
                 + m_sin_theta_ialpha * py;
            pyout = m_neg_alpha_sin_theta * y
                  + m_cos_theta           * py;

            tout = t + m_ds_bg2 * pt;
            ptout = pt;

            // assign intermediate momenta
            px = pxout;
            py = pyout;
            pt = ptout;

            // advance positions and momenta using map for rotation
            x = m_cos_theta * xout
              + m_sin_theta * yout;
            pxout = m_cos_theta * px
                  + m_sin_theta * py;

            y = -m_sin_theta * xout
               + m_cos_theta * yout;
            pyout = -m_sin_theta * px
                   + m_cos_theta * py;

            t = tout;
            ptout = pt;

            // assign updated values
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

            // initialize linear map matrix elements
            Map6x6 R = Map6x6::Identity();
            Map6x6 Rrotation = Map6x6::Identity();
            Map6x6 Rfocusing = Map6x6::Identity();

            // access reference particle values to find beta*gamma^2
            amrex::ParticleReal const betgam2 = amrex::Math::powi<2>(refpart.pt) - 1.0_prt;

            // compute phase advance per unit length (in rad/m) and
            // rotation angle (in rad)
            amrex::ParticleReal const alpha = m_ks / 2.0_prt;
            amrex::ParticleReal const theta = alpha*slice_ds;
            auto const [sin_theta, cos_theta] = amrex::Math::sincos(theta);

            amrex::ParticleReal const sin_theta_ialpha = sin_theta / alpha;
            amrex::ParticleReal const neg_alpha_sin_theta = -alpha * sin_theta;
            amrex::ParticleReal const ds_bg2 = slice_ds / betgam2;

            Rfocusing(1,1) = cos_theta;
            Rfocusing(1,2) = sin_theta_ialpha;
            Rfocusing(2,1) = neg_alpha_sin_theta;
            Rfocusing(2,2) = cos_theta;
            Rfocusing(3,3) = cos_theta;
            Rfocusing(3,4) = sin_theta_ialpha;
            Rfocusing(4,3) = neg_alpha_sin_theta;
            Rfocusing(4,4) = cos_theta;
            Rfocusing(5,6) = ds_bg2;

            Rrotation(1,1) = cos_theta;
            Rrotation(1,3) = sin_theta;
            Rrotation(2,2) = cos_theta;
            Rrotation(2,4) = sin_theta;
            Rrotation(3,1) = -sin_theta;
            Rrotation(3,3) = cos_theta;
            Rrotation(4,2) = -sin_theta;
            Rrotation(4,4) = cos_theta;

            R = Rrotation * Rfocusing;

            return R;
        }

        amrex::ParticleReal m_ks; //! solenoid strength in 1/m

      private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_ds_bg2, m_cos_theta, m_sin_theta, m_sin_theta_ialpha, m_neg_alpha_sin_theta;
    };

} // namespace impactx

#endif // IMPACTX_SOL_H
