/* Copyright 2022-2024 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENT_LINEAR_MAP_H
#define IMPACTX_ELEMENT_LINEAR_MAP_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/alignment.H"
#include "mixin/beamoptic.H"
#include "mixin/lineartransport.H"
#include "mixin/named.H"
#include "mixin/nofinalize.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <stdexcept>


namespace impactx::elements
{
    struct LinearMap
    : public mixin::Named,
      public mixin::BeamOptic<LinearMap>,
      public mixin::LinearTransport<LinearMap>,
      public mixin::Alignment,
      public mixin::NoFinalize
    {
        static constexpr auto type = "LinearMap";
        using PType = ImpactXParticleContainer::ParticleType;

        /** A thin element that applies a user-provided linear transport map R
         *  to the 6-vector of phase space coordinates (x,px,y,py,t,pt).
         *  Thus  x_final = R(1,1)*x + R(1,2)*px + R(1,3)*y + ...,
         *       px_final = R(2,1)*x + R(2,2)*px + R(2,3)*y + ..., etc.
         *
         * @param R user-provided transport map
         * @param ds Segment length in m
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param name a user defined and not necessarily unique name of the element
         */
        LinearMap (
            Map6x6 const & R,
            amrex::ParticleReal ds = 0,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0,
            std::optional<std::string> name = std::nullopt
        )
        : Named(std::move(name)),
          Alignment(dx, dy, rotation_degree)
        {
            m_transport_map = R;
            m_ds = ds;
        }

        /** Push all particles */
        using BeamOptic::operator();

        /** This is a LinearMap functor, so that a variable of this type can be used like a
         *  LinearMap function.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t (unused)
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t (unused)
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
            [[maybe_unused]] uint64_t & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // shift due to alignment errors of the element
            shift_in(x, y, px, py);

            // input and output phase space vectors
            amrex::SmallVector<amrex::ParticleReal, 6, 1> vectorin{
                x, px, y, py, t, pt
            };

            amrex::SmallVector<amrex::ParticleReal, 6, 1> const vectorout = m_transport_map * vectorin;

            // assign updated values
            x = vectorout(1);
            px = vectorout(2);
            y = vectorout(3);
            py = vectorout(4);
            t = vectorout(5);
            pt = vectorout(6);

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
            if (m_ds > 0)  // Drift
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
                amrex::ParticleReal const step = slice_ds /std::sqrt(std::pow(pt,2)-1.0_prt);

                // advance position and momentum (drift)
                refpart.x = x + step*px;
                refpart.y = y + step*py;
                refpart.z = z + step*pz;
                refpart.t = t - step*pt;

                // advance integrated path length
                refpart.s = s + slice_ds;
            }
            // else nothing to do for a zero-length element
        }

        /** Number of slices used for the application of space charge
         *
         * @return one, because we do not support slicing of this element
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        int nslice () const
        {
            return 1;
        }

        /** Return the segment length
         *
         * @return by default zero, but users can set a corresponding ds for bookkeeping
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal ds () const
        {
            return m_ds;
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

        Map6x6 m_transport_map;  // 6x6 transport map
        amrex::ParticleReal m_ds;  // finite ds allowed for bookkeeping, but we do not allow slicing
    };

} // namespace impactx

#endif // IMPACTX_ELEMENT_LINEAR_MAP_H
