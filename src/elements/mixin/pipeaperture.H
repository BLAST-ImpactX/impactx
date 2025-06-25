/* Copyright 2022-2025 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell, Alexander Sinn
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_MIXIN_PIPE_APERTURE_H
#define IMPACTX_ELEMENTS_MIXIN_PIPE_APERTURE_H

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
#include <AMReX_Particle.H>
#include <AMReX_REAL.H>


namespace impactx::elements::mixin
{
    /** This is a helper class for applying a transverse aperture restriction to thick lattice elements
     */
    struct PipeAperture
    {
        /** A finite-length element with a constant elliptical aperture over s.
         *
         * @param aperture_x horizontal half-aperture size [m]
         * @param aperture_y vertical half-aperture size [m]
         */
        PipeAperture (
            amrex::ParticleReal aperture_x,
            amrex::ParticleReal aperture_y
        )
        {
            using namespace amrex::literals;  // for _prt

            m_inv_aperture_x = aperture_x > 0_prt ? 1_prt / aperture_x : 0_prt;
            m_inv_aperture_y = aperture_y > 0_prt ? 1_prt / aperture_y : 0_prt;
        }

        PipeAperture () = default;
        PipeAperture (PipeAperture const &) = default;
        PipeAperture& operator= (PipeAperture const &) = default;
        PipeAperture (PipeAperture&&) = default;
        PipeAperture& operator= (PipeAperture&& rhs) = default;

        ~PipeAperture () = default;

        /** Apply the transverse aperture
         *
         * @param[inout] x horizontal position relative to reference particle
         * @param[inout] y vertical position relative to reference particle
         * @param idcpu particle global index
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void apply_aperture (
            amrex::ParticleReal & AMREX_RESTRICT x,
            amrex::ParticleReal & AMREX_RESTRICT y,
            uint64_t & AMREX_RESTRICT idcpu
        ) const {
            using namespace amrex::literals; // for _rt and _prt

            // skip aperture application if aperture_x <= 0 or aperture_y <= 0
            if (m_inv_aperture_x > 0_prt && m_inv_aperture_y > 0_prt) {

               // scale horizontal and vertical coordinates
               amrex::ParticleReal const u = x * m_inv_aperture_x;
               amrex::ParticleReal const v = y * m_inv_aperture_y;

               // compare against the aperture boundary
               if (amrex::Math::powi<2>(u) + amrex::Math::powi<2>(v) > 1_prt) {
                  amrex::ParticleIDWrapper{idcpu}.make_invalid();
               }
            }
        }

        /** Horizontal aperture size
         *
         * @return horizontal aperture size [m]
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal aperture_x () const
        {
            using namespace amrex::literals; // for _prt
            return 1_prt / m_inv_aperture_x;
        }

        /** Vertical aperture size
         *
         * @return vertical aperture size [m]
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal aperture_y () const
        {
            using namespace amrex::literals; // for _prt
            return 1_prt / m_inv_aperture_y;
        }

        amrex::ParticleReal m_inv_aperture_x = 0; //! inverse of the horizontal aperture size [1/m]
        amrex::ParticleReal m_inv_aperture_y = 0; //! inverse of the vertical aperture size [1/m]
    };

} // namespace impactx::elements::mixin

#endif // IMPACTX_ELEMENTS_MIXIN_PIPE_APERTURE_H
