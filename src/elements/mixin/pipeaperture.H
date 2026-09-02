/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
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
#include <AMReX_SIMD.H>


namespace impactx::elements::mixin
{
    /** This is a helper class for applying a transverse aperture restriction to thick lattice elements
     */
    struct PipeAperture
    {
        static constexpr amrex::ParticleReal DEFAULT_aperture_x = 0;
        static constexpr amrex::ParticleReal DEFAULT_aperture_y = 0;

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
            set_aperture_x(aperture_x);
            set_aperture_y(aperture_y);
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
        template<typename T_Real=amrex::ParticleReal, typename T_IdCpu=uint64_t>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void apply_aperture (
            T_Real & AMREX_RESTRICT x,
            T_Real & AMREX_RESTRICT y,
            T_IdCpu & AMREX_RESTRICT idcpu
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;
            namespace stdx = amrex::simd::stdx;

            // skip aperture application if aperture_x <= 0 or aperture_y <= 0
            if (aperture_active()) {

                // scale horizontal and vertical coordinates
                // performance optimization: we intentionally store the inverse of
                //   the aperture, so we can do a quick multiplication (instead of a
                //   costly division) here.
                T_Real const u = x * m_inv_aperture_x;
                T_Real const v = y * m_inv_aperture_y;

                // compare against the aperture boundary
                auto const mask = u*u + v*v > 1_prt;
                amrex::ParticleIDWrapper<T_IdCpu>{idcpu}.make_invalid(mask);
            }
        }

        /** Is the aperture restriction active?
         *
         * A non-positive aperture in either plane removes the restriction,
         * @see apply_aperture.
         *
         * @return true if particles can be lost in this aperture
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool aperture_active () const
        {
            using namespace amrex::literals; // for _prt
            return m_inv_aperture_x > 0_prt && m_inv_aperture_y > 0_prt;
        }

        /** Horizontal aperture size
         *
         * @return horizontal aperture size [m]
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal aperture_x () const
        {
            using namespace amrex::literals; // for _prt
            return m_inv_aperture_x > 0_prt ? 1_prt / m_inv_aperture_x : 0_prt;
        }

        /** Vertical aperture size
         *
         * @return vertical aperture size [m]
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal aperture_y () const
        {
            using namespace amrex::literals; // for _prt
            return m_inv_aperture_y > 0_prt ? 1_prt / m_inv_aperture_y : 0_prt;
        }

        /** Set the horizontal aperture size
         *
         * A value of zero or less removes the horizontal aperture restriction.
         *
         * @param aperture_x horizontal aperture size [m]
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void set_aperture_x (amrex::ParticleReal aperture_x)
        {
            using namespace amrex::literals; // for _prt

            // performance optimization: we intentionally store the inverse of
            //   the aperture, @see m_inv_aperture_x
            m_inv_aperture_x = aperture_x > 0_prt ? 1_prt / aperture_x : 0_prt;
        }

        /** Set the vertical aperture size
         *
         * A value of zero or less removes the vertical aperture restriction.
         *
         * @param aperture_y vertical aperture size [m]
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void set_aperture_y (amrex::ParticleReal aperture_y)
        {
            using namespace amrex::literals; // for _prt

            // performance optimization: we intentionally store the inverse of
            //   the aperture, @see m_inv_aperture_y
            m_inv_aperture_y = aperture_y > 0_prt ? 1_prt / aperture_y : 0_prt;
        }

        // performance optimization: we intentionally store the inverse of
        //   the aperture, so we can do a quick multiplication instead of a costly
        //   division in the comparison operation per particle
        amrex::ParticleReal m_inv_aperture_x = 0; //! inverse of the horizontal aperture size [1/m]
        amrex::ParticleReal m_inv_aperture_y = 0; //! inverse of the vertical aperture size [1/m]
    };

} // namespace impactx::elements::mixin

#endif // IMPACTX_ELEMENTS_MIXIN_PIPE_APERTURE_H
