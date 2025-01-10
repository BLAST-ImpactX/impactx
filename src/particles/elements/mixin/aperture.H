/* Copyright 2022-2025 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_MIXIN_APERTURE_H
#define IMPACTX_ELEMENTS_MIXIN_APERTURE_H

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
#include <AMReX_Particle.H>
#include <AMReX_REAL.H>


namespace impactx::elements
{
    /** This is a helper class for applying a transverse aperture restriction to thick lattice elements
     */
    struct Aperture
    {

        /** A finite-length element
         *
         * @param x_aperture horizontal half-aperture size in m
         * @param y_aperture vertical half-aperture size in m
         */
        Aperture (
            amrex::ParticleReal x_aperture,
            amrex::ParticleReal y_aperture
        )
        : m_x_aperture(x_aperture), m_y_aperture(y_aperture)
        {
        }

        Aperture () = default;
        Aperture (Aperture const &) = default;
        Aperture& operator= (Aperture const &) = default;
        Aperture (Aperture&&) = default;
        Aperture& operator= (Aperture&& rhs) = default;

        ~Aperture () = default;

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

            // skip aperture application if x_aperture <= 0 or y_aperture <= 0
            if (m_x_aperture > 0 && m_y_aperture > 0) {

               // scale horizontal and vertical coordinates
               amrex::ParticleReal const u = x / m_x_aperture;
               amrex::ParticleReal const v = y / m_y_aperture;

               // compare against the aperture boundary
               if (std::pow(u,2) + std::pow(v,2) > 1_prt) {
                  amrex::ParticleIDWrapper{idcpu}.make_invalid();
               }

            }
        }

        /** Horizontal aperture size
         *
         * @return horizontal aperture size in m
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal x_aperture () const
        {
            return m_x_aperture;
        }

        /** Vertical aperture size
         *
         * @return vertical aperture size in m
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal y_aperture () const
        {
            return m_y_aperture;
        }

        amrex::ParticleReal m_x_aperture = 0; //! horizontal aperture size [m]
        amrex::ParticleReal m_y_aperture = 0; //! vertical aperture size [m]
    };

} // namespace impactx::elements

#endif // IMPACTX_ELEMENTS_MIXIN_APERTURE_H
