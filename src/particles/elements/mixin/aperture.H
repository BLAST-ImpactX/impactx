/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_MIXIN_APERTURE_H
#define IMPACTX_ELEMENTS_MIXIN_APERTURE_H

#include "particles/ImpactXParticleContainer.H"

#include <ablastr/constant.H>

#include <AMReX_Math.H>
#include <AMReX_Extension.H>
#include <AMReX_REAL.H>


namespace impactx::elements
{
    /** This is a helper class for applying a transverse aperture restriction to thick lattice elements
     */
    struct Aperture
    {

        /** A finite-length element
         *
         * @param xmax horizontal aperture size in m
         * @param ymax vertical aperture size in m
         */
        Aperture (
            amrex::ParticleReal xmax,
            amrex::ParticleReal ymax
        )
        : m_xmax(xmax), m_ymax(ymax)
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

            // scale horizontal and vertical coordinates
            amrex::ParticleReal const u = x / m_xmax;
            amrex::ParticleReal const v = y / m_ymax;

            // compare against the aperture boundary
            if (std::pow(u,2) + std::pow(v,2) > 1_prt) {
               amrex::ParticleIDWrapper{idcpu}.make_invalid();
            }

        }

        /** Horizontal aperture size
         *
         * @return horizontal aperture size in m
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal xmax () const
        {
            return m_xmax;
        }

        /** Vertical aperture size
         *
         * @return vertical aperture size in m
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal ymax () const
        {
            return m_ymax;
        }

        amrex::ParticleReal m_xmax = 0; //! horizontal aperture size [m]
        amrex::ParticleReal m_ymax = 0; //! vertical aperture size [m]
    };

} // namespace impactx::elements

#endif // IMPACTX_ELEMENTS_MIXIN_APERTURE_H
