/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_MIXIN_THICK_H
#define IMPACTX_ELEMENTS_MIXIN_THICK_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>


namespace impactx::elements::mixin
{
    /** This is a helper class for lattice elements with finite length
     */
    struct Thick
    {
        /** Whether the element's external-field transport map is an exact constant-field,
         *  one-parameter flow in the slice length, i.e. ``M(ds/2) M(ds/2) == M(ds)`` exactly
         *  and independent of how the element is sliced (closed-form map, constant field,
         *  no acceleration).
         *
         *  Only such elements are eligible for the cheaper drift-outer (MKM) Strang split
         *  of the collective effects (\see impactx::track_lattice_particles), which realizes
         *  the half-slice transports by temporarily doubling ``nslice``. Elements whose map
         *  is s-dependent / numerically integrated (soft, ``Exact{Quad,Multipole,CFbend}``,
         *  ``RFCavity``) or accelerating (``ChrAcc``) must keep this ``false``, as doubling
         *  ``nslice`` would refine their field integration. Eligible elements opt in by
         *  shadowing this with ``static constexpr bool is_constant_field_flow = true;``.
         */
        static constexpr bool is_constant_field_flow = false;

        /** A finite-length element
         *
         * @param ds Segment length in m
         * @param nslice number of slices used for the application of space charge
         */
        Thick (
            amrex::ParticleReal ds,
            int nslice
        )
        : m_ds(ds), m_nslice(nslice)
        {
        }

        /** Number of slices used for the application of space charge
         *
         * @return positive integer
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        int nslice () const
        {
            return m_nslice;
        }

        /** Return the segment length
         *
         * @return value in meters
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal ds () const
        {
            return m_ds;
        }

        /** Reverse the element (negate ds)
         *
         * After calling this method, pushing particles through the element
         * reverses the effect of the original element.
         */
        void reverse ()
        {
            m_ds = -m_ds;
        }

        amrex::ParticleReal m_ds; //! segment length in m
        int m_nslice; //! number of slices used for the application of space charge
    };

} // namespace impactx::elements::mixin

#endif // IMPACTX_ELEMENTS_MIXIN_THICK_H
