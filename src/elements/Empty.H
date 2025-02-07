/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENT_EMPTY_H
#define IMPACTX_ELEMENT_EMPTY_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/thin.H"
#include "mixin/lineartransport.H"
#include "mixin/named.H"
#include "mixin/nofinalize.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>


namespace impactx::elements
{
    struct Empty
    : public mixin::Thin,
      public mixin::LinearTransport<Empty>,
      public mixin::NoFinalize
    {
        static constexpr auto type = "None";
        using PType = ImpactXParticleContainer::ParticleType;

        /** This element does nothing.
         */
        Empty ()
        {
        }

        /** Push all particles - nothing to do here */
        void operator() (
            ImpactXParticleContainer & /* pc */,
            int /* step */,
            int /* period */
        ) {
            // nothing to do
        }

        /** Push all particles - nothing to do here */
        void operator() (
            ImpactXParticleContainer::iterator & /* pti */,
            RefPart & AMREX_RESTRICT /* ref_part */
        ) {
            // nothing to do
        }

        /** Dos nothing to the reference particle.
         *
         * @param[in,out] refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() ([[maybe_unused]] RefPart & AMREX_RESTRICT refpart) const
        {
            // nothing to do
        }

        /** This function returns the linear transport map.
         *
         * @param[in] refpart reference particle
         * @returns 6x6 transport matrix
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        Map6x6
        transport_map ([[maybe_unused]] RefPart const & AMREX_RESTRICT refpart) const
        {
            // nothing to do
            return Map6x6::Identity();
        }

        /** Does nothing to a particle.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param idcpu particle global index (unused)
         * @param refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            [[maybe_unused]] amrex::ParticleReal & AMREX_RESTRICT x,
            [[maybe_unused]] amrex::ParticleReal & AMREX_RESTRICT y,
            [[maybe_unused]] amrex::ParticleReal & AMREX_RESTRICT t,
            [[maybe_unused]] amrex::ParticleReal & AMREX_RESTRICT px,
            [[maybe_unused]] amrex::ParticleReal & AMREX_RESTRICT py,
            [[maybe_unused]] amrex::ParticleReal & AMREX_RESTRICT pt,
            [[maybe_unused]] uint64_t & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & refpart
        ) const
        {
            // nothing to do
        }

        /** This pushes the reference particle. */
        using Thin::operator();

        /** This pushes the covariance matrix. */
        using LinearTransport::operator();
    };

} // namespace impactx

#endif // IMPACTX_ELEMENT_EMPTY_H
