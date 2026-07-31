/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_SOURCE_H
#define IMPACTX_SOURCE_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/lineartransport.H"
#include "mixin/named.H"
#include "mixin/nofinalize.H"
#include "mixin/thin.H"

#include <stdexcept>
#include <string>

namespace impactx::elements
{
    struct Source
    : public mixin::Named,
      public mixin::LinearTransport<Source>,
      public mixin::Thin,
      public mixin::NoFinalize
    {
        static constexpr auto type = "Source";
        using PType = ImpactXParticleContainer::ParticleType;

        static constexpr bool DEFAULT_active_once = true;
        static constexpr bool DEFAULT_load_ref_particle = true;

        /** A particle source
         *
         * @param distribution Must read "openPMD" for this constructor
         * @param openpmd_path path to openPMD series as accepted by openPMD::Series
         * @param active_once inject particles only for the first lattice period
         * @param load_ref_particle restore the reference particle from the species metadata
         *        of the openPMD file (particle tracking only)
         * @param name a user defined and not necessarily unique name of the element
         */
        Source (
            std::string distribution,
            std::string openpmd_path,
            bool active_once = DEFAULT_active_once,
            bool load_ref_particle = DEFAULT_load_ref_particle,
            std::optional<std::string> name = DEFAULT_name
        )
        : Named(std::move(name)),
          m_distribution(distribution),
          m_series_name(std::move(openpmd_path)),
          m_active_once(active_once),
          m_load_ref_particle(load_ref_particle)
        {
            if (distribution != "openPMD") {
                throw std::runtime_error("Only 'openPMD' distribution is supported if openpmd_path is provided!");
            }
        }

        /** Reverse the element (not supported for Source) */
        void reverse ()
        {
            throw std::runtime_error("Source: reverse() is not supported.");
        }

        /** Initialize/Load particles.
         *
         * Particles are relative to the reference particle.
         * If load_ref_particle is set, the reference particle is restored from
         * the species metadata of the file before the particles are added.
         *
         * @param[in,out] pc particle container to push
         * @param[in] step global step for diagnostics
         * @param[in] period for periodic lattices, this is the current period (turn or cycle)
         */
        void operator() (
            ImpactXParticleContainer & pc,
            [[maybe_unused]] int step,
            int period
        );

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

        /** This does nothing to the reference particle. */
        using Thin::operator();

        /** This pushes the covariance matrix. */
        using LinearTransport::operator();

        std::string m_distribution; //! Distribution type of particles in the source
        std::string m_series_name; //! openPMD filename
        bool m_active_once = true; //! Inject particles only for the first lattice period
        bool m_load_ref_particle = true; //! Restore the reference particle from the file's species metadata
    };

} // namespace impactx

#endif  // IMPACTX_SOURCE_H
