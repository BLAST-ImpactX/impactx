/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_AMR_CORE_DATA_H
#define IMPACTX_AMR_CORE_DATA_H

#include "AmrCoreData_fwd.H"
#include "particles/ImpactXParticleContainer.H"
#include "particles/CovarianceMatrix.H"

#include <AMReX_AmrCore.H>
#include <AMReX_AmrMesh.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_TagBox.H>

#include <optional>
#include <string>
#include <type_traits>
#include <unordered_map>


namespace impactx::initialization
{
    /** Runtime precision tag for the beam particles / fields. */
    enum class Precision
    {
        single,
        double_
    };

    /** Precision-agnostic base of the AmrCore data.
     *
     * This holds everything that does not depend on the runtime-selected beam
     * precision: the AMR mesh hierarchy (amrex::AmrCore, always double geometry),
     * the envelope/reference tracking data, and a precision-agnostic accessor to
     * the beam particle container. The concrete, precision-parametric data lives
     * in AmrCoreData<T_PT>.
     *
     * This wrapper is also needed because AmrCore cannot be constructed directly,
     * due to pure virtuals in it.
     */
    class AmrCoreDataBase
        : public amrex::AmrCore
    {
      public:
        AmrCoreDataBase (
            amrex::Geometry const& level_0_geom,
            amrex::AmrInfo const& amr_info
        );
        AmrCoreDataBase (
            amrex::RealBox const & rb,
            int max_level_in,
            amrex::Vector<int> const & n_cell_in,
            int coord,
            amrex::Vector<amrex::IntVect> const & ref_ratios,
            amrex::Array<int,AMREX_SPACEDIM> const & is_per
        );

        AmrCoreDataBase (AmrCoreDataBase&& rhs) = default;
        AmrCoreDataBase& operator= (AmrCoreDataBase&& rhs) = default;

        AmrCoreDataBase (const AmrCoreDataBase& rhs) = delete;
        AmrCoreDataBase& operator= (const AmrCoreDataBase& rhs) = delete;

        ~AmrCoreDataBase() override = default;

        /** The runtime precision of the beam particles in this data. */
        [[nodiscard]] virtual Precision precision () const noexcept = 0;

        /** Precision-agnostic access to the beam particle container. */
        [[nodiscard]] virtual ImpactXParticleContainerBase & particle_container () = 0;
        [[nodiscard]] virtual ImpactXParticleContainerBase const & particle_container () const = 0;

        /** Precision-agnostic access to the lost-particle container. */
        [[nodiscard]] virtual ImpactXParticleContainerBase & lost_particle_container () = 0;
        [[nodiscard]] virtual ImpactXParticleContainerBase const & lost_particle_container () const = 0;

        /** Precision-agnostic access to the space-charge mesh fields.
         *
         * These arrays are always double precision, independent of the beam
         * storage precision, so they are exposed on the agnostic base.
         */
        [[nodiscard]] virtual std::unordered_map<int, amrex::MultiFab> & rho () = 0;
        [[nodiscard]] virtual std::unordered_map<int, amrex::MultiFab> & phi () = 0;
        [[nodiscard]] virtual std::unordered_map<int, std::unordered_map<std::string, amrex::MultiFab> > & space_charge_field () = 0;

        /** Data needed for the beam envelop (covariance matrix) tracking mode */
        struct TrackModeEnvelope
        {
            /** The reference particle */
            std::optional<RefPart> m_ref;

            /** The 6x6 covariance matrix for the beam envelope, relative to the reference particle */
            std::optional<Envelope> m_env;

        } track_envelope;

        /** Data needed for the reference orbit tracking mode */
        struct TrackModeReference
        {
            /** The reference particle */
            std::optional<RefPart> m_ref;

        } track_reference;

        // precision-independent AMR callbacks
        void ErrorEst (
            [[maybe_unused]] int lev,
            [[maybe_unused]] amrex::TagBoxArray& tags,
            [[maybe_unused]] amrex::Real time,
            [[maybe_unused]] int ngrow) override;

        void MakeNewLevelFromCoarse (
            [[maybe_unused]] int lev,
            [[maybe_unused]] amrex::Real time,
            [[maybe_unused]] const amrex::BoxArray& ba,
            [[maybe_unused]] const amrex::DistributionMapping& dm) override;

        void RemakeLevel (
            [[maybe_unused]] int lev,
            [[maybe_unused]] amrex::Real time,
            [[maybe_unused]] const amrex::BoxArray& ba,
            [[maybe_unused]] const amrex::DistributionMapping& dm) override;

        // MakeNewLevelFromScratch and ClearLevel create/erase the precision-
        // dependent mesh fields and are implemented in AmrCoreData<T_PT>.
    };

    /** Data in AmrCore at a concrete beam precision T_PT.
     *
     * Adds the precision-dependent beam particle container (and lost particles)
     * to AmrCoreDataBase. The space-charge mesh fields are kept in double
     * precision for now (the particle<->mesh deposit/gather bridges the
     * precisions); making the field arrays single precision is a follow-up.
     *
     * @tparam T_PT the beam particle (SoA storage) precision, e.g. float or double
     */
    template <class T_PT=amrex::ParticleReal>
    class AmrCoreData final
        : public AmrCoreDataBase
    {
      public:
        using AmrCoreDataBase::AmrCoreDataBase;

        ~AmrCoreData() override = default;

        [[nodiscard]] Precision precision () const noexcept override
        {
            return std::is_same_v<T_PT, float> ? Precision::single : Precision::double_;
        }

        [[nodiscard]] ImpactXParticleContainerBase & particle_container () override
        {
            return *track_particles.m_particle_container;
        }
        [[nodiscard]] ImpactXParticleContainerBase const & particle_container () const override
        {
            return *track_particles.m_particle_container;
        }

        [[nodiscard]] ImpactXParticleContainerBase & lost_particle_container () override
        {
            return *track_particles.m_particles_lost;
        }
        [[nodiscard]] ImpactXParticleContainerBase const & lost_particle_container () const override
        {
            return *track_particles.m_particles_lost;
        }

        [[nodiscard]] std::unordered_map<int, amrex::MultiFab> & rho () override
        {
            return track_particles.m_rho;
        }
        [[nodiscard]] std::unordered_map<int, amrex::MultiFab> & phi () override
        {
            return track_particles.m_phi;
        }
        [[nodiscard]] std::unordered_map<int, std::unordered_map<std::string, amrex::MultiFab> > & space_charge_field () override
        {
            return track_particles.m_space_charge_field;
        }

        /** Data needed for the symplectic particle tracking mode */
        struct TrackModeParticles
        {
            /** these are the physical/beam particles of the simulation */
            std::unique_ptr<impactx::ImpactXParticleContainerT<T_PT>> m_particle_container;

            /** former beam particles that got lost in apertures, the wall, etc. */
            std::unique_ptr<impactx::ImpactXParticleContainerT<T_PT>> m_particles_lost;

            /** charge density per level */
            std::unordered_map<int, amrex::MultiFab> m_rho;
            /** scalar potential per level */
            std::unordered_map<int, amrex::MultiFab> m_phi;
            /** space charge field (vector) per level */
            std::unordered_map<int, std::unordered_map<std::string, amrex::MultiFab> > m_space_charge_field;

        } track_particles;

        void MakeNewLevelFromScratch (
            int lev,
            amrex::Real time,
            const amrex::BoxArray& ba,
            const amrex::DistributionMapping& dm) override;

        void ClearLevel (int lev) override;
    };

    /** Dispatch on the runtime beam precision and invoke a callable with the
     * concrete AmrCoreData<T_PT>.
     *
     * This is the single place where the runtime precision selection crosses
     * from agnostic code into compile-time templated code. The callable @p f is
     * invoked with a reference to the concrete AmrCoreData<T_PT> that matches
     * data.precision().
     *
     * @param data the precision-agnostic AmrCore data
     * @param f a generic callable, e.g. a generic lambda [&](auto& d){ ... }
     */
    template <class F>
    decltype(auto)
    visit_amr_data ([[maybe_unused]] AmrCoreDataBase & data, F && f)
    {
        // Dispatch to the concrete AmrCoreData<T_PT> for the runtime precision.
        // Only the compiled precisions have a case (and an instantiation); a
        // single-precision build collapses this to one branch.
#if defined(IMPACTX_COMPILE_SINGLE) && defined(IMPACTX_COMPILE_DOUBLE)
        if (data.precision() == Precision::single)
        {
            return f(static_cast<AmrCoreData<float>&>(data));
        }
        return f(static_cast<AmrCoreData<double>&>(data));
#elif defined(IMPACTX_COMPILE_SINGLE)
        return f(static_cast<AmrCoreData<float>&>(data));
#else
        return f(static_cast<AmrCoreData<double>&>(data));
#endif
    }

} // namespace impactx::initialization

#endif // IMPACTX_AMR_CORE_DATA_H
