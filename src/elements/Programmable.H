/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_PROGRAMMABLE_H
#define IMPACTX_ELEMENTS_PROGRAMMABLE_H

#include "particles/ImpactXParticleContainer.H"
#include "particles/CovarianceMatrix.H"
#include "mixin/lineartransport.H"
#include "mixin/named.H"
#include "mixin/thick.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>
// TODO: #include <AMReX_SIMD.H>
//   SIMD will require better complex math support in vir-simd / C++26
//   https://github.com/mattkretz/vir-simd/issues/42

#include <functional>
#include <stdexcept>


namespace impactx::elements
{
    struct Programmable
    : public mixin::Named,
      public mixin::LinearTransport<Programmable, false>
      // TODO: public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "Programmable";
        using PType = ImpactXParticleContainer::ParticleType;

        /** This element can be programmed
         *
         * @param ds Segment length in m
         * @param nslice number of slices used for the application of space charge
         * @param name a user defined and not necessarily unique name of the element
         */
        Programmable (
            amrex::ParticleReal ds=0.0,
            int nslice=1,
            std::optional<std::string> name = std::nullopt
        )
          : Named(std::move(name)), m_ds(ds), m_nslice(nslice)
        {}

        /** Reverse the element (not supported for Programmable) */
        void reverse ()
        {
            throw std::runtime_error("Programmable: reverse() is not supported.");
        }

        /** Push all particles relative to the reference particle
         *
         * @param[in,out] pc particle container to push
         * @param[in] step global step for diagnostics
         * @param[in] period for periodic lattices, this is the current period (turn or cycle)
         */
        void operator() (
            ImpactXParticleContainer & pc,
            int step,
            int period
        ) const;

        /** Push all particles relative to the reference particle */
        void operator() (
            ImpactXParticleContainer::iterator & pti,
            RefPart & ref_part,
            [[maybe_unused]] bool spin
        ) const;

        /** This pushes the reference particle.
         *
         * @param[in,out] refpart reference particle
         */
        void operator() (
            RefPart & refpart
        ) const;

        /** This pushes the covariance matrix. */
        using LinearTransport::operator();

        /** This function returns the linear transport map.
         *
         * By default, @c Programmable provides no linear map because its
         * particle push is an arbitrary user callback. Users who know the
         * Jacobian of their push around the reference trajectory can opt
         * in by setting the @c m_transport_map hook: a callable that
         * returns a 6x6 transport matrix as a function of the current
         * reference particle. When set, this hook is used by envelope
         * tracking, @c lattice.transfer_map, @c map_trace, and Twiss.
         *
         * If the hook is not set, this function throws a runtime error
         * pointing the user to the hook (or, for purely linear pushes,
         * to the @c LinearMap element).
         *
         * @param[in] refpart reference particle
         * @returns 6x6 transport matrix from the user-supplied hook
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        Map6x6
        transport_map (RefPart const & AMREX_RESTRICT refpart) const
        {
            if (m_transport_map == nullptr) {
                throw std::runtime_error(
                    "Programmable: no linear transport map hook was set. "
                    "Set Programmable.transport_map_hook = <callable(refpart)->6x6> "
                    "for envelope/linear-optics tracking, or use the LinearMap "
                    "element for a purely linear map."
                );
            }
            return m_transport_map(refpart);
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

        /** Close and deallocate all data and handles.
         */
        void
        finalize ();

        amrex::ParticleReal m_ds = 0.0; //! segment length in m
        int m_nslice = 1; //! number of slices used for the application of space charge

        /** Allow threading via OpenMP for the particle iterator loop
         *
         * This will only affect threading if the OMP backend is active.
         *
         * The default value is false to do the safe thing by default. Users
         * must opt-in their guarantee that their code is thread-safe.
         */
        bool m_threadsafe = false;

        std::function<void(ImpactXParticleContainer *, int, int)> m_push; //! hook for push of whole container (pc, step, period)
        std::function<void(ImpactXParticleContainer::iterator *, RefPart &)> m_beam_particles; //! hook for beam particles (pti, ref)
        std::function<void(RefPart &)> m_ref_particle; //! hook for reference particle
        std::function<Map6x6(RefPart const &)> m_transport_map; //! optional hook for linear transport map (refpart -> 6x6)
        std::function<void()> m_finalize; //! hook for finalize cleanup
    };

} // namespace impactx

#endif // IMPACTX_ELEMENTS_PROGRAMMABLE_H
