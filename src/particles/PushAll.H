/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_PUSH_ALL_H
#define IMPACTX_PUSH_ALL_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_BLProfiler.H>


namespace impactx
{
    /** Push all particles in a particle container.
     *
     * This element pushes first the reference particle, then all other particles.
     * All particles are pushed independently with the same logic.
     * Particles are relative to the reference particle.
     *
     * @param[in,out] pc particle container to push
     * @param[in,out] element the beamline element
     * @param[in] step global step for diagnostics
     * @param[in] period for periodic lattices, this is the current period (turn or cycle)
     * @param[in] omp_parallel allow threading via OpenMP for the particle iterator loop (note: if OMP backend is active)
     */
    template<typename T_Element, class T_PC>
    void push_all (
        T_PC & pc,
        T_Element & element,
        [[maybe_unused]] int step,
        [[maybe_unused]] int period,
        [[maybe_unused]] bool omp_parallel = true
    )
    {
        // performance profiling per element
        std::string const profile_name = "impactx::push::" + std::string(T_Element::type);
        BL_PROFILE(profile_name);

        // is spin tracking enabled?
        bool spin = false;
        amrex::ParmParse pp_algo("algo");
        pp_algo.query("spin", spin);

        // preparing to access reference particle data: RefPart
        RefPart & ref_part = pc.GetRefParticle();

        // push reference particle in global coordinates
        {
            BL_PROFILE("impactx::push::RefPart");
            element(ref_part);
        }

        // loop over refinement levels
        int const nLevel = pc.finestLevel();
        for (int lev = 0; lev <= nLevel; ++lev)
        {
            // loop over all particle boxes
            using ParIt = typename T_PC::iterator;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion() && omp_parallel)
#endif
            for (ParIt pti(pc, lev); pti.isValid(); ++pti) {
                // push beam particles relative to reference particle
                element(pti, ref_part, spin);
            } // end loop over all particle boxes
        } // env mesh-refinement level loop
    }

} // namespace impactx

/** Suppress implicit instantiation of push_all<T> — use in element headers. */
#define IMPACTX_PUSH_EXTERN_TEMPLATE(ElementType) \
    extern template void impactx::push_all<ElementType>( \
        impactx::ImpactXParticleContainer &, \
        ElementType &, \
        int, int, bool);

/** Explicit instantiation of push_all<T> — use in element .cpp files. */
#define IMPACTX_PUSH_INSTANTIATE(ElementType) \
    template void impactx::push_all<ElementType>( \
        impactx::ImpactXParticleContainer &, \
        ElementType &, \
        int, int, bool);

#endif // IMPACTX_PUSH_ALL_H
