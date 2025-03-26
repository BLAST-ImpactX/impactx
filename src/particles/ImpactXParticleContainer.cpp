/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Remi Lehe
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactXParticleContainer.H"

#include "initialization/AmrCoreData.H"

#include <ablastr/constant.H>
#include <ablastr/particles/ParticleMoments.H>
#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particle.H>

#include <algorithm>
#include <iterator>
#include <stdexcept>


namespace
{
    bool do_omp_dynamic ()
    {
        bool do_dynamic = true;
        amrex::ParmParse const pp_impactx("impactx");
        pp_impactx.query("do_dynamic_scheduling", do_dynamic);
        return do_dynamic;
    }
}

namespace impactx
{
    ParIterSoA::ParIterSoA (ContainerType& pc, int level)
        : amrex::ParIterSoA<RealSoA::nattribs, IntSoA::nattribs>(pc, level,
                   amrex::MFItInfo().SetDynamic(do_omp_dynamic())) {}

    ParIterSoA::ParIterSoA (ContainerType& pc, int level, amrex::MFItInfo& info)
        : amrex::ParIterSoA<RealSoA::nattribs, IntSoA::nattribs>(pc, level,
              info.SetDynamic(do_omp_dynamic())) {}

    ParConstIterSoA::ParConstIterSoA (ContainerType& pc, int level)
        : amrex::ParConstIterSoA<RealSoA::nattribs, IntSoA::nattribs>(pc, level,
              amrex::MFItInfo().SetDynamic(do_omp_dynamic())) {}

    ParConstIterSoA::ParConstIterSoA (ContainerType& pc, int level, amrex::MFItInfo& info)
        : amrex::ParConstIterSoA<RealSoA::nattribs, IntSoA::nattribs>(pc, level,
              info.SetDynamic(do_omp_dynamic())) {}

    ImpactXParticleContainer::ImpactXParticleContainer (initialization::AmrCoreData* amr_core)
        : amrex::ParticleContainerPureSoA<RealSoA::nattribs, IntSoA::nattribs>(amr_core->GetParGDB())
    {
        SetParticleSize();
        SetSoACompileTimeNames(
            {RealSoA::names_s.begin(), RealSoA::names_s.end()},
            {IntSoA::names_s.begin(), IntSoA::names_s.end()}
        );
    }

    void
    ImpactXParticleContainer::SetLostParticleContainer (ImpactXParticleContainer * lost_pc)
    {
        m_particles_lost = lost_pc;
    }

    ImpactXParticleContainer *
    ImpactXParticleContainer::GetLostParticleContainer ()
    {
        if (m_particles_lost == nullptr)
        {
            throw std::logic_error(
                    "ImpactXParticleContainer::GetLostParticleContainer No lost particle container registered yet.");
        } else {
            return m_particles_lost;
        }
    }

    void ImpactXParticleContainer::SetParticleShape (int order) {
        if (m_particle_shape.has_value())
        {
            throw std::logic_error(
                "ImpactXParticleContainer::SetParticleShape This was already called before and cannot be changed.");
        } else
        {
            if (order < 1 || order > 3) {
                amrex::Abort("algo.particle_shape order can be only 1, 2, or 3");
            }
            m_particle_shape = order;
        }
    }

    void ImpactXParticleContainer::SetParticleShape ()
    {
        amrex::ParmParse const pp_algo("algo");
        int v = 0;
        bool const has_shape = pp_algo.queryWithParser("particle_shape", v);
        if (!has_shape)
            throw std::runtime_error("particle_shape is not set, cannot initialize grids with guard cells.");
        SetParticleShape(v);
    }

    void
    ImpactXParticleContainer::prepare ()
    {
        // make sure we have at least one box with enough tiles for each OpenMP thread

        // make sure level 0, grid 0 exists
        int lid = 0, gid = 0;
        {
            const auto& pmap = ParticleDistributionMap(lid).ProcessorMap();
            auto it = std::find(pmap.begin(), pmap.end(), amrex::ParallelDescriptor::MyProc());
            if (it == std::end(pmap)) {
                amrex::Abort("Particle container needs to have at least one grid.");
            } else {
                gid = *it;
            }
        }

        int nthreads = 1;
#if defined(AMREX_USE_OMP)
        nthreads = omp_get_max_threads();
#endif

        // When running without space charge and OpenMP parallelization,
        // we need to make sure we have enough tiles on level 0, grid 0
        // to thread over the available tiles. We try to set the tile_size
        // appropriately here. The tiles start as (very large number, 8, 8)
        // in (x, y, z). So we try to reduce the tile size in the y and z
        // directions alternately until there are enough tiles for the number of threads.

        const auto& ba = ParticleBoxArray(lid);
        auto n_logical = numTilesInBox(ba[gid], true, tile_size);

        int ntry = 0;
        constexpr int max_tries = 6;
        while ((n_logical < nthreads) && (ntry++ < max_tries)) {
            int idim = (ntry % 2) + 1;  // alternate between 1 and 2
            tile_size[idim] /= 2;
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(tile_size[idim] > 0,
                                             "Failed to set proper tile size for the number of OpenMP threads. "
                                             "Consider lowering the number of OpenMP threads via the environment variable OMP_NUM_THREADS."
                                             );
            n_logical = numTilesInBox(ba[gid], true, tile_size);
        }

        if (n_logical < nthreads) {
            amrex::Abort("ImpactParticleContainer::prepare() "
                "could not find good tile size for the number of OpenMP threads. "
                "Consider lowering the number of OpenMP threads via the environment variable OMP_NUM_THREADS."
            );
        }

        reserveData();
        resizeData();
    }

    void
    ImpactXParticleContainer::AddNParticles (
        amrex::Gpu::DeviceVector<amrex::ParticleReal> const & x,
        amrex::Gpu::DeviceVector<amrex::ParticleReal> const & y,
        amrex::Gpu::DeviceVector<amrex::ParticleReal> const & t,
        amrex::Gpu::DeviceVector<amrex::ParticleReal> const & px,
        amrex::Gpu::DeviceVector<amrex::ParticleReal> const & py,
        amrex::Gpu::DeviceVector<amrex::ParticleReal> const & pt,
        amrex::ParticleReal qm,
        amrex::ParticleReal bchchg
    )
    {
        BL_PROFILE("ImpactX::AddNParticles");

        AMREX_ALWAYS_ASSERT(x.size() == y.size());
        AMREX_ALWAYS_ASSERT(x.size() == t.size());
        AMREX_ALWAYS_ASSERT(x.size() == px.size());
        AMREX_ALWAYS_ASSERT(x.size() == py.size());
        AMREX_ALWAYS_ASSERT(x.size() == pt.size());

        // number of particles to add
        int const np = x.size();

        // we add particles to lev 0, grid 0
        int lid = 0, gid = 0;
        {
            const auto& pmap = ParticleDistributionMap(lid).ProcessorMap();
            auto it = std::find(pmap.begin(), pmap.end(), amrex::ParallelDescriptor::MyProc());
            if (it == std::end(pmap)) {
                amrex::Abort("Attempting to add particles to box that does not exist.");
            } else {
                gid = *it;
            }
        }

        int nthreads = 1;
#if defined(AMREX_USE_OMP)
        nthreads = omp_get_max_threads();
#endif

        // split up particles over nthreads tiles
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(numTilesInBox(ParticleBoxArray(lid)[gid], true, tile_size) >= nthreads, "Not enough tiles for the number of OpenMP threads - please try reducing particles.tile_size or increasing the number of cells in the domain.");
        for (int ithr = 0; ithr < nthreads; ++ithr) {
            DefineAndReturnParticleTile(lid, gid, ithr);
        }

        int pid = ParticleType::NextID();
        ParticleType::NextID(pid + np);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            static_cast<amrex::Long>(pid) + static_cast<amrex::Long>(np) < amrex::LongParticleIds::LastParticleID,
            "ERROR: overflow on particle id numbers"
        );

#if defined(AMREX_USE_OMP)
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        {
            int tid = 0;
#if defined(AMREX_USE_OMP)
            tid = omp_get_thread_num();
#endif

            // we split up the np particles onto multiple tiles.
            // if nthreads evenly divides np, each thread will
            // get get n_regular particles. If there are some
            // leftovers, however, the first n_leftover threads
            // will get one extra
            int n_regular  = np / nthreads;
            int n_leftover = np - n_regular*nthreads;

            int num_to_add = 0; /* how many particles this tile will get*/
            int my_offset = 0; /*  offset into global arrays x, y, etc. for this thread*/

            if (tid < n_leftover) { // get n_regular+1 items
                my_offset = tid * (n_regular + 1);
                num_to_add = n_regular + 1;
            } else {         // get n_regular items
                my_offset = tid * n_regular + n_leftover;
                num_to_add = n_regular;
            }

            auto& particle_tile = ParticlesAt(lid, gid, tid);
            auto old_np = particle_tile.numParticles();
            auto new_np = old_np + num_to_add;
            particle_tile.resize(new_np);

            const int cpuid = amrex::ParallelDescriptor::MyProc();

            auto & soa = particle_tile.GetStructOfArrays().GetRealData();
            amrex::ParticleReal * const AMREX_RESTRICT x_arr = soa[RealSoA::x].dataPtr();
            amrex::ParticleReal * const AMREX_RESTRICT y_arr = soa[RealSoA::y].dataPtr();
            amrex::ParticleReal * const AMREX_RESTRICT t_arr = soa[RealSoA::t].dataPtr();
            amrex::ParticleReal * const AMREX_RESTRICT px_arr = soa[RealSoA::px].dataPtr();
            amrex::ParticleReal * const AMREX_RESTRICT py_arr = soa[RealSoA::py].dataPtr();
            amrex::ParticleReal * const AMREX_RESTRICT pt_arr = soa[RealSoA::pt].dataPtr();
            amrex::ParticleReal * const AMREX_RESTRICT qm_arr = soa[RealSoA::qm].dataPtr();
            amrex::ParticleReal * const AMREX_RESTRICT w_arr  = soa[RealSoA::w ].dataPtr();

            uint64_t * const AMREX_RESTRICT idcpu_arr = particle_tile.GetStructOfArrays().GetIdCPUData().dataPtr();

            amrex::ParticleReal const * const AMREX_RESTRICT x_ptr = x.data();
            amrex::ParticleReal const * const AMREX_RESTRICT y_ptr = y.data();
            amrex::ParticleReal const * const AMREX_RESTRICT t_ptr = t.data();
            amrex::ParticleReal const * const AMREX_RESTRICT px_ptr = px.data();
            amrex::ParticleReal const * const AMREX_RESTRICT py_ptr = py.data();
            amrex::ParticleReal const * const AMREX_RESTRICT pt_ptr = pt.data();

            amrex::ParallelFor(num_to_add,
                [=] AMREX_GPU_DEVICE (int i) noexcept
            {
                idcpu_arr[old_np+i] = amrex::SetParticleIDandCPU(pid + my_offset + i, cpuid);

                x_arr[old_np+i] = x_ptr[my_offset+i];
                y_arr[old_np+i] = y_ptr[my_offset+i];
                t_arr[old_np+i] = t_ptr[my_offset+i];

                px_arr[old_np+i] = px_ptr[my_offset+i];
                py_arr[old_np+i] = py_ptr[my_offset+i];
                pt_arr[old_np+i] = pt_ptr[my_offset+i];

                qm_arr[old_np+i] = qm;
                w_arr[old_np+i]  = bchchg/ablastr::constant::SI::q_e/np;
            });
        }

        // safety first: in case passed attribute arrays were temporary, we
        // want to make sure the ParallelFor has ended here
        amrex::Gpu::streamSynchronize();
    }

    void
    ImpactXParticleContainer::SetRefParticle (RefPart const & refpart)
    {
        m_refpart = refpart;
    }

    RefPart &
    ImpactXParticleContainer::GetRefParticle ()
    {
        return m_refpart;
    }

    RefPart const &
    ImpactXParticleContainer::GetRefParticle () const
    {
        return m_refpart;
    }

    void
    ImpactXParticleContainer::SetRefParticleEdge ()
    {
        m_refpart.sedge = m_refpart.s;
    }

    std::tuple<
            amrex::ParticleReal, amrex::ParticleReal,
            amrex::ParticleReal, amrex::ParticleReal,
            amrex::ParticleReal, amrex::ParticleReal>
    ImpactXParticleContainer::MinAndMaxPositions ()
    {
        BL_PROFILE("ImpactXParticleContainer::MinAndMaxPositions");
        return ablastr::particles::MinAndMaxPositions(*this);
    }

    std::tuple<
            amrex::ParticleReal, amrex::ParticleReal,
            amrex::ParticleReal, amrex::ParticleReal,
            amrex::ParticleReal, amrex::ParticleReal>
    ImpactXParticleContainer::MeanAndStdPositions ()
    {
        BL_PROFILE("ImpactXParticleContainer::MeanAndStdPositions");
        return ablastr::particles::MeanAndStdPositions<
            ImpactXParticleContainer, RealSoA::w
        >(*this);
    }

    CoordSystem
    ImpactXParticleContainer::GetCoordSystem () const
    {
        return m_coordsystem;
    }

    void
    ImpactXParticleContainer::SetCoordSystem (CoordSystem coord_system)
    {
        m_coordsystem = coord_system;
    }
} // namespace impactx
