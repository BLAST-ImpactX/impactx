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
#include "initialization/Algorithms.H"
#include "diagnostics/ReducedBeamCharacteristics.H"
#include "particles/SplitEqually.H"

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
#include <sstream>


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
    template <class T_PT>
    ParIterSoAT<T_PT>::ParIterSoAT (ContainerType& pc, int level)
        : base_t(pc, level, amrex::MFItInfo().SetDynamic(do_omp_dynamic())) {}

    template <class T_PT>
    ParIterSoAT<T_PT>::ParIterSoAT (ContainerType& pc, int level, amrex::MFItInfo& info)
        : base_t(pc, level, info.SetDynamic(do_omp_dynamic())) {}

    template <class T_PT>
    ParConstIterSoAT<T_PT>::ParConstIterSoAT (ContainerType& pc, int level)
        : base_t(pc, level, amrex::MFItInfo().SetDynamic(do_omp_dynamic())) {}

    template <class T_PT>
    ParConstIterSoAT<T_PT>::ParConstIterSoAT (ContainerType& pc, int level, amrex::MFItInfo& info)
        : base_t(pc, level, info.SetDynamic(do_omp_dynamic())) {}

    template <class T_PT>
    ImpactXParticleContainerT<T_PT>::ImpactXParticleContainerT (initialization::AmrCoreData* amr_core)
        : base_t(amr_core->GetParGDB())
    {
        this->SetArena(amrex::The_Arena());
        this->SetParticleSize();
        this->SetSoACompileTimeNames(
            {RealSoA::names_s.begin(), RealSoA::names_s.end()},
            {IntSoA::names_s.begin(), IntSoA::names_s.end()}
        );
    }

    template <class T_PT>
    void
    ImpactXParticleContainerT<T_PT>::SetLostParticleContainer (ImpactXParticleContainerT<T_PT> * lost_pc)
    {
        m_particles_lost = lost_pc;
    }

    template <class T_PT>
    ImpactXParticleContainerT<T_PT> *
    ImpactXParticleContainerT<T_PT>::GetLostParticleContainer ()
    {
        if (m_particles_lost == nullptr)
        {
            throw std::logic_error(
                    "ImpactXParticleContainer::GetLostParticleContainer No lost particle container registered yet.");
        } else {
            return m_particles_lost;
        }
    }

    template <class T_PT>
    void ImpactXParticleContainerT<T_PT>::SetParticleShape (int order) {
        if (m_particle_shape.has_value())
        {
            throw std::logic_error(
                "ImpactXParticleContainer::SetParticleShape This was already called before and cannot be changed.");
        } else
        {
            if (order < 0 || order > 3) {
                amrex::Abort("algo.particle_shape order can be only 0, 1, 2, or 3");
            }
            m_particle_shape = order;
        }
    }

    template <class T_PT>
    void ImpactXParticleContainerT<T_PT>::SetParticleShape ()
    {
        amrex::ParmParse const pp_algo("algo");
        int v = 0;
        bool const has_shape = pp_algo.queryWithParser("particle_shape", v);
        if (!has_shape) {
            bool csr = false;
            pp_algo.query("csr", csr);
            auto space_charge = get_space_charge_algo();
            std::string track = "particles";
            pp_algo.query("track", track);
            if (csr ||
                (space_charge != SpaceChargeAlgo::False && track == "particles"))
            {
                throw std::runtime_error("particle_shape is not set, cannot initialize grids with guard cells for collective effects.");
            }
        }
        SetParticleShape(v);
    }

    template <class T_PT>
    void ImpactXParticleContainerT<T_PT>::SetBucketLength (amrex::ParticleReal bucket_length)
    {
        m_bucket_length = bucket_length;
    }

    template <class T_PT>
    amrex::ParticleReal
    ImpactXParticleContainerT<T_PT>::GetBucketLength ()
    {
        return m_bucket_length;
    }

    template <class T_PT>
    void
    ImpactXParticleContainerT<T_PT>::prepare ()
    {
        // make sure we have at least one box with enough tiles for each OpenMP thread

        // make sure level 0, grid 0 exists
        int lid = 0, gid = 0;
        {
            const auto& pmap = this->ParticleDistributionMap(lid).ProcessorMap();
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

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            nthreads == 1 || this->do_tiling == true,
            "OpenMP threads are used for parallelism but particles.do_tiling is set to false. "
            "This would create invalid particle data."
        );

        // When running without space charge and OpenMP parallelization,
        // we need to make sure we have enough tiles on level 0, grid 0
        // to thread over the available tiles. We try to set the tile_size
        // appropriately here. We start by setting the tile_size to the size of box 0.
        // Then we try to reduce the tile size in the (z, y, x) directions alternately
        // until there are enough tiles for the number of threads.

        // if the user has set particles.tile_size, do not override this choice
        // unlike particles.do_tiling, we do not add this to the table in init,
        // so if it's there the user set it
        bool user_set = false;
        {
            amrex::ParmParse pp_particles("particles");
            amrex::Vector<int> tilesize(AMREX_SPACEDIM);
            user_set = pp_particles.queryarr("tile_size", tilesize, 0, AMREX_SPACEDIM);
        }

        int n_logical = 0;
        const auto& ba = this->ParticleBoxArray(lid);
        if (!user_set)
        {
            this->tile_size = ba[gid].size();
            n_logical = amrex::numTilesInBox(ba[gid], true, this->tile_size);
            int ntry = 0;
            constexpr int max_tries = 10;
            while ((n_logical < nthreads) && (ntry++ < max_tries))
            {
                int idim = 2 - (ntry % 3);  // alternate between (2, 1, 0)
                if (this->tile_size[idim] == 1) { continue; }
                this->tile_size[idim] /= 2;
                n_logical = amrex::numTilesInBox(ba[gid], true, this->tile_size);
            }
        } else
        {
            n_logical = amrex::numTilesInBox(ba[gid], true, this->tile_size);
        }

        if (n_logical < nthreads)
        {
            std::string warning_message;
            if (!user_set) {
                warning_message.append(
                    "Could not find a good tile size for the requested " +
                    std::to_string(nthreads) + " OpenMP threads. "
                );
            }
            std::stringstream sstr_tile_size;
            sstr_tile_size << this->tile_size;
            warning_message.append(
                "The number of available tiles is " +
                std::to_string(n_logical) + " with each a size of " + sstr_tile_size.str() +
                ". Lowering the number of threads to match the available tiles now. This may result "
                "in poorer performance than expected. "
            );
            warning_message.append("You may want to try ");
            if (user_set) {
                warning_message.append("decreasing the particles.tile_size or ");
            }
            warning_message.append("increasing the blocking factor and max grid size.");
            ablastr::warn_manager::WMRecordWarning(
                "ImpactXParticleContainer::prepare",
                warning_message,
                ablastr::warn_manager::WarnPriority::medium
            );

#if defined(AMREX_USE_OMP)
            omp_set_num_threads(n_logical);
#endif
        }

        this->reserveData();
        this->resizeData();
    }

    template <class T_PT>
    void
    ImpactXParticleContainerT<T_PT>::clear (bool keep_mass, bool keep_charge)
    {
        this->clearParticles();
        this->reset_beam_moments_history();
        m_refpart.reset(keep_mass, keep_charge);
    }

    template <class T_PT>
    void
    ImpactXParticleContainerT<T_PT>::AddNParticles (
        amrex::Gpu::DeviceVector<amrex::ParticleReal> const & x,
        amrex::Gpu::DeviceVector<amrex::ParticleReal> const & y,
        amrex::Gpu::DeviceVector<amrex::ParticleReal> const & t,
        amrex::Gpu::DeviceVector<amrex::ParticleReal> const & px,
        amrex::Gpu::DeviceVector<amrex::ParticleReal> const & py,
        amrex::Gpu::DeviceVector<amrex::ParticleReal> const & pt,
        amrex::ParticleReal qm,
        std::optional<amrex::ParticleReal> bunch_charge,
        std::optional<amrex::Gpu::DeviceVector<amrex::ParticleReal>> w,
        std::optional<amrex::Gpu::DeviceVector<amrex::ParticleReal>> sx,
        std::optional<amrex::Gpu::DeviceVector<amrex::ParticleReal>> sy,
        std::optional<amrex::Gpu::DeviceVector<amrex::ParticleReal>> sz
    )
    {
        BL_PROFILE("ImpactX::AddNParticles");

        using namespace amrex::literals;  // for _rt and _prt

        // number of particles to add
        std::size_t const np_s = x.size();

        // input validation
        bool const has_w = w.has_value();
        if (!(bunch_charge.has_value() ^ has_w))
        {
            throw std::runtime_error("AddNParticles: Exactly one of bunch_charge or w must be provided!");
        }

        bool const has_spin = sx.has_value();

        AMREX_ALWAYS_ASSERT(np_s == y.size());
        AMREX_ALWAYS_ASSERT(np_s == t.size());
        AMREX_ALWAYS_ASSERT(np_s == px.size());
        AMREX_ALWAYS_ASSERT(np_s == py.size());
        AMREX_ALWAYS_ASSERT(np_s == pt.size());
        if (has_spin) {
            AMREX_ALWAYS_ASSERT(sy.has_value());
            AMREX_ALWAYS_ASSERT(sz.has_value());
            AMREX_ALWAYS_ASSERT(np_s == sy.value().size());
            AMREX_ALWAYS_ASSERT(np_s == sz.value().size());
        }
        if (has_w) { AMREX_ALWAYS_ASSERT(np_s == w->size()); }

        // we add particles to lev 0, grid 0
        int lid = 0, gid = 0;
        {
            const auto& pmap = this->ParticleDistributionMap(lid).ProcessorMap();
            auto it = std::find(pmap.begin(), pmap.end(), amrex::ParallelDescriptor::MyProc());
            if (it == std::end(pmap)) {
                throw std::runtime_error("AddNParticles: Attempting to add particles to box that does not exist.");
            } else {
                gid = *it;
            }
        }

        int nthreads = 1;
#if defined(AMREX_USE_OMP)
        nthreads = omp_get_max_threads();
#endif

        // split up particles over nthreads tiles
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(amrex::numTilesInBox(this->ParticleBoxArray(lid)[gid], true, this->tile_size) >= nthreads, "Not enough tiles for the number of OpenMP threads - please try reducing particles.tile_size or increasing the number of cells in the domain.");
        for (int ithr = 0; ithr < nthreads; ++ithr) {
            this->DefineAndReturnParticleTile(lid, gid, ithr);
        }

        amrex::Long const pid = base_t::ParticleType::NextID();
        amrex::Long const np = np_s;
        base_t::ParticleType::NextID(pid + np);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            pid + np < amrex::LongParticleIds::LastParticleID,
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
            ParticleChunk thread_chunk = split_equally(
                np,
                tid,
                nthreads
            );
            amrex::Long const my_offset = thread_chunk.offset;
            amrex::Long const num_to_add = thread_chunk.size;

            auto& particle_tile = this->ParticlesAt(lid, gid, tid);
            auto old_np = particle_tile.numParticles();
            auto new_np = old_np + num_to_add;
            particle_tile.resize(new_np);

            const amrex::Long cpuid = amrex::ParallelDescriptor::MyProc();

            auto & soa = particle_tile.GetStructOfArrays().GetRealData();
            // SoA storage real type is T_PT (single/double); deduce the pointer type
            auto * const AMREX_RESTRICT x_arr = soa[RealSoA::x].dataPtr();
            auto * const AMREX_RESTRICT y_arr = soa[RealSoA::y].dataPtr();
            auto * const AMREX_RESTRICT t_arr = soa[RealSoA::t].dataPtr();
            auto * const AMREX_RESTRICT px_arr = soa[RealSoA::px].dataPtr();
            auto * const AMREX_RESTRICT py_arr = soa[RealSoA::py].dataPtr();
            auto * const AMREX_RESTRICT pt_arr = soa[RealSoA::pt].dataPtr();
            auto * const AMREX_RESTRICT sx_arr = soa[RealSoA::sx].dataPtr();
            auto * const AMREX_RESTRICT sy_arr = soa[RealSoA::sy].dataPtr();
            auto * const AMREX_RESTRICT sz_arr = soa[RealSoA::sz].dataPtr();
            auto * const AMREX_RESTRICT qm_arr = soa[RealSoA::qm].dataPtr();
            auto * const AMREX_RESTRICT w_arr  = soa[RealSoA::w ].dataPtr();

            uint64_t * const AMREX_RESTRICT idcpu_arr = particle_tile.GetStructOfArrays().GetIdCPUData().dataPtr();

            amrex::ParticleReal const * const AMREX_RESTRICT x_ptr = x.data();
            amrex::ParticleReal const * const AMREX_RESTRICT y_ptr = y.data();
            amrex::ParticleReal const * const AMREX_RESTRICT t_ptr = t.data();
            amrex::ParticleReal const * const AMREX_RESTRICT px_ptr = px.data();
            amrex::ParticleReal const * const AMREX_RESTRICT py_ptr = py.data();
            amrex::ParticleReal const * const AMREX_RESTRICT pt_ptr = pt.data();
            amrex::ParticleReal const * const AMREX_RESTRICT sx_ptr = has_spin ? sx.value().data() : nullptr;
            amrex::ParticleReal const * const AMREX_RESTRICT sy_ptr = has_spin ? sy.value().data() : nullptr;
            amrex::ParticleReal const * const AMREX_RESTRICT sz_ptr = has_spin ? sz.value().data() : nullptr;
            amrex::ParticleReal const * const AMREX_RESTRICT w_ptr = has_w ? w->data() : nullptr;
            amrex::ParticleReal const bunch_charge_value = has_w ? 0_prt : bunch_charge.value();

            amrex::ParallelFor(num_to_add,
                [=] AMREX_GPU_DEVICE (amrex::Long i) noexcept
            {
                idcpu_arr[old_np+i] = amrex::SetParticleIDandCPU(pid + my_offset + i, cpuid);

                x_arr[old_np+i] = x_ptr[my_offset+i];
                y_arr[old_np+i] = y_ptr[my_offset+i];
                t_arr[old_np+i] = t_ptr[my_offset+i];

                px_arr[old_np+i] = px_ptr[my_offset+i];
                py_arr[old_np+i] = py_ptr[my_offset+i];
                pt_arr[old_np+i] = pt_ptr[my_offset+i];

                if (has_spin) {
                    sx_arr[old_np+i] = sx_ptr[my_offset+i];
                    sy_arr[old_np+i] = sy_ptr[my_offset+i];
                    sz_arr[old_np+i] = sz_ptr[my_offset+i];
                } else {
                    sx_arr[old_np+i] = 0_prt;
                    sy_arr[old_np+i] = 0_prt;
                    sz_arr[old_np+i] = 0_prt;
                }

                qm_arr[old_np+i] = qm;
                w_arr[old_np+i]  = has_w ? w_ptr[my_offset+i] : bunch_charge_value / ablastr::constant::SI::q_e/np;
            });
        }

        // safety first: in case passed attribute arrays were temporary, we
        // want to make sure the ParallelFor has ended here
        amrex::Gpu::streamSynchronize();
    }

    template <class T_PT>
    void
    ImpactXParticleContainerT<T_PT>::SetRefParticle (RefPart const & refpart)
    {
        m_refpart = refpart;
    }

    template <class T_PT>
    RefPart &
    ImpactXParticleContainerT<T_PT>::GetRefParticle ()
    {
        return m_refpart;
    }

    template <class T_PT>
    RefPart const &
    ImpactXParticleContainerT<T_PT>::GetRefParticle () const
    {
        return m_refpart;
    }

    template <class T_PT>
    void
    ImpactXParticleContainerT<T_PT>::SetRefParticleEdge ()
    {
        m_refpart.sedge = m_refpart.s;
    }

    template <class T_PT>
    std::tuple<
            amrex::ParticleReal, amrex::ParticleReal,
            amrex::ParticleReal, amrex::ParticleReal,
            amrex::ParticleReal, amrex::ParticleReal>
    ImpactXParticleContainerT<T_PT>::MinAndMaxPositions ()
    {
        BL_PROFILE("ImpactXParticleContainer::MinAndMaxPositions");
        return ablastr::particles::MinAndMaxPositions(*this);
    }

    template <class T_PT>
    std::tuple<
            amrex::ParticleReal, amrex::ParticleReal,
            amrex::ParticleReal, amrex::ParticleReal,
            amrex::ParticleReal, amrex::ParticleReal>
    ImpactXParticleContainerT<T_PT>::MeanAndStdPositions ()
    {
        BL_PROFILE("ImpactXParticleContainer::MeanAndStdPositions");
        return ablastr::particles::MeanAndStdPositions<
            ImpactXParticleContainerT<T_PT>, RealSoA::w
        >(*this);
    }

    template <class T_PT>
    CoordSystem
    ImpactXParticleContainerT<T_PT>::GetCoordSystem () const
    {
        return m_coordsystem;
    }

    template <class T_PT>
    void
    ImpactXParticleContainerT<T_PT>::SetCoordSystem (CoordSystem coord_system)
    {
        m_coordsystem = coord_system;
    }

    template <class T_PT>
    void
    ImpactXParticleContainerT<T_PT>::record_beam_moments ()
    {
        BL_PROFILE("ImpactXParticleContainer::record_beam_moments");

        auto rbc = diagnostics::reduced_beam_characteristics(*this);
        amrex::ParticleReal const s = this->GetRefParticle().s;
        rbc["s"] = s;

        m_beam_moments.push_back(rbc);
    }

    template <class T_PT>
    std::unordered_map<std::string, amrex::ParticleReal>
    ImpactXParticleContainerT<T_PT>::beam_moments ()
    {
        BL_PROFILE("ImpactXParticleContainer::beam_moments");

        auto rbc = diagnostics::reduced_beam_characteristics(*this);
        amrex::ParticleReal const s = this->GetRefParticle().s;
        rbc["s"] = s;

        return rbc;
    }

    // explicit template instantiations for the precisions ImpactX provides
    template class ParIterSoAT<IMPACTX_PARTICLE_REAL>;
    template class ParConstIterSoAT<IMPACTX_PARTICLE_REAL>;
    template class ImpactXParticleContainerT<IMPACTX_PARTICLE_REAL>;
} // namespace impactx
