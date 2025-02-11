/* Copyright 2022-2025 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "elements/Source.H"

#include <AMReX_REAL.H>

#ifdef ImpactX_USE_OPENPMD
#   include "elements/diagnostics/openPMD.H"
#   include <openPMD/openPMD.hpp>
namespace io = openPMD;
#else
#  include <stdexcept>
#endif


namespace impactx::elements
{
    void
    Source::operator() (
        ImpactXParticleContainer & pc,
        int,
        int
    )
    {
#ifdef ImpactX_USE_OPENPMD
        auto series = io::Series(m_series_name, io::Access::READ_ONLY
#   if openPMD_HAVE_MPI==1
            , amrex::ParallelDescriptor::Communicator()
#   endif
        );

        // read the last iteration (highest openPMD iteration)
        // TODO: later we can make this an option
        int const read_iteration = std::prev(series.iterations.end())->first;
        io::Iteration iteration = series.iterations[read_iteration];

        // TODO: later we can make the particle species name an option
        std::string const species_name = "beam";
        io::ParticleSpecies beam = iteration.particles[species_name];
        // TODO: later we can make the bunch charge an option (i.e., allow rescaling a distribution)
        amrex::ParticleReal bunch_charge = beam.getAttribute("charge_C").get<amrex::ParticleReal>();

        auto const scalar = openPMD::RecordComponent::SCALAR;
        auto const getComponentRecord = [&beam](std::string comp_name) {
            return diagnostics::detail::get_component_record(beam, std::move(comp_name));
        };

        int const npart = beam["id"][scalar].getExtent()[0];  // how many particles to read total

        // TODO: read reference particle (optional?)

        // read the particles

        // Logic: We initialize 1/Nth of particles, independent of their
        // position, per MPI rank. We then measure the distribution's spatial
        // extent, create a grid, resize it to fit the beam, and then
        // redistribute particles so that they reside on the correct MPI rank.
        int const myproc = amrex::ParallelDescriptor::MyProc();
        int const nprocs = amrex::ParallelDescriptor::NProcs();
        int const navg = npart / nprocs;  // note: integer division
        int const nleft = npart - navg * nprocs;
        std::uint64_t const npart_this_proc = (myproc < nleft) ? navg+1 : navg;  // add 1 to each proc until distributed
        std::uint64_t npart_before_this_proc = (myproc < nleft) ? (navg+1) * myproc : navg * myproc + nleft;
        auto const rel_part_this_proc =
            amrex::ParticleReal(npart_this_proc) / amrex::ParticleReal(npart);

        // alloc data for particle attributes
        std::map<std::string, amrex::Gpu::PinnedVector<amrex::ParticleReal>> pinned_SoA;
        std::vector<std::string> real_soa_names = pc.GetRealSoANames();
        for (auto real_idx = 0; real_idx < pc.NumRealComps(); real_idx++) {
            pinned_SoA[real_soa_names.at(real_idx)].resize(npart_this_proc);
        }

        // read from file
        // idcpu: TODO
        //amrex::Gpu::PinnedVector<std::uint64_t> pinned_idcpu(npart_this_proc);
        //beam["id"][scalar].loadChunkRaw(pinned_idcpu.data(), {npart_before_this_proc}, {npart_this_proc});
        // SoA: Real
        {
            for (auto real_idx = 0; real_idx < pc.NumRealComps(); real_idx++) {
                auto const component_name = real_soa_names.at(real_idx);
                getComponentRecord(component_name).loadChunkRaw(
                    pinned_SoA[component_name].data(),
                    {npart_before_this_proc},
                    {npart_this_proc}
                );
            }
        }
        // SoA: Int
        std::vector<std::string> int_soa_names = pc.GetIntSoANames();
        static_assert(IntSoA::nattribs == 0); // not yet used
        if (!int_soa_names.empty())
            throw std::runtime_error("BeamMonitor: int_soa_names output not yet implemented!");

        series.flush();
        series.close();

        // copy to device
        std::map<std::string, amrex::Gpu::DeviceVector<amrex::ParticleReal>> d_SoA;
        for (auto real_idx = 0; real_idx < pc.NumRealComps(); real_idx++) {
            auto const component_name = real_soa_names.at(real_idx);
            d_SoA[component_name].resize(npart_this_proc);
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, pinned_SoA[component_name].begin(), pinned_SoA[component_name].end(), d_SoA[component_name].begin());
        }

        // finalize distributions and deallocate temporary device global memory
        amrex::Gpu::streamSynchronize();

        // TODO: at this point, we ignore the "id", "qm" and "weighting" in the file. Could be improved
        pc.AddNParticles(d_SoA["position_x"], d_SoA["position_y"], d_SoA["position_t"],
                         d_SoA["momentum_x"], d_SoA["momentum_y"], d_SoA["momentum_t"],
                         pc.GetRefParticle().qm_ratio_SI(),
                         bunch_charge * rel_part_this_proc);

#else  // ImpactX_USE_OPENPMD
        amrex::ignore_unused(pc);
        throw std::runtime_error("BeamMonitor: openPMD not compiled");
#endif  // ImpactX_USE_OPENPMD
    }

} // namespace impactx
