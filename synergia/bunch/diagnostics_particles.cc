#include "synergia/foundation/physical_constants.h"

#include "synergia/bunch/bunch.h"
#include "synergia/bunch/bunch_particles.h"
#include "synergia/bunch/diagnostics.h"
#include "synergia/bunch/diagnostics_particles.h"

Diagnostics_particles::Diagnostics_particles(std::string const& filename,
                                             int num_part,
                                             int offset,
                                             int num_spec_part,
                                             int spec_offset)
#ifdef SYNERGIA_HAVE_OPENPMD
    // use a single file and iterations with openPMD
    : Diagnostics("diagnostics_particles", filename, true)
#else
    // use a new file for each write with old HDF5 backend
    : Diagnostics("diagnostics_particles", filename, false)
#endif
    , bunch_ref(std::nullopt)
    , num_part(num_part)
    , offset(offset)
    , num_spec_part(num_spec_part)
    , spec_offset(spec_offset)
{}

void
Diagnostics_particles::do_first_write(io_device& file)
{
#ifdef SYNERGIA_HAVE_OPENPMD
    assert(bunch_ref.has_value());
    auto const& ref_part = bunch_ref.value().get().get_reference_particle();
    file.setAttribute("charge", ref_part.get_charge());
    file.setAttribute("mass", ref_part.get_mass());
    file.setAttribute("s", ref_part.get_s());
    file.setAttribute("s_n", ref_part.get_s_n());
    file.setAttribute("repetition", ref_part.get_repetition());
    file.setAttribute("abs_offset", ref_part.get_bunch_abs_offset());
    file.setAttribute("abs_time", ref_part.get_bunch_abs_time());
    file.setAttribute("pz", ref_part.get_momentum());
    file.flush();

    auto bunch_parts =
        bunch_ref.value().get().get_bunch_particles(ParticleGroup::regular);
    auto bunch_spec_parts =
        bunch_ref.value().get().get_bunch_particles(ParticleGroup::spectator);

    std::tie(local_num, local_offset) =
        bunch_ref.value().get().get_local_particle_count_in_range(
            ParticleGroup::regular, num_part, offset);

    if (local_num < 0 || local_offset < 0 ||
        local_num + local_offset > bunch_parts.num_active()) {
        throw std::runtime_error("invalid num_part or offset for "
                                 "diagnostics_particles!");
    }

    std::tie(spec_local_num, spec_local_offset) =
        bunch_ref.value().get().get_local_particle_count_in_range(
            ParticleGroup::spectator, num_spec_part, spec_offset);
    if (spec_local_num < 0 || spec_local_offset < 0 ||
        spec_local_num + spec_local_offset > bunch_spec_parts.num_active()) {
        throw std::runtime_error("invalid num_part or offset for "
                                 "diagnostics_particles!");
    }

    if (MPI_Scan(&local_num,
                 &file_offset,
                 1,
                 MPI_SIZE_T,
                 MPI_SUM,
                 bunch_ref.value().get().get_comm()) != MPI_SUCCESS) {
        std::runtime_error("Error in MPI_Scan in diagnostics-particles!");
    }

    if (MPI_Scan(&spec_local_num,
                 &spec_file_offset,
                 1,
                 MPI_SIZE_T,
                 MPI_SUM,
                 bunch_ref.value().get().get_comm()) != MPI_SUCCESS) {
        std::runtime_error("Error in MPI_Scan in diagnostics-particles!");
    }

    parts_subset = bunch_particles_t<double>::host_parts_t(
        "diag_particles_parts", local_num);
    masks_subset = bunch_particles_t<double>::host_masks_t(
        "diag_particles_masks", local_num);

    spec_parts_subset = bunch_particles_t<double>::host_parts_t(
        "diag_particles_spec_parts", spec_local_num);
    spec_masks_subset = bunch_particles_t<double>::host_masks_t(
        "diag_particles_spec_masks", spec_local_num);

#else
    // nothing to do if using the old HDF5 backend!
#endif
    return;
}

void
Diagnostics_particles::do_write(io_device& file, const size_t iteration)
{
    assert(bunch_ref.has_value());
    auto const& ref_part = bunch_ref.value().get().get_reference_particle();

#ifdef SYNERGIA_HAVE_OPENPMD

    // writing particles
    if (local_num > 0) {

        if (num_part == -1) {
            num_part =
                bunch_ref.value().get().get_total_num(ParticleGroup::regular);
        }

        openPMD::ParticleSpecies& protons =
            file.iterations[iteration].particles["bunch_particles"];

        // write mass in SI units!
        protons.setAttribute(
            "mass", bunch_ref.value().get().get_mass() / pconstants::kg_to_GeV);
        protons.setAttribute(
            "beta_ref",
            (bunch_ref.value().get().get_reference_particle()).get_beta());
        protons.setAttribute(
            "gamma_ref",
            (bunch_ref.value().get().get_reference_particle()).get_gamma());

        // OpenPMD2.0: EXT_BeamPhysics
        protons.setAttribute("numParticles",
                             bunch_ref.value().get().get_total_num());
        protons.setAttribute("speciesType", std::string("proton"));

        openPMD::Datatype datatype = openPMD::determineDatatype<double>();
        openPMD::Extent global_extent = {static_cast<size_t>(num_part)};
        openPMD::Dataset dataset = openPMD::Dataset(datatype, global_extent);

        openPMD::Datatype masks_datatype =
            openPMD::determineDatatype<uint8_t>();
        openPMD::Dataset masks_dataset =
            openPMD::Dataset(masks_datatype, global_extent);

        bunch_ref.value()
            .get()
            .get_bunch_particles(ParticleGroup::regular)
            .get_particles_in_range(
                parts_subset, masks_subset, local_num, local_offset);

        auto parts_x = Kokkos::subview(parts_subset, Kokkos::ALL, Bunch::x);
        auto parts_y = Kokkos::subview(parts_subset, Kokkos::ALL, Bunch::y);
        auto parts_z = Kokkos::subview(parts_subset, Kokkos::ALL, Bunch::cdt);

        auto moments_x = Kokkos::subview(parts_subset, Kokkos::ALL, Bunch::xp);
        auto moments_y = Kokkos::subview(parts_subset, Kokkos::ALL, Bunch::yp);
        auto moments_dpop =
            Kokkos::subview(parts_subset, Kokkos::ALL, Bunch::dpop);

        auto parts_idx = Kokkos::subview(parts_subset, Kokkos::ALL, Bunch::id);

        protons["position"]["x"].resetDataset(dataset);
        protons["position"]["y"].resetDataset(dataset);
        protons["position"]["z"].resetDataset(dataset);

        protons["positionOffset"]["x"].resetDataset(dataset);
        protons["positionOffset"]["y"].resetDataset(dataset);
        protons["positionOffset"]["z"].resetDataset(dataset);

        protons["moments"]["x"].resetDataset(dataset);
        protons["moments"]["y"].resetDataset(dataset);
        protons["moments"]["z"].resetDataset(dataset);

        protons["id"][openPMD::RecordComponent::SCALAR].resetDataset(dataset);
        protons["masks"][openPMD::RecordComponent::SCALAR].resetDataset(
            masks_dataset);

        openPMD::Offset chunk_offset = {file_offset - local_num};
        openPMD::Extent chunk_extent = {local_num};

        protons["position"]["x"].storeChunkRaw(
            parts_x.data(), chunk_offset, chunk_extent);
        protons["position"]["y"].storeChunkRaw(
            parts_y.data(), chunk_offset, chunk_extent);
        protons["position"]["z"].storeChunkRaw(
            parts_z.data(), chunk_offset, chunk_extent);

        protons["moments"]["x"].storeChunkRaw(
            moments_x.data(), chunk_offset, chunk_extent);
        protons["moments"]["y"].storeChunkRaw(
            moments_y.data(), chunk_offset, chunk_extent);
        protons["moments"]["z"].storeChunkRaw(
            moments_dpop.data(), chunk_offset, chunk_extent);

        protons["id"][openPMD::RecordComponent::SCALAR].storeChunkRaw(
            parts_idx.data(), chunk_offset, chunk_extent);

        protons["masks"][openPMD::RecordComponent::SCALAR].storeChunkRaw(
            masks_subset.data(), chunk_offset, chunk_extent);

        protons["positionOffset"]["x"].makeConstant(0);
        protons["positionOffset"]["y"].makeConstant(0);
        protons["positionOffset"]["z"].makeConstant(0);

        file.flush();
    }

    // writing spectator particles
    if (spec_local_num > 0) {

        if (num_spec_part == -1) {
            num_spec_part =
                bunch_ref.value().get().get_total_num(ParticleGroup::spectator);
        }

        openPMD::ParticleSpecies& protons =
            file.iterations[iteration].particles["bunch_spectator_particles"];

        openPMD::Datatype datatype = openPMD::determineDatatype<double>();
        openPMD::Extent global_extent = {static_cast<size_t>(num_spec_part)};
        openPMD::Dataset dataset = openPMD::Dataset(datatype, global_extent);

        openPMD::Datatype masks_datatype =
            openPMD::determineDatatype<uint8_t>();
        openPMD::Dataset masks_dataset =
            openPMD::Dataset(masks_datatype, global_extent);

        bunch_ref.value()
            .get()
            .get_bunch_particles(ParticleGroup::spectator)
            .get_particles_in_range(spec_parts_subset,
                                    spec_masks_subset,
                                    spec_local_num,
                                    spec_local_offset);
        auto parts_x =
            Kokkos::subview(spec_parts_subset, Kokkos::ALL, Bunch::x);
        auto parts_y =
            Kokkos::subview(spec_parts_subset, Kokkos::ALL, Bunch::y);
        auto parts_z =
            Kokkos::subview(spec_parts_subset, Kokkos::ALL, Bunch::cdt);

        auto moments_x =
            Kokkos::subview(spec_parts_subset, Kokkos::ALL, Bunch::xp);
        auto moments_y =
            Kokkos::subview(spec_parts_subset, Kokkos::ALL, Bunch::yp);
        auto moments_dpop =
            Kokkos::subview(spec_parts_subset, Kokkos::ALL, Bunch::dpop);

        auto parts_idx = Kokkos::subview(parts_subset, Kokkos::ALL, Bunch::id);

        protons["position"]["x"].resetDataset(dataset);
        protons["position"]["y"].resetDataset(dataset);
        protons["position"]["z"].resetDataset(dataset);

        protons["positionOffset"]["x"].resetDataset(dataset);
        protons["positionOffset"]["y"].resetDataset(dataset);
        protons["positionOffset"]["z"].resetDataset(dataset);

        protons["moments"]["x"].resetDataset(dataset);
        protons["moments"]["y"].resetDataset(dataset);
        protons["moments"]["z"].resetDataset(dataset);

        protons["id"][openPMD::RecordComponent::SCALAR].resetDataset(dataset);
        protons["masks"][openPMD::RecordComponent::SCALAR].resetDataset(
            masks_dataset);

        openPMD::Offset chunk_offset = {spec_file_offset - spec_local_num};
        openPMD::Extent chunk_extent = {spec_local_num};

        protons["position"]["x"].storeChunkRaw(
            parts_x.data(), chunk_offset, chunk_extent);
        protons["position"]["y"].storeChunkRaw(
            parts_y.data(), chunk_offset, chunk_extent);
        protons["position"]["z"].storeChunkRaw(
            parts_z.data(), chunk_offset, chunk_extent);

        protons["moments"]["x"].storeChunkRaw(
            moments_x.data(), chunk_offset, chunk_extent);
        protons["moments"]["y"].storeChunkRaw(
            moments_y.data(), chunk_offset, chunk_extent);
        protons["moments"]["z"].storeChunkRaw(
            moments_dpop.data(), chunk_offset, chunk_extent);

        protons["id"][openPMD::RecordComponent::SCALAR].storeChunkRaw(
            parts_idx.data(), chunk_offset, chunk_extent);

        protons["masks"][openPMD::RecordComponent::SCALAR].storeChunkRaw(
            masks_subset.data(), chunk_offset, chunk_extent);

        protons["positionOffset"]["x"].makeConstant(0);
        protons["positionOffset"]["y"].makeConstant(0);
        protons["positionOffset"]["z"].makeConstant(0);

        file.flush();
    }

#else
    file.write("charge", ref_part.get_charge());
    file.write("mass", ref_part.get_four_momentum().get_mass());

    file.write("s", ref_part.get_s());
    file.write("s_n", ref_part.get_s_n());
    file.write("repetition", ref_part.get_repetition());
    file.write("abs_offset", ref_part.get_bunch_abs_offset());
    file.write("abs_time", ref_part.get_bunch_abs_time());
    file.write("pz", ref_part.get_momentum());

    bunch_ref.value().get().write_file(
        file, num_part, offset, num_spec_part, spec_offset);
#endif
    // reset bunch_ref
    bunch_ref = std::nullopt;
    return;
}
