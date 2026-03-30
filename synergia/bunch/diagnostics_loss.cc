#include "synergia/bunch/diagnostics_loss.h"
#include "synergia/bunch/bunch.h"
#include "synergia/utils/simple_timer.h"
#include <stdexcept>

void
Diagnostics_loss::do_update(Bunch const& bunch)
{
    scoped_simple_timer timer("diag_loss_update");

    auto ref = bunch.get_reference_particle();
    repetition = ref.get_repetition();
    s_ref_particle = ref.get_s();
    sn_ref_particle = ref.get_s_n();
    bucket_index = bunch.get_bucket_index();

    coords = bunch.get_particles_last_discarded();
#ifdef SYNERGIA_HAVE_OPENPMD
    discards_local_num = coords.extent(0);

    if (MPI_Allreduce(&discards_local_num,
                      &discards_total_num,
                      1,
                      MPI_SIZE_T,
                      MPI_SUM,
                      bunch.get_comm()) != MPI_SUCCESS) {
        std::runtime_error("Error in MPI_Allreduce in diagnostics-loss!");
    }

    if (MPI_Scan(&discards_local_num,
                 &discards_offset_num,
                 1,
                 MPI_SIZE_T,
                 MPI_SUM,
                 bunch.get_comm()) != MPI_SUCCESS) {
        std::runtime_error("Error in MPI_Scan in diagnostics-loss!");
    }

#endif
}

void
Diagnostics_loss::do_first_write(io_device& file)
{
#ifdef SYNERGIA_HAVE_OPENPMD
    file.setAttribute("bucket_index", bucket_index);
#else
    file.append("bucket_index", bucket_index);
#endif
    return;
}

void
Diagnostics_loss::do_write(io_device& file, size_t iteration)
{
    scoped_simple_timer timer("diag_loss_write");
#ifdef SYNERGIA_HAVE_OPENPMD
    auto i = file.iterations[iteration];
    i.setAttribute("repetition", repetition);
    i.setAttribute("s", s_ref_particle);
    i.setAttribute("s_n", sn_ref_particle);

    openPMD::ParticleSpecies& protons =
        file.iterations[iteration].particles["bunch_discards"];

    openPMD::Datatype datatype = openPMD::determineDatatype<double>();
    openPMD::Extent global_extent = {discards_total_num};
    openPMD::Dataset dataset = openPMD::Dataset(datatype, global_extent);

    auto parts_x = Kokkos::subview(coords, Kokkos::ALL, Bunch::x);
    auto parts_y = Kokkos::subview(coords, Kokkos::ALL, Bunch::y);
    auto parts_z = Kokkos::subview(coords, Kokkos::ALL, Bunch::cdt);

    auto moments_x = Kokkos::subview(coords, Kokkos::ALL, Bunch::xp);
    auto moments_y = Kokkos::subview(coords, Kokkos::ALL, Bunch::yp);
    auto moments_dpop = Kokkos::subview(coords, Kokkos::ALL, Bunch::dpop);

    protons["position"]["x"].resetDataset(dataset);
    protons["position"]["y"].resetDataset(dataset);
    protons["position"]["z"].resetDataset(dataset);

    protons["moments"]["x"].resetDataset(dataset);
    protons["moments"]["y"].resetDataset(dataset);
    protons["moments"]["z"].resetDataset(dataset);

    file.flush();

    openPMD::Offset chunk_offset = {discards_offset_num - discards_local_num};
    openPMD::Extent chunk_extent = {discards_local_num};

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

    file.flush();

#else
    file.append("repetition", repetition);
    file.append("s", s_ref_particle);
    file.append("s_n", sn_ref_particle);

    file.write_collective("coordinates_" + std::to_string(s_ref_particle),
                          coords);
#endif

    return;
}
