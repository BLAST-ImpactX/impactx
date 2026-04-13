#include "diagnostics_bulk_track.h"
#include "synergia/bunch/bunch.h"
#include "synergia/utils/parallel_utils.h"
#include "synergia/utils/simple_timer.h"
#include <iostream>
#include <string>

Diagnostics_bulk_track::Diagnostics_bulk_track(std::string const& filename,
                                               int num_tracks,
                                               int offset,
                                               ParticleGroup pg)
    : Diagnostics("diagnostis_bulk_track", filename, true)
    , total_num_tracks(num_tracks)
    , local_num_tracks(0)
    , offset(offset)
    , local_offset(0)
    , setup(false)
    , track_coords("local_coords", 0, 0)
    , pg(pg)
{}

void
Diagnostics_bulk_track::do_update(Bunch const& bunch)
{
    scoped_simple_timer("diag_bulk_track_update");

    auto const& ref = bunch.get_reference_particle();

    if (!setup) {
        ref_charge = ref.get_charge();
        ref_mass = ref.get_four_momentum().get_mass();
        ref_pz = ref.get_four_momentum().get_momentum();


        auto const& comm = bunch.get_comm();

        local_num_tracks = decompose_1d_local(comm, total_num_tracks);
        local_offset = decompose_1d_local(comm, offset);

        if (local_num_tracks + local_offset > bunch.size(pg))
            local_num_tracks = bunch.size(pg) - local_offset;

        setup = true;
    }

    pz = ref.get_momentum();
    s = ref.get_s();
    s_n = ref.get_s_n();
    repetition = ref.get_repetition();
    bunch_abs_offset = ref.get_bunch_abs_offset();
    bunch_abs_time = ref.get_bunch_abs_time();


    track_coords =
        bunch.get_particles_in_range_row(local_offset, local_num_tracks, pg);

#ifdef SYNERGIA_HAVE_OPENPMD
    track_parts_local_num = track_coords.extent(0);

    if (MPI_Allreduce(&track_parts_local_num,
                      &track_parts_total_num,
                      1,
                      MPI_SIZE_T,
                      MPI_SUM,
                      bunch.get_comm()) != MPI_SUCCESS) {
        std::runtime_error("Error in MPI_Allreduce in diagnostics-bulk-track!");
    }

    if (MPI_Scan(&track_parts_local_num,
                 &track_parts_offset_num,
                 1,
                 MPI_SIZE_T,
                 MPI_SUM,
                 bunch.get_comm()) != MPI_SUCCESS) {
        std::runtime_error("Error in MPI_Scan in diagnostics-bulk-track!");
    }

#endif
}

void
Diagnostics_bulk_track::do_first_write(io_device& file)
{
#ifdef SYNERGIA_HAVE_OPENPMD
    file.setAttribute("charge", ref_charge);
    file.setAttribute("mass", ref_mass);
    file.setAttribute("pz", ref_pz);
    file.setAttribute("bunch_abs_offset", bunch_abs_offset);
#else
    file.write("charge", ref_charge);
    file.write("mass", ref_mass);
    file.write("pz", ref_pz);
    file.write("bunch_abs_offset", bunch_abs_offset);
#endif
}

void
Diagnostics_bulk_track::do_write(io_device& file, const size_t iteration)
{
    scoped_simple_timer("diag_bulk_track_write");

#ifdef SYNERGIA_HAVE_OPENPMD
    auto i = file.iterations[iteration];
    i.setAttribute("track_pz", pz);
    i.setAttribute("track_s", s);
    i.setAttribute("track_s_n", s_n);
    i.setAttribute("track_repitition", repetition);
    i.setAttribute("bunch_abs_time", bunch_abs_time);

    openPMD::ParticleSpecies& protons =
        file.iterations[iteration].particles["track_coords"];

    openPMD::Datatype datatype = openPMD::determineDatatype<double>();
    openPMD::Extent global_extent = {track_parts_total_num};
    openPMD::Dataset dataset = openPMD::Dataset(datatype, global_extent);

    auto parts_x = Kokkos::subview(track_coords, Kokkos::ALL, Bunch::x);
    auto parts_y = Kokkos::subview(track_coords, Kokkos::ALL, Bunch::y);
    auto parts_z = Kokkos::subview(track_coords, Kokkos::ALL, Bunch::cdt);

    auto moments_x = Kokkos::subview(track_coords, Kokkos::ALL, Bunch::xp);
    auto moments_y = Kokkos::subview(track_coords, Kokkos::ALL, Bunch::yp);
    auto moments_dpop = Kokkos::subview(track_coords, Kokkos::ALL, Bunch::dpop);

    protons["position"]["x"].resetDataset(dataset);
    protons["position"]["y"].resetDataset(dataset);
    protons["position"]["z"].resetDataset(dataset);

    protons["moments"]["x"].resetDataset(dataset);
    protons["moments"]["y"].resetDataset(dataset);
    protons["moments"]["z"].resetDataset(dataset);

    file.flush();

    openPMD::Offset chunk_offset = {track_parts_offset_num -
                                    track_parts_local_num};
    openPMD::Extent chunk_extent = {track_parts_local_num};

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
    // write serial
    file.append_single("track_pz", pz);
    file.append_single("track_s", s);
    file.append_single("track_s_n", s_n);
    file.append_single("track_repetition", repetition);
    file.append_single("bunch_abs_time", bunch_abs_time);

    // write collective from all ranks
    file.append_collective("track_coords", track_coords);
#endif
}
