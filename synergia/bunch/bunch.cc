
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/core_diagnostics.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/utils/parallel_utils.h"
#include "synergia/utils/synergia_config.h"

#include <iostream>
#include <synergia_version.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

template <>
Bunch::bunch_t(Reference_particle const& reference_particle,
               int total_num,
               double real_num,
               Commxx bunch_comm,
               int total_spectator_num,
               int bunch_index,
               int bucket_index,
               int array_index,
               int train_index)
    : comm(std::make_shared<Commxx>(bunch_comm))
    , boundary(LB::open)
    , boundary_param(0.0)
    , ref_part(reference_particle)
    , design_ref_part(reference_particle)
    , particle_charge(reference_particle.get_charge())
    , real_num(real_num)
    , parts{BunchParticles(PG::regular, total_num, -1, *comm),
            BunchParticles(PG::spectator, total_spectator_num, -1, *comm)}
    , bunch_index(bunch_index)
    , bucket_index(bucket_index)
    , array_index(array_index)
    , train_index(train_index)
{}

template <>
Bunch::bunch_t()
    : comm(new Commxx())
    , boundary(LB::open)
    , boundary_param(0.0)
    , ref_part()
    , design_ref_part()
    , particle_charge(ref_part.get_charge())
    , real_num(1.0)
    , parts{BunchParticles(PG::regular, 0, 0, *comm),
            BunchParticles(PG::spectator, 0, 0, *comm)}
    , bunch_index(0)
    , bucket_index(0)
    , array_index(0)
    , train_index(0)
{}

template <>
void
Bunch::inject(Bunch const& o)
{
    const double weight_tolerance = 1.0e-14;
    const double particle_tolerance = 1.0e-14;

    auto const& ref = ref_part;
    auto const& oref = o.ref_part;

    // The charge and mass of the bunch particles must match
    if (particle_charge != o.particle_charge) {
        throw std::runtime_error(
            "Bunch.inject: bunch particle charges do not match.");
    }

    if (std::abs(ref.get_mass() / oref.get_mass() - 1.0) > particle_tolerance) {
        throw std::runtime_error(
            "Bunch:inject: bunch particle masses do not match.");
    }

    // total num for regualr particles
    int total_num = parts[0].num_total();

    // can only check particle weight if total_num is nonzero
    if (total_num != 0) {
        double wgt1 = real_num / total_num;
        double wgt2 = o.get_real_num() / o.get_total_num();

        if (std::abs(wgt1 / wgt2 - 1.0) > weight_tolerance) {
            throw std::runtime_error("Bunch.inject: macroparticle weight of "
                                     "injected bunch does not match.");
        }
    }

    // prepare
    double target_momentum = ref.get_momentum();
    double injected_momentum = oref.get_momentum();
    double pdiff = injected_momentum / target_momentum;

    karray1d_dev ref_st_diff("", 6);
    karray1d_dev tgt_st("", 6);
    karray1d_dev inj_st("", 6);

    auto h_ref_st_diff = Kokkos::create_mirror_view(ref_st_diff);
    auto h_tgt_st = Kokkos::create_mirror_view(tgt_st);
    auto h_inj_st = Kokkos::create_mirror_view(inj_st);

    for (int i = 0; i < 6; ++i) {
        h_tgt_st(i) = ref.get_state()[i];
        h_inj_st(i) = oref.get_state()[i];
        h_ref_st_diff(i) = h_inj_st(i) - h_tgt_st(i);
    }

    Kokkos::deep_copy(ref_st_diff, h_ref_st_diff);
    Kokkos::deep_copy(tgt_st, h_tgt_st);
    Kokkos::deep_copy(inj_st, h_inj_st);

    // regular
    get_bunch_particles(PG::regular)
        .inject(o.get_bunch_particles(PG::regular),
                ref_st_diff,
                tgt_st,
                inj_st,
                pdiff);

    // spectator
    get_bunch_particles(PG::spectator)
        .inject(o.get_bunch_particles(PG::spectator),
                ref_st_diff,
                tgt_st,
                inj_st,
                pdiff);

    // update total number, for both real and spectator particles
    int old_total = update_total_num();

    // target bunch is empty.  Set the weights from the injected bunch
    if (old_total == 0) real_num = o.get_real_num();
}

template <>
void
Bunch::print_statistics(Logger& logger) const
{
    using LV = LoggerV;

    logger(LV::DEBUG) << "Bunch statistics: "
                      << "num_valid = " << get_local_num()
                      << ", size = " << size() << ", capacity = " << capacity()
                      << ", total_num = " << get_total_num()
                      << "\nMean and std: ";

    // print particles after propagate
    auto mean = Core_diagnostics::calculate_mean(*this);
    auto std = Core_diagnostics::calculate_std(*this, mean);

    logger(LV::DEBUG) << std::resetiosflags(std::ios::fixed)
                      << std::setprecision(16)
                      << std::setiosflags(std::ios::showpos |
                                          std::ios::scientific)
                      << "\n"
        //<< "\nmean\tstd\n"
        ;

    for (int i = 0; i < 6; ++i)
        logger(LV::DEBUG) << mean[i] << ", " << std[i] << "\n";

    logger(LV::DEBUG) << std::resetiosflags(std::ios::showpos |
                                            std::ios::scientific)
                      << "\n";
}

template <>
void
Bunch::adjust_bunch_particles_reference_energy(double newE)
{
  FF_adjust_bunch_ref_coords::apply(newE, *this);
}

template <>
void
Bunch::offset_bunch_particles_cdt(double offset)
{
  FF_offset_cdt::apply(offset, *this);
}


#if defined SYNERGIA_HAVE_OPENPMD
template <>
void
Bunch::read_openpmd_file(std::string const& filename,
                         std::optional<std::size_t> idx)
{
    auto series = openPMD::Series(filename, openPMD::Access::READ_ONLY);

    // get the last iteration from the series

    std::size_t index_to_read;

    if (idx.has_value()) {
        index_to_read = idx.value();
    } else {
        auto iters = series.iterations;
        // C++ map.end is one-past the end!
        index_to_read = (--iters.cend())->first;
    }

    auto iteration = series.iterations[index_to_read];

    // determine which software created the OpenPMD file
    std::string software_name;
    try {
        software_name = series.software();
    }
    catch (openPMD::error::NoSuchAttribute) {
        std::cerr
            << "Could not determine which software generated the input! \n";
        return;
    }

    if (software_name == "synergia3") {
        // input was created by synergia3!

        openPMD::ParticleSpecies& protons =
            iteration.particles["bunch_particles"];
        openPMD::ParticleSpecies& masks =
            iteration.particles["bunch_particles_masks"];
        auto parts_ids = protons["id"][openPMD::RecordComponent::SCALAR];
        auto extent = parts_ids.getExtent();
        //  The extent should be a vector of one element!
        auto num_part = extent[0];
        assert((extent.size() == 1) &&
               "Extent of particle_ids extent should be 1-dimensional!");

        double mass = protons.getAttribute("mass").get<double>();
        // convert mass to GeV!
        mass = mass * pconstants::kg_to_GeV;
        double beta_ref = protons.getAttribute("beta_ref").get<double>();
        double gamma_ref = protons.getAttribute("gamma_ref").get<double>();

        Reference_particle& ref_part = this->get_reference_particle();
        Four_momentum fm = Four_momentum(mass);
        fm.set_beta(beta_ref);
        fm.set_gamma(gamma_ref);
        ref_part.set_four_momentum(fm);

        // reserver required space
        this->get_bunch_particles(PG::regular)
            .reserve(num_part, this->get_comm());

        size_t local_num, local_offset, file_offset;
        std::tie(local_num, local_offset) =
            this->get_local_particle_count_in_range(
                ParticleGroup::regular, num_part, 0);
        if (MPI_Scan(&local_num,
                     &file_offset,
                     1,
                     MPI_SIZE_T,
                     MPI_SUM,
                     this->get_comm()) != MPI_SUCCESS) {
            std::runtime_error("Error in MPI_Scan in bunch-read-file-openpmd!");
        }

        // drain to reset
        this->get_bunch_particles(PG::regular).drain();

        bunch_particles_t<double>::host_parts_t parts_subset;
        bunch_particles_t<double>::host_masks_t masks_subset;
        parts_subset = bunch_particles_t<double>::host_parts_t(
            "bunch_read_openpmd_" +
                this->get_bunch_particles(PG::regular).get_label(),
            local_num);
        masks_subset = bunch_particles_t<double>::host_masks_t(
            "bunch_read_openpmd_" +
                this->get_bunch_particles(PG::regular).get_label() + "_masks",
            local_num);
        this->get_bunch_particles(PG::regular)
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

        openPMD::Offset chunk_offset = {file_offset - local_num};
        openPMD::Extent chunk_extent = {local_num};

        protons["position"]["x"].loadChunkRaw(
            parts_x.data(), chunk_offset, chunk_extent);
        protons["position"]["y"].loadChunkRaw(
            parts_y.data(), chunk_offset, chunk_extent);
        protons["position"]["z"].loadChunkRaw(
            parts_z.data(), chunk_offset, chunk_extent);
        protons["moments"]["x"].loadChunkRaw(
            moments_x.data(), chunk_offset, chunk_extent);
        protons["moments"]["y"].loadChunkRaw(
            moments_y.data(), chunk_offset, chunk_extent);
        protons["moments"]["z"].loadChunkRaw(
            moments_dpop.data(), chunk_offset, chunk_extent);
        protons["id"][openPMD::RecordComponent::SCALAR].loadChunkRaw(
            parts_idx.data(), chunk_offset, chunk_extent);

        masks["id"][openPMD::RecordComponent::SCALAR].loadChunkRaw(
            masks_subset.data(), chunk_offset, chunk_extent);

        series.flush();

        this->get_bunch_particles(PG::regular)
            .put_particles_in_range(parts_subset,
                                    masks_subset,
                                    local_num,
                                    local_offset,
                                    this->get_comm());

        return;
    }

    else if (software_name == "ImpactX") {

        openPMD::ParticleSpecies& protons = iteration.particles["beam"];
        auto parts_ids = protons["id"][openPMD::RecordComponent::SCALAR];
        auto extent = parts_ids.getExtent();
        //  The extent should be a vector of one element!
        auto num_part = extent[0];
        assert((extent.size() == 1) &&
               "Extent of particle_ids extent should be 1-dimensional!");

        double mass = protons.getAttribute("mass").get<double>();
        // convert mass to GeV!
        mass = mass * pconstants::kg_to_GeV;
        double beta_ref = protons.getAttribute("beta_ref").get<double>();
        double gamma_ref = protons.getAttribute("gamma_ref").get<double>();

        Reference_particle& ref_part = this->get_reference_particle();
        Four_momentum fm = Four_momentum(mass);
        fm.set_beta(beta_ref);
        fm.set_gamma(gamma_ref);
        ref_part.set_four_momentum(fm);

        // reserver required space
        this->get_bunch_particles(PG::regular)
            .reserve(num_part, this->get_comm());

        size_t local_num, local_offset, file_offset;
        std::tie(local_num, local_offset) =
            this->get_local_particle_count_in_range(
                ParticleGroup::regular, num_part, 0);
        if (MPI_Scan(&local_num,
                     &file_offset,
                     1,
                     MPI_SIZE_T,
                     MPI_SUM,
                     this->get_comm()) != MPI_SUCCESS) {
            std::runtime_error("Error in MPI_Scan in bunch-read-file-openpmd!");
        }

        // drain to reset
        this->get_bunch_particles(PG::regular).drain();

        bunch_particles_t<double>::host_parts_t parts_subset;
        bunch_particles_t<double>::host_masks_t masks_subset;
        parts_subset = bunch_particles_t<double>::host_parts_t(
            "bunch_read_openpmd_" +
                this->get_bunch_particles(PG::regular).get_label(),
            local_num);
        masks_subset = bunch_particles_t<double>::host_masks_t(
            "bunch_read_openpmd_" +
                this->get_bunch_particles(PG::regular).get_label() + "_masks",
            local_num);
        this->get_bunch_particles(PG::regular)
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

        Kokkos::View<uint64_t*> idx("bunch_openpmd_read_idx", num_part);

        openPMD::Offset chunk_offset = {file_offset - local_num};
        openPMD::Extent chunk_extent = {local_num};

        protons["position"]["x"].loadChunkRaw(
            parts_x.data(), chunk_offset, chunk_extent);
        protons["position"]["y"].loadChunkRaw(
            parts_y.data(), chunk_offset, chunk_extent);
        protons["position"]["t"].loadChunkRaw(
            parts_z.data(), chunk_offset, chunk_extent);
        protons["momentum"]["x"].loadChunkRaw(
            moments_x.data(), chunk_offset, chunk_extent);
        protons["momentum"]["y"].loadChunkRaw(
            moments_y.data(), chunk_offset, chunk_extent);
        protons["momentum"]["t"].loadChunkRaw(
            moments_dpop.data(), chunk_offset, chunk_extent);
        protons["id"][openPMD::RecordComponent::SCALAR].loadChunkRaw(
            idx.data(), chunk_offset, chunk_extent);

        series.flush();

        // as of Jul 2023, all particles in ImpactX output
        // are valid particles!
        Kokkos::parallel_for(
            "set_particle_mask_ids", local_num, KOKKOS_LAMBDA(const int i) {
                masks_subset[i] = 1;
            });
        // OpenPMD cannot implicitly convert uint64_t to float64, so do
        // that conversion here while transferring input read into
        // temporary idx view into the parts_subset view
        Kokkos::parallel_for(
            "set_particle_ids", local_num, KOKKOS_LAMBDA(const int i) {
                parts_subset(i, 6) = static_cast<double>(idx[i]);
            });

        Kokkos::parallel_for(
            "convert_dpt_to_E", local_num, KOKKOS_LAMBDA(const int i) {
                parts_subset(i, 5) =
                    mass *
                    (gamma_ref - (parts_subset(i, 5) * beta_ref * gamma_ref));
            });
        Kokkos::parallel_for(
            "convert_E_to_dpop", local_num, KOKKOS_LAMBDA(const int i) {
                auto p0 = beta_ref * gamma_ref * mass;

                parts_subset(i, 5) =
                    (
                        /* sqrt (E^2 - m^2) */
                        Kokkos::sqrt((parts_subset(i, 5) * parts_subset(i, 5)) -
                                     (mass * mass)) -
                        /* subtract p0 */
                        p0) *
                    /* scale by p0 */
                    (1 / p0);
            });
        Kokkos::fence();

        this->get_bunch_particles(PG::regular)
            .put_particles_in_range(parts_subset,
                                    masks_subset,
                                    local_num,
                                    local_offset,
                                    this->get_comm());

        return;

    }

    else {
        std::runtime_error("Reading from OpenPMD files is only implemented for "
                           "output generated by synergia3 and ImpactX! \n");
    }
}

// num_part = -1 means write all particles
template <>
void
Bunch::write_openpmd_file(std::string const& filename,
                          int num_part,
                          int offset,
                          int num_part_spec,
                          int offset_spec) const
{
    size_t local_num, local_offset, file_offset;
    size_t spec_local_num, spec_local_offset, spec_file_offset;

    std::tie(local_num, local_offset) = this->get_local_particle_count_in_range(
        ParticleGroup::regular, num_part, offset);

    std::tie(spec_local_num, spec_local_offset) =
        this->get_local_particle_count_in_range(
            ParticleGroup::spectator, num_part_spec, offset_spec);

    if (MPI_Scan(&local_num,
                 &file_offset,
                 1,
                 MPI_SIZE_T,
                 MPI_SUM,
                 this->get_comm()) != MPI_SUCCESS) {
        std::runtime_error("Error in MPI_Scan in bunch-write-openpmd-file!");
    }

    if (MPI_Scan(&spec_local_num,
                 &spec_file_offset,
                 1,
                 MPI_SIZE_T,
                 MPI_SUM,
                 this->get_comm()) != MPI_SUCCESS) {
        std::runtime_error("Error in MPI_Scan in bunch-write-openpmd-file!");
    }

    auto io_device =
        openPMD::Series(filename, openPMD::Access::CREATE, this->get_comm());

    std::string version =
        std::to_string(synergia_version::synergia_version_year);
    version.append(".");
    version.append(std::to_string(synergia_version::synergia_version_month));
    version.append(".");
    version.append(std::to_string(synergia_version::synergia_version_day));
    version.append("-");
    version.append(synergia_version::synergia_git_hash);

    io_device.setSoftware("synergia3", version);

    if (local_num > 0) {
        if (num_part == -1) {
            if (MPI_Allreduce(&local_num,
                              &num_part,
                              1,
                              MPI_INT,
                              MPI_SUM,
                              this->get_comm()) != MPI_SUCCESS) {
                std::runtime_error(
                    "Error in MPI_Allreduce in bunch-write-file-openpmd!");
            };
        }
        auto label = (this->get_bunch_particles(PG::regular)).get_label();

        size_t iteration = 0;
        bunch_particles_t<double>::host_parts_t parts_subset;
        bunch_particles_t<double>::host_masks_t masks_subset;
        parts_subset = bunch_particles_t<double>::host_parts_t(
            "bunch_write_" + label, local_num);
        masks_subset = bunch_particles_t<double>::host_masks_t(
            "bunch_write_" + label + "_masks", local_num);

        openPMD::ParticleSpecies& protons =
            io_device.iterations[iteration].particles["bunch_" + label];
        openPMD::ParticleSpecies& masks =
            io_device.iterations[iteration]
                .particles["bunch_" + label + "_masks"];

        // write mass in SI units!
        protons.setAttribute("mass", this->get_mass() / pconstants::kg_to_GeV);
        protons.setAttribute("beta_ref",
                             (this->get_reference_particle()).get_beta());
        protons.setAttribute("gamma_ref",
                             (this->get_reference_particle()).get_gamma());

        openPMD::Datatype datatype = openPMD::determineDatatype<double>();
        openPMD::Extent global_extent = {static_cast<size_t>(num_part)};
        openPMD::Dataset dataset = openPMD::Dataset(datatype, global_extent);

        openPMD::Datatype masks_datatype =
            openPMD::determineDatatype<uint8_t>();
        openPMD::Dataset masks_dataset =
            openPMD::Dataset(masks_datatype, global_extent);

        this->get_bunch_particles(ParticleGroup::regular)
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

        masks["id"][openPMD::RecordComponent::SCALAR].resetDataset(
            masks_dataset);

        protons["id"][openPMD::RecordComponent::SCALAR].resetDataset(dataset);
        protons["position"]["x"].resetDataset(dataset);
        protons["position"]["y"].resetDataset(dataset);
        protons["position"]["z"].resetDataset(dataset);
        protons["moments"]["x"].resetDataset(dataset);
        protons["moments"]["y"].resetDataset(dataset);
        protons["moments"]["z"].resetDataset(dataset);

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

        masks["id"][openPMD::RecordComponent::SCALAR].storeChunkRaw(
            masks_subset.data(), chunk_offset, chunk_extent);

        io_device.flush();
    }
    return;
}

#endif
