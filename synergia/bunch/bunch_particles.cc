
#include <iomanip>
#include <iostream>

#include "synergia/bunch/bunch_particles.h"
#include "synergia/utils/hdf5_file.h"
#include "synergia/utils/parallel_utils.h"

// static init
namespace bunch_particles_impl {
    int pid_offset::offset = 0;
}

namespace {
    // Kokkos functors

    struct particle_copier_many_row {
        ConstParticles src;
        karray2d_row_dev dst;

        ConstParticleMasks s_masks;
        ParticleMasks d_masks;

        int idx;

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            dst(i, 0) = src(idx + i, 0);
            dst(i, 1) = src(idx + i, 1);
            dst(i, 2) = src(idx + i, 2);
            dst(i, 3) = src(idx + i, 3);
            dst(i, 4) = src(idx + i, 4);
            dst(i, 5) = src(idx + i, 5);
            dst(i, 6) = src(idx + i, 6);

            d_masks(i) = s_masks(idx + i);
        }
    };

    struct particle_copier_one {
        ConstParticles src;
        karray1d_row_dev dst;

        ConstParticleMasks s_masks;
        ParticleMasks d_masks;

        int idx;

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            dst(0) = src(idx + i, 0);
            dst(1) = src(idx + i, 1);
            dst(2) = src(idx + i, 2);
            dst(3) = src(idx + i, 3);
            dst(4) = src(idx + i, 4);
            dst(5) = src(idx + i, 5);
            dst(6) = src(idx + i, 6);

            d_masks(0) = s_masks(idx + i);
        }
    };

    struct particle_id_checker {
        ConstParticles parts;
        int idx;
        int pid;

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i, int& match) const
        {
            match = (((int)parts(idx, 6)) == pid);
        }
    };

    struct particle_finder {
        ConstParticles parts;
        int pid;

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i, int& idx) const
        {
            if (((int)parts(i, 6)) == pid) idx = i;
        }
    };

    struct particle_zeroer {
        Particles parts;
        int offset;

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            for (int j = 0; j < 7; ++j)
                parts(offset + i, j) = 0.0;
        }
    };

    struct masks_zeroer {
        ParticleMasks masks;
        int offset;

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            masks(offset + i) = 0.0;
        }
    };

    struct mask_reducer {
        ConstParticleMasks masks;

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i, int& valid) const
        {
            if (masks(i)) ++valid;
        }
    };

    struct particle_injector {
        Particles dst;
        ParticleMasks dst_masks;

        ConstParticles src;
        ConstParticleMasks src_masks;

        karray1d_dev ref_st_diff;
        karray1d_dev tgt_st;
        karray1d_dev inj_st;

        int off;
        double pdiff;

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            // space-like coordinates
            dst(off + i, 0) = src(i, 0); // + ref_st_diff(0);
            dst(off + i, 2) = src(i, 2); // + ref_st_diff(2);
            dst(off + i, 4) = src(i, 4); // + ref_st_diff(4);

            // npx and npy coordinates are scaled with p_ref which can be
            // different for different bunches
            dst(off + i, 1) = pdiff * (src(i, 1) - inj_st(1)) + tgt_st(1);
            dst(off + i, 3) = pdiff * (src(i, 3) - inj_st(3)) + tgt_st(3);

            // ndp coordinate is delta-p scaled with pref
            dst(off + i, 5) =
                pdiff * (1.0 + src(i, 5) - inj_st(5)) + tgt_st(5) - 1.0;

            // particle id
            dst(off + i, 6) = src(i, 6);

            dst_masks(off + i) = src_masks(i);
        }
    };

    struct fixed_z_to_t_converter {
        Particles parts;
        ConstParticleMasks masks;

        double p_ref;
        double beta;

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            if (masks(i)) {
                double p = p_ref + parts(i, 5) * p_ref;
                double px = parts(i, 1) * p_ref;
                double py = parts(i, 3) * p_ref;
                double pz2 = p * p - px * px - py * py;
                double pz = sqrt(pz2);

                parts(i, 4) = -parts(i, 4) * beta;
                parts(i, 5) = pz / p_ref;
            }
        }
    };

    struct fixed_t_to_z_converter {
        Particles parts;
        ConstParticleMasks masks;

        double p_ref;
        double beta;

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            if (masks(i)) {
                double px = parts(i, 1) * p_ref;
                double py = parts(i, 3) * p_ref;
                double pz = parts(i, 5) * p_ref;
                double p = sqrt(px * px + py * py + pz * pz);

                parts(i, 4) = -parts(i, 4) / beta;
                parts(i, 5) = (p - p_ref) / p_ref;
            }
        }
    };

    struct discarded_particle_mover {
        Kokkos::View<int*> counter;

        ConstParticles parts;
        ParticleMasks discards;
        karray2d_row_dev discarded_parts;

        discarded_particle_mover(ConstParticles const& parts,
                                 ParticleMasks const& discards,
                                 karray2d_row_dev const& discarded_parts)
            : counter("counter", 1)
            , parts(parts)
            , discards(discards)
            , discarded_parts(discarded_parts)
        {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            if (discards(i)) {
                int pos = Kokkos::atomic_fetch_add(counter.data(), 1);

                discarded_parts(pos, 0) = parts(i, 0);
                discarded_parts(pos, 1) = parts(i, 1);
                discarded_parts(pos, 2) = parts(i, 2);
                discarded_parts(pos, 3) = parts(i, 3);
                discarded_parts(pos, 4) = parts(i, 4);
                discarded_parts(pos, 5) = parts(i, 5);
                discarded_parts(pos, 6) = parts(i, 6);
            }
        }
    };
}

template <>
void
bunch_particles_t<double>::reserve_local(int r)
{
    if (r <= n_reserved) return;

#ifdef NO_PADDING
    Kokkos::resize(parts, r);
#else
    Kokkos::resize(Kokkos::AllowPadding, parts, r);
#endif
    n_reserved = parts.stride(1);

    Kokkos::resize(masks, n_reserved);
    Kokkos::resize(discards, n_reserved);

    hparts = Kokkos::create_mirror_view(parts);
    hmasks = Kokkos::create_mirror_view(masks);
    hdiscards = Kokkos::create_mirror_view(discards);

    // update n_reserved to reflect the increased size of particle views
    n_reserved = r;
}

template <>
void
bunch_particles_t<double>::reserve(int n, Commxx const& comm)
{
    int r = decompose_1d_local(comm, n);
    reserve_local(r);
}

template <>
void
bunch_particles_t<double>::assign_ids(int train_idx, int bunch_idx)
{
    // each bunch is assined a range in the global id space
    //
    // here we assume the max number of particles that a single
    // bunch can hold is 2^32.
    // a bunch can have 2^2, or 4 groups of particles
    // a train can have 2^16, or 65536 bunches.
    // a bunch simulator can have 2^2, or 4 trains
    //
    // so the overall bits the ids take is 32+2+16+2 = 52, less
    // than 53 which is the max integer number that can be
    // represented in double

    int64_t base = 0;
    base |= (int64_t)train_idx << 50;
    base |= (int64_t)bunch_idx << 34;
    base |= (int64_t)group << 32;

    // particle_id_assigner defined in header
    using namespace bunch_particles_impl;

    pid_assigner<parts_t> pia{parts, base + poffset};
    Kokkos::parallel_for(n_active, pia);
}

template <>
void
bunch_particles_t<double>::drain()
{
    // reset the particle array
    particle_zeroer pz{parts, 0};
    masks_zeroer mz{masks, 0};

    Kokkos::parallel_for(n_reserved, pz);
    Kokkos::parallel_for(n_reserved, mz);

    // reset the pointers
    n_valid = 0;
    n_total = 0;
    n_active = 0;

    return;
}

template <>
void
bunch_particles_t<double>::inject(bunch_particles_t const& o,
                                  karray1d_dev const& ref_st_diff,
                                  karray1d_dev const& tgt_st,
                                  karray1d_dev const& inj_st,
                                  double pdiff)
{
    // no nothing for empty inj bunch
    if (!o.n_active) return;

    // expand the particle array
    reserve_local(n_active + o.n_active);

    // inject
    particle_injector pi{parts,
                         masks,
                         o.parts,
                         o.masks,
                         ref_st_diff,
                         tgt_st,
                         inj_st,
                         n_active,
                         pdiff};

    Kokkos::parallel_for(o.n_active, pi);

    // update number of particles
    n_active += o.n_active;
    n_valid += o.n_valid;
}

template <>
void
bunch_particles_t<double>::convert_to_fixed_t_lab(double p_ref, double beta)
{
    fixed_z_to_t_converter alg{parts, masks, p_ref, beta};
    Kokkos::parallel_for(n_active, alg);
}

template <>
void
bunch_particles_t<double>::convert_to_fixed_z_lab(double p_ref, double beta)
{
    fixed_t_to_z_converter alg{parts, masks, p_ref, beta};
    Kokkos::parallel_for(n_active, alg);
}

template <>
int
bunch_particles_t<double>::update_valid_num()
{
    int old_valid_num = n_valid;
    mask_reducer mr{masks};
    Kokkos::parallel_reduce(n_active, mr, n_valid);
    return old_valid_num;
}

template <>
int
bunch_particles_t<double>::update_total_num(Commxx const& comm)
{
    int old_total_num = n_total;
    if (MPI_Allreduce(&n_valid, &n_total, 1, MPI_INT, MPI_SUM, comm) !=
        MPI_SUCCESS) {
        std::runtime_error(
            "Error in MPI_Allreduce in bunch_particles::update_total_num!");
    }
    return old_total_num;
}

template <>
void
bunch_particles_t<double>::put_particles_in_range(host_parts_t subset_parts,
                                                  host_masks_t subset_masks,
                                                  size_t local_num,
                                                  size_t local_idx,
                                                  Commxx const& comm)
{
    // index out of range
    if (local_idx == particle_index_null ||
        local_idx + local_num > this->n_reserved)
        throw std::runtime_error("Bunch::get_particle() index out of range");

    if (subset_parts.size() != local_num * 7) {
        std::cerr << "subset_parts size is : " << subset_parts.size()
                  << ", local_num * 7 is : " << local_num * 7 << "\n";
        throw std::runtime_error("Bunch::put_particles_in_range() subset_parts "
                                 "has inconsistent size!");
    }
    if (subset_masks.size() != local_num) {
        std::cerr << "subset_masks size is : " << subset_masks.size()
                  << ", local_num is : " << local_num << "\n";
        throw std::runtime_error("Bunch::put_particles_in_range() subset_masks "
                                 "has inconsistent size!");
    }

    ParticlesSubView parts_in_range =
        Kokkos::subview(parts,
                        Kokkos::make_pair(local_idx, local_idx + local_num),
                        Kokkos::ALL);
    ParticleMasksSubView masks_in_range = Kokkos::subview(
        masks, Kokkos::make_pair(local_idx, local_idx + local_num));

    Kokkos::deep_copy(parts_in_range, subset_parts);
    Kokkos::deep_copy(masks_in_range, subset_masks);

    // update num active
    this->n_active += local_num;

    // update local num from parts
    this->update_valid_num();

    // now we have the actual local num, update the total number
    this->update_total_num(comm);

    return;
}

template <>
void
bunch_particles_t<double>::get_particles_in_range(host_parts_t subset_parts,
                                                  host_masks_t subset_masks,
                                                  size_t local_num,
                                                  size_t local_idx) const
{
    // index out of range
    if (local_idx == particle_index_null ||
        local_idx + local_num > this->n_reserved)
        throw std::runtime_error("Bunch::get_particle() index out of range");

    ParticlesSubView parts_in_range =
        Kokkos::subview(parts,
                        Kokkos::make_pair(local_idx, local_idx + local_num),
                        Kokkos::ALL);
    ParticleMasksSubView masks_in_range = Kokkos::subview(
        masks, Kokkos::make_pair(local_idx, local_idx + local_num));

    for (int i = 0; i < 7; i++) {
        auto parts_col_dev = Kokkos::subview(parts_in_range, Kokkos::ALL, i);
        auto parts_col_host = Kokkos::subview(subset_parts, Kokkos::ALL, i);
        Kokkos::deep_copy(parts_col_host, parts_col_dev);
    }
    Kokkos::deep_copy(subset_masks, masks_in_range);

    return;
}

template <>
std::pair<karray2d_row, HostParticleMasks>
bunch_particles_t<double>::get_particles_in_range_row(int idx, int n) const
{
    // index out of range
    if (idx == particle_index_null || idx < 0 || idx + n > n_active)
        throw std::runtime_error("Bunch::get_particle() index out of range");

    karray2d_row_dev p("sub_p", n, 7);
    ParticleMasks pm("masks", n);

    particle_copier_many_row pc{parts, p, masks, pm, idx};
    Kokkos::parallel_for(n, pc);

    karray2d_row hp = create_mirror_view(p);
    Kokkos::deep_copy(hp, p);

    HostParticleMasks hpm = create_mirror_view(pm);
    Kokkos::deep_copy(hpm, pm);

    return std::make_pair(hp, hpm);
}

template <>
std::pair<karray1d_row, bool>
bunch_particles_t<double>::get_particle(int idx) const
{
    // index out of range
    if (idx == particle_index_null || idx < 0 || idx > n_active)
        throw std::runtime_error("Bunch::get_particle() index out of range");

    karray1d_row_dev p("particle", 7);
    ParticleMasks pm("mask", 1);

    particle_copier_one pc{parts, p, masks, pm, idx};
    Kokkos::parallel_for(1, pc);

    karray1d_row hp = create_mirror_view(p);
    Kokkos::deep_copy(hp, p);

    HostParticleMasks hpm = create_mirror_view(pm);
    Kokkos::deep_copy(hpm, pm);

    return std::make_pair(hp, hpm(0));
}

template <>
int
bunch_particles_t<double>::search_particle(int pid, int last_idx) const
{
    if (last_idx != particle_index_null) {
        int match = 0;
        particle_id_checker pic{parts, last_idx, pid};
        Kokkos::parallel_reduce(1, pic, match);

        if (match) return last_idx;
    }

    int idx = particle_index_null;
    particle_finder pf{parts, pid};
    Kokkos::parallel_reduce(n_active, pf, idx);

    return idx;
}

template <>
karray2d_row
bunch_particles_t<double>::get_particles_last_discarded() const
{
    karray2d_row_dev discarded("discarded", n_last_discarded, 7);
    karray2d_row hdiscarded = Kokkos::create_mirror_view(discarded);

    discarded_particle_mover dpm(parts, discards, discarded);
    Kokkos::parallel_for(n_active, dpm);

    Kokkos::deep_copy(hdiscarded, discarded);
    return hdiscarded;
}

template <>
void
bunch_particles_t<double>::check_pz2_positive()
{
    checkout_particles();

    for (int p = 0; p < n_active; ++p) {
        if (!hmasks(p)) continue;

        double pzop2 = (1. + hparts(p, 5)) * (1. + hparts(p, 5)) -
                       hparts(p, 1) * hparts(p, 1) -
                       hparts(p, 3) * hparts(p, 3);

        if (pzop2 < 0.0) {
            std::cout << "pzop^2 = " << pzop2 << std::endl;
            throw std::runtime_error(
                " check pz2:  pz square cannot be negative!");
        }
    }
}

template <>
void
bunch_particles_t<double>::read_file_legacy(Hdf5_file const& file,
                                            Commxx const& comm)
{
    auto dims = file.get_dims("particles");
    if (dims.size() != 2 || dims[1] != 7) {
        throw std::runtime_error(
            "BunchParticle::read_file_legacy(): wrong data dimensions in file");
    }

    int file_total = dims[0];
    int file_num = decompose_1d_local(comm, file_total);

    // size check
    if (n_active != file_num) {
        throw std::runtime_error(
            " local_num incompatibility when initializing the bunch");
    }

    // read
    auto read_particles = file.read<karray2d_row>("particles", n_active);

    // transpose: read_particles is row major, hparts is col major
    for (int part = 0; part < n_active; ++part)
        for (int i = 0; i < 7; ++i)
            hparts(part, i) = read_particles(part, i);

    // check in to device mem
    checkin_particles();

    // all particles in the legacy file are valid
    // do the init after checkin because masks_initializer is performed
    // on the device memory
    n_valid = n_active;

    // particle_masks_initializer defined in header
    using namespace bunch_particles_impl;

    particle_masks_initializer<masks_t> pmi{masks, n_valid};
    Kokkos::parallel_for(n_reserved, pmi);
}

template <>
void
bunch_particles_t<double>::read_file(Hdf5_file const& file, Commxx const& comm)
{
    auto dims = file.get_dims(label);
    if (dims.size() != 2 || dims[1] != 7) {
        throw std::runtime_error(
            "BunchParticle::read_file(): wrong data dimensions in file");
    }

    int file_total = dims[0];
    int file_num = decompose_1d_local(comm, file_total);

    // allocation if capacity is smaller
    if (n_reserved < file_num) {
#ifdef NO_PADDING
        auto alloc = Kokkos::view_alloc(label);
#else
        auto alloc = Kokkos::view_alloc(label, Kokkos::AllowPadding);
#endif
        parts = Particles(alloc, file_num);
        n_reserved = parts.stride(1);

        masks = ParticleMasks(label + "_masks", n_reserved);
        discards = ParticleMasks(label + "_discards", n_reserved);

        hparts = Kokkos::create_mirror_view(parts);
        hmasks = Kokkos::create_mirror_view(masks);
        hdiscards = Kokkos::create_mirror_view(discards);
    }

    // reset the pointers
    n_valid = 0;
    n_total = 0;
    n_active = file_num;

    // read from file
    auto read_particles = file.read<karray2d_row>(label, file_num);
    auto read_masks = file.read<HostParticleMasks>(label + "_masks", file_num);

    // transpose: read_particles is row major, hparts is col major
    for (int part = 0; part < file_num; ++part) {
        for (int i = 0; i < 7; ++i)
            hparts(part, i) = read_particles(part, i);

        hmasks(part) = read_masks(part);
        if (hmasks(part)) ++n_valid;
    }

    // now we have the actual local num, update the total number
    update_total_num(comm);

    // check in to device mem
    checkin_particles();
}

template <>
void
bunch_particles_t<double>::write_file(Hdf5_file const& file,
                                      int num_part,
                                      int offset,
                                      Commxx const& comm) const
{
    int local_num_part = 0;
    int local_offset = 0;

    if (num_part == -1) {
        local_num_part = n_active;
        local_offset = 0;
    } else {
        local_num_part = decompose_1d_local(comm, num_part);
        local_offset = decompose_1d_local(comm, offset);
    }

    if (local_num_part < 0 || local_offset < 0 ||
        local_num_part + local_offset > n_active) {
        throw std::runtime_error(
            "invalid num_part or offset for bunch_particles_t::write_file()");
    }

    auto parts = get_particles_in_range_row(local_offset, local_num_part);
    file.write_collective(label, parts.first);
    file.write_collective(label + "_masks", parts.second);
}

template <>
void
bunch_particles_t<double>::print_particle(size_t idx, Logger& logger) const
{
    logger(LoggerV::DEBUG) << std::showpos << std::scientific
                           << std::setprecision(8) << std::setw(12)
                           << hparts(idx, 0) << ", " << std::setw(12)
                           << hparts(idx, 1) << ", " << std::setw(12)
                           << hparts(idx, 2) << ", " << std::setw(12)
                           << hparts(idx, 3) << ", " << std::setw(12)
                           << hparts(idx, 4) << ", " << std::setw(12)
                           << hparts(idx, 5) << "\n"
                           << std::defaultfloat << std::noshowpos;
}

template <>
void
bunch_particles_t<double>::save_checkpoint_particles(Hdf5_file& file, int idx)
{
    checkout_particles();

    std::stringstream ss;
    ss << "bunch_particles_" << label << "_parts_" << idx;
    file.write(ss.str(), hparts.data(), hparts.span(), true);

    ss.str("");
    ss << "bunch_particles_" << label << "_masks_" << idx;
    file.write(ss.str(), hmasks.data(), hmasks.span(), true);
}

template <>
void
bunch_particles_t<double>::load_checkpoint_particles(Hdf5_file& file, int idx)
{
    std::stringstream ss;
    ss << "bunch_particles_" << label << "_parts_" << idx;
    file.read(ss.str(), hparts.data(), hparts.span());

    ss.str("");
    ss << "bunch_particles_" << label << "_masks_" << idx;
    file.read(ss.str(), hmasks.data(), hmasks.span());

    checkin_particles();
}
