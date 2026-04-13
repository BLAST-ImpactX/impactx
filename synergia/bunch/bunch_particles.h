#ifndef BUNCH_PARTICLES_H
#define BUNCH_PARTICLES_H

#include "synergia/foundation/trigon_traits.h"
#include "synergia/utils/commxx.h"
#include "synergia/utils/gsvector.h"
#include "synergia/utils/hdf5_file.h"
#include "synergia/utils/kokkos_views.h"
#include "synergia/utils/logger.h"

#include <cereal/cereal.hpp>
#include <string>
#include <utility>

#if defined SYNERGIA_HAVE_OPENPMD
#include <openPMD/openPMD.hpp>
#endif

// To keep track of the memory location for the particles
enum class MemoryLocation { Host, Device };

enum class ParticleGroup { regular = 0, spectator = 1 };

using Particles = Kokkos::View<double* [7],
                               Kokkos::LayoutLeft,
                               Kokkos::DefaultExecutionSpace::memory_space>;
using ConstParticles =
    Kokkos::View<const double* [7],
                 Kokkos::LayoutLeft,
                 Kokkos::DefaultExecutionSpace::memory_space>;

using ParticlesSubView = decltype(Kokkos::subview(std::declval<Particles>(),
                                                  std::pair(1u, 1u),
                                                  Kokkos::ALL));
using ConstParticlesSubView =
    decltype(Kokkos::subview(std::declval<ConstParticles>(),
                             std::pair(1u, 1u),
                             Kokkos::ALL));

using ParticleMasks =
    Kokkos::View<uint8_t*, Kokkos::DefaultExecutionSpace::memory_space>;
using ConstParticleMasks =
    Kokkos::View<const uint8_t*, Kokkos::DefaultExecutionSpace::memory_space>;

using ParticleMasksSubView =
    decltype(Kokkos::subview(std::declval<ParticleMasks>(), std::pair(1u, 1u)));
using ConstParticleMasksSubView =
    decltype(Kokkos::subview(std::declval<ConstParticleMasks>(),
                             std::pair(1u, 1u)));

using HostParticles = Particles::HostMirror;
using ConstHostParticles = ConstParticles::HostMirror;
using HostParticlesSubView = ParticlesSubView::HostMirror;
using ConstHostParticlesSubView = ConstParticlesSubView::HostMirror;

using HostParticleMasks = ParticleMasks::HostMirror;
using ConstHostParticleMasks = ConstParticleMasks::HostMirror;
using HostParticleMasksSubView = ParticleMasksSubView::HostMirror;
using ConstHostParticleMasksSubView = ConstParticleMasksSubView::HostMirror;

// serialization
namespace cereal {
    // particles
    template <class AR>
    void
    save(AR& ar, Particles const& p)
    {
        ar(p.label(), p.stride(1));
    }

    template <class AR>
    void
    load(AR& ar, Particles& p)
    {
        std::string label;
        int slots;
        ar(label, slots);

#ifdef NO_PADDING
        auto alloc = Kokkos::view_alloc(label);
#else
        auto alloc = Kokkos::view_alloc(label, Kokkos::AllowPadding);
#endif

        p = Particles(alloc, slots);

        if (p.stride(1) != slots)
            throw std::runtime_error(
                "inconsistent padding while loading particles");
    }

    // masks
    template <class AR>
    void
    save(AR& ar, ParticleMasks const& p)
    {
        ar(p.label(), p.extent(0));
    }

    template <class AR>
    void
    load(AR& ar, ParticleMasks& p)
    {
        std::string label;
        int slots;

        ar(label, slots);
        p = ParticleMasks(label, slots);
    }
}

template <class PART>
class bunch_particles_t {
    using PG = ParticleGroup;

  public:
    constexpr static const int particle_index_null = -1;

    using part_t = PART;

    using default_memspace =
        typename Kokkos::DefaultExecutionSpace::memory_space;

    using memspace = typename std::conditional<is_trigon<PART>::value,
                                               Kokkos::HostSpace,
                                               default_memspace>::type;

    using parts_t =
        typename Kokkos::View<PART* [7], Kokkos::LayoutLeft, memspace>;

    using masks_t = typename Kokkos::View<uint8_t*, memspace>;

    using host_parts_t = typename parts_t::HostMirror;

    using host_masks_t = typename masks_t::HostMirror;

    using const_parts_t =
        typename Kokkos::View<const PART* [7], Kokkos::LayoutLeft, memspace>;

    using const_masks_t = typename Kokkos::View<const uint8_t*, memspace>;

    using const_host_parts_t = typename const_parts_t::HostMirror;

    using const_host_masks_t = typename const_masks_t::HostMirror;

    using gsv_t = typename std::
        conditional<is_trigon<PART>::value, GSVec<PART>, GSVector>::type;

    using exec_space = typename parts_t::execution_space;

  private:
    /*
     * Local Particle Array Memory Layout:
     *
     *   P: regular particle,
     *   I: invalid (lost) particle
     *   R: reserved particle
     *
     *    part       mask
     *   +=====+    +=====+
     *   |  P  |    |  1  |
     *   +-----+    +-----+
     *   |  I  |    |  0  |
     *   +-----+    +-----+
     *   |  P  |    |  1  |
     *   +-----+    +-----+
     *   |  P  |    |  1  |
     *   +-----+    +-----+
     *   |  P  |    |  1  |
     *   +-----+    +-----+
     *   |  I  |    |  0  |
     *   +-----+    +-----+
     *   |  I  |    |  0  |
     *   +-----+    +-----+
     *   |  P  |    |  1  |
     *   +-----+    +-----+
     *   |  R  |    |  0  |
     *   +-----+    +-----+
     *   |  R  |    |  0  |
     *   +-----+    +-----+
     *   |  R  |    |  0  |
     *   +-----+    +-----+
     *   |  R  |    |  0  |
     *   +=====+    +=====+
     *
     *   num_valid    = 5   -- num_valid()
     *   num_active   = 8   -- num_active() / size()
     *   num_reserved = 12  -- num_reserved() / capacity()
     *
     *   when allocating an particle array, the padding flag of the
     *   Kokkos::View object is turned on. So the 'num_reserved' is
     *   the actual array size in the first dimension.
     *
     *   'num_valid' is the number of remainig local particles in
     *   the bunch
     *
     *   'num_active' is the number of local particles including lost
     *   ones. This one should be used when looping through the
     *   particles array
     *
     *   'memory_location' indicates whether the particles are currently
     *   on the host or on the device.
     *
     */

    // particle group (regular or spectator)
    ParticleGroup group;
    std::string label;

    // see the memory layout
    int n_valid;
    int n_active;
    int n_reserved;
    MemoryLocation memory_location;

    // sum of n_valid on all ranks
    int n_total;

    // number of particles discarded from the most recent aperture apply.
    int n_last_discarded;

    // particle offset in cases where a single bunch span across
    // multiple ranks
    int poffset;

  public:
    parts_t parts;
    masks_t masks;
    masks_t discards;

    host_parts_t hparts;
    host_masks_t hmasks;
    host_masks_t hdiscards;

  public:
    // particles array will be allocated to the size of reserved_num
    // if reserved is less than total, it will be set to the total_num
    bunch_particles_t(ParticleGroup pg,
                      int total_num,
                      int reserved_num,
                      Commxx const& comm);

    int
    num_valid() const
    {
        return n_valid;
    }
    int
    num_active() const
    {
        return n_active;
    }
    int
    num_reserved() const
    {
        return n_reserved;
    }
    int
    num_total() const
    {
        return n_total;
    }
    int
    num_last_discarded() const
    {
        return n_last_discarded;
    }

    MemoryLocation
    get_memory_location() const
    {
        return memory_location;
    }

    std::pair<size_t, size_t>
    get_local_particle_count_in_range(int num_part, int offset) const
    {
        size_t local_num_part = 0;
        size_t local_offset = 0;
        int n_active = this->num_active();

        if (num_part == -1) {
            local_num_part = n_active;
            local_offset = 0;
        } else {
            local_num_part = decompose_1d_local(*(this->comm), num_part);

            local_offset = decompose_1d_local(*(this->comm), offset);
        }

        if (local_num_part < 0 || local_offset < 0 ||
            local_num_part + local_offset > n_active) {
            throw std::runtime_error("invalid num_part or offset for "
                                     "bunch_particles_t::write_file()");
        }

        return std::make_pair(local_num_part, local_offset);
    }

    // getters with names more consistent with std containers
    int
    size() const
    {
        return n_active;
    }
    int
    capacity() const
    {
        return n_reserved;
    }
    int
    size_in_gsv() const
    {
        return ceil(1.0 * n_active / gsv_t::size());
    }

    // copy particles/masks between host and device memories
    void
    checkin_particles()
    {
        Kokkos::deep_copy(parts, hparts);
        Kokkos::deep_copy(masks, hmasks);
        memory_location = MemoryLocation::Device;
    }

    void
    checkout_particles()
    {
        Kokkos::deep_copy(hparts, parts);
        Kokkos::deep_copy(hmasks, masks);
        memory_location = MemoryLocation::Host;
    }

    // change capacity (can only increase)
    // n is the new cap across the entire bunch
    void reserve(int n, Commxx const& comm);

    // change local capacity (can only increase)
    // n is the new local capacity
    void reserve_local(int n);

    // drain all particles and set n_active/n_valid to 0
    // n_reserved does not change
    void drain();

    // inject with
    void inject(bunch_particles_t const& o,
                karray1d_dev const& ref_st_diff,
                karray1d_dev const& tgt_st,
                karray1d_dev const& inj_st,
                double pdiff);

    // convert between fixed z lab and fixed t lab
    void convert_to_fixed_t_lab(double p_ref, double beta);
    void convert_to_fixed_z_lab(double p_ref, double beta);

    // update the valid num from the masks and return the old valid num
    int update_valid_num();

    // update total num across the ranks and returns the old total number
    int update_total_num(Commxx const& comm);

    // apply aperture operation
    template <typename AP>
    int apply_aperture(AP const& ap);

    // search/get particle(s)
    int search_particle(int pid, int last_idx) const;

    std::pair<karray1d_row, bool> get_particle(int idx) const;

    // hostview is allocated by the caller and copied from by this routine
    // this routine internally calls update valid and total num to
    // ensure that the state is consistent!
    void put_particles_in_range(host_parts_t subset_parts,
                                host_masks_t subset_masks,
                                size_t local_num,
                                size_t local_offset,
                                Commxx const& comm);

    // hostview is allocated by the caller and filled by this routine
    void get_particles_in_range(host_parts_t subset_parts,
                                host_masks_t subset_masks,
                                size_t num,
                                size_t offset) const;

    std::pair<karray2d_row, HostParticleMasks> get_particles_in_range_row(
        int idx,
        int num) const;

    karray2d_row get_particles_last_discarded() const;

    void check_pz2_positive();
    void print_particle(size_t idx, Logger& logger) const;

    std::string
    get_label() const
    {
        return this->label;
    };

    // read from a hdf5 file. total_num of current bunch must be the same
    // as the one stored in the particle file
    void read_file_legacy(Hdf5_file const& file, Commxx const& comm);

    void read_file(Hdf5_file const& file, Commxx const& comm);
    void write_file(Hdf5_file const& file,
                    int num_part,
                    int offset,
                    Commxx const& comm) const;

    // checkpoint save/load
    void save_checkpoint_particles(Hdf5_file& file, int idx);
    void load_checkpoint_particles(Hdf5_file& file, int idx);

    // assign ids cooperatively
    void assign_ids(int train_idx, int bunch_idx);

    // only available for trigons
    template <class U = PART>
    karray2d_row
    get_jacobian(int idx) const
    {
        static_assert(is_trigon<U>::value,
                      "X::get_jacoabian<T> is only defined when U is a trigon");

        karray2d_row res("jacobian", 6, 6);
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j)
                res(i, j) = hparts(idx, i).template get_subpower<1>().terms[j];
        return res;
    }

  private:
    void default_ids(int local_offset, Commxx const& comm);

    // serialization
    friend class cereal::access;

    template <class AR>
    void
    save(AR& ar) const
    {
        ar(CEREAL_NVP(label));
        ar(CEREAL_NVP(memory_location));
        ar(CEREAL_NVP(n_valid));
        ar(CEREAL_NVP(n_active));
        ar(CEREAL_NVP(n_reserved));
        ar(CEREAL_NVP(n_total));

        ar(CEREAL_NVP(parts));
        ar(CEREAL_NVP(masks));
        ar(CEREAL_NVP(discards));
    }

    template <class AR>
    void
    load(AR& ar)
    {
        ar(CEREAL_NVP(label));
        ar(CEREAL_NVP(memory_location));
        ar(CEREAL_NVP(n_valid));
        ar(CEREAL_NVP(n_active));
        ar(CEREAL_NVP(n_reserved));
        ar(CEREAL_NVP(n_total));

        ar(CEREAL_NVP(parts));
        ar(CEREAL_NVP(masks));
        ar(CEREAL_NVP(discards));

        // construct the host array first
        hparts = Kokkos::create_mirror_view(parts);
        hmasks = Kokkos::create_mirror_view(masks);
        hdiscards = Kokkos::create_mirror_view(discards);
    }
};

// instantiated with a double is still the default BunchParticles
typedef bunch_particles_t<double> BunchParticles;

// implementations
namespace bunch_particles_impl {
    template <class AP>
    struct discard_applier {
        typedef int value_type;

        AP ap;
        ConstParticles parts;
        ParticleMasks masks;
        ParticleMasks discards;

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i, int& discarded) const
        {
            discards(i) = 0;

            if (masks(i) && ap.discard(parts, masks, i)) {
                discards(i) = 1;
                masks(i) = 0;
                ++discarded;
            }
        }
    };
}

template <>
template <typename AP>
inline int
bunch_particles_t<double>::apply_aperture(AP const& ap)
{
    using namespace bunch_particles_impl;

    int ndiscarded = 0;

    // go through each particle to see which one is been filtered out
    discard_applier<AP> da{ap, parts, masks, discards};
    Kokkos::parallel_reduce(n_active, da, ndiscarded);

    // std::cout << "      discarded = " << ndiscarded << "\n";
    n_last_discarded = ndiscarded;
    n_valid -= ndiscarded;

    return ndiscarded;
}

#include "synergia/utils/parallel_utils.h"

namespace bunch_particles_impl {
    struct pid_offset {
        static int offset;

        static int
        get(int request_num, Commxx const& comm)
        {
            MPI_Bcast((void*)&offset, 1, MPI_INT, 0, comm);
            int old_offset = offset;
            int total_num;
            MPI_Reduce((void*)&request_num,
                       (void*)&total_num,
                       1,
                       MPI_INT,
                       MPI_SUM,
                       0,
                       comm);
            offset += total_num;
            return old_offset;
        }
    };

    template <typename parts_t>
    struct pid_assigner {
        parts_t parts;
        int64_t offset;

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            parts(i, 6) = i + offset;
        }
    };

    template <typename masks_t>
    struct particle_masks_initializer {
        masks_t masks;
        const int num;

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            masks(i) = i < num ? 1 : 0;
        }
    };
}

// default ids only for double typed bunche_particle object
template <typename PART>
inline void
bunch_particles_t<PART>::default_ids(int local_offset, Commxx const& comm)
{
    using namespace bunch_particles_impl;

    int request_num = (comm.rank() == 0) ? n_total : 0;
    int global_offset = pid_offset::get(request_num, comm);

    auto range = Kokkos::RangePolicy<exec_space>(0, n_active);
    pid_assigner<parts_t> pia{parts, local_offset + global_offset};
    Kokkos::parallel_for(range, pia);
}

template <typename PART>
inline bunch_particles_t<PART>::bunch_particles_t(ParticleGroup pg,
                                                  int total,
                                                  int reserved,
                                                  Commxx const& comm)
    : group(pg)
    , label(pg == PG::regular ? "particles" : "spectators")
    , memory_location(MemoryLocation::Host)
    , n_valid(0)
    , n_active(0)
    , n_reserved(0)
    , n_total(total)
    , n_last_discarded(0)
    , poffset(0)
    , parts()
    , masks()
    , discards()
    , hparts(Kokkos::create_mirror_view(parts))
    , hmasks(Kokkos::create_mirror_view(masks))
    , hdiscards(Kokkos::create_mirror_view(discards))
{
    using namespace bunch_particles_impl;

    // minimum reserved
    if (reserved < total) reserved = total;

    if (!comm.is_null() && reserved) {
        int mpi_size = comm.size();
        int mpi_rank = comm.rank();

        std::vector<int> offsets_t(mpi_size);
        std::vector<int> counts_t(mpi_size);
        decompose_1d(comm, total, offsets_t, counts_t);

        std::vector<int> offsets_r(mpi_size);
        std::vector<int> counts_r(mpi_size);
        decompose_1d(comm, reserved, offsets_r, counts_r);

        // local_num
        n_valid = counts_t[mpi_rank];
        n_active = counts_t[mpi_rank];
        n_reserved = counts_r[mpi_rank];

        if (n_active % gsv_t::size()) {
            int padded = n_active + gsv_t::size() - n_active % gsv_t::size();
            if (n_reserved < padded) n_reserved = padded;
        }

        // local_num offset
        for (int i = 0; i < mpi_rank; ++i)
            poffset += counts_t[i];

            // allocate
#ifdef NO_PADDING
        auto alloc = Kokkos::view_alloc(label);
#else
        auto alloc = Kokkos::view_alloc(label, Kokkos::AllowPadding);
#endif
        parts = parts_t(alloc, n_reserved);

        // with possible paddings
        n_reserved = parts.stride(1);

        masks = masks_t(label + "_masks", n_reserved);
        discards = masks_t(label + "_discards", n_reserved);

        hparts = Kokkos::create_mirror_view(parts);
        hmasks = Kokkos::create_mirror_view(masks);
        hdiscards = Kokkos::create_mirror_view(discards);

        // set default ids
        default_ids(offsets_t[mpi_rank], comm);

        // valid particles
        auto range = Kokkos::RangePolicy<exec_space>(0, n_reserved);
        particle_masks_initializer<masks_t> pmi{masks, n_valid};
        Kokkos::parallel_for(range, pmi);

        // sync host arrays with device arrays
        checkout_particles();
    }
}

#endif
