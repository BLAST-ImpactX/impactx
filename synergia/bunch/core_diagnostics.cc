#include <cmath>
#include <stdexcept>

#include "core_diagnostics.h"

namespace core_diagnostics_impl {
    struct mean_tag {
        constexpr static int size = 6;
    };
    struct abs_mean_tag {
        constexpr static int size = 6;
    };
    struct z_mean_tag {
        constexpr static int size = 1;
    };
    struct std_tag {
        constexpr static int size = 6;
    };
    struct min_tag {
        constexpr static int size = 3;
    };
    struct max_tag {
        constexpr static int size = 3;
    };
    struct spatial_mean_stddev_tag {
        constexpr static int size = 6;
    };
    struct mom2_tag {
        constexpr static int size = 36;
    };

    template <typename F>
    struct particle_reducer {
        typedef double value_type[];

        const int value_count = F::size;
        const ConstParticles p;
        const ConstParticleMasks masks;
        const karray1d_dev dev_mean;

        particle_reducer(ConstParticles const& p,
                         ConstParticleMasks const& masks,
                         karray1d const& mean = karray1d("mean", 6))
            : p(p), masks(masks), dev_mean("dev_mean", 6)
        {
            Kokkos::deep_copy(dev_mean, mean);
        }

        KOKKOS_INLINE_FUNCTION
        void
        init(value_type dst) const
        {
            for (int j = 0; j < value_count; ++j)
                dst[j] = 0.0;
        }

        KOKKOS_INLINE_FUNCTION
        void
        join(value_type dst, const value_type src) const
        {
            for (int j = 0; j < value_count; ++j)
                dst[j] += src[j];
        }

        KOKKOS_INLINE_FUNCTION
        void operator()(const int i, value_type sum) const;
    };

    // init
    template <>
    KOKKOS_INLINE_FUNCTION void
    particle_reducer<min_tag>::init(value_type dst) const
    {
        for (int j = 0; j < value_count; ++j)
            dst[j] = 1e100;
    }

    template <>
    KOKKOS_INLINE_FUNCTION void
    particle_reducer<max_tag>::init(value_type dst) const
    {
        for (int j = 0; j < value_count; ++j)
            dst[j] = -1e100;
    }

    // join
    template <>
    KOKKOS_INLINE_FUNCTION void
    particle_reducer<min_tag>::join(value_type dst, const value_type src) const
    {
        for (int j = 0; j < value_count; ++j)
            if (dst[j] > src[j]) dst[j] = src[j];
    }

    template <>
    KOKKOS_INLINE_FUNCTION void
    particle_reducer<max_tag>::join(value_type dst, const value_type src) const
    {
        for (int j = 0; j < value_count; ++j)
            if (dst[j] < src[j]) dst[j] = src[j];
    }

    // operator()
    template <>
    KOKKOS_INLINE_FUNCTION void
    particle_reducer<mean_tag>::operator()(const int i, value_type sum) const
    {
        if (masks(i))
            for (int j = 0; j < value_count; ++j)
                sum[j] += p(i, j);
    }

    template <>
    KOKKOS_INLINE_FUNCTION void
    particle_reducer<abs_mean_tag>::operator()(const int i,
                                               value_type sum) const
    {
        if (masks(i))
            for (int j = 0; j < value_count; ++j)
                sum[j] += fabs(p(i, j));
    }

    template <>
    KOKKOS_INLINE_FUNCTION void
    particle_reducer<z_mean_tag>::operator()(const int i, value_type sum) const
    {
        if (masks(i)) sum[0] += p(i, 4);
    }

    template <>
    KOKKOS_INLINE_FUNCTION void
    particle_reducer<std_tag>::operator()(const int i, value_type sum) const
    {
        if (masks(i))
            for (int j = 0; j < value_count; ++j)
                sum[j] += (p(i, j) - dev_mean(j)) * (p(i, j) - dev_mean(j));
    }

    template <>
    KOKKOS_INLINE_FUNCTION void
    particle_reducer<spatial_mean_stddev_tag>::operator()(const int i,
                                                          value_type sum) const

    {
        if (masks(i)) {
            sum[0] += p(i, 0);
            sum[1] += p(i, 2);
            sum[2] += p(i, 4);
            sum[3] += p(i, 0) * p(i, 0);
            sum[4] += p(i, 2) * p(i, 2);
            sum[5] += p(i, 4) * p(i, 4);
        }
    }

    template <>
    KOKKOS_INLINE_FUNCTION void
    particle_reducer<min_tag>::operator()(const int i, value_type min) const
    {
        if (masks(i)) {
            min[0] = (p(i, 0) < min[0]) ? p(i, 0) : min[0];
            min[1] = (p(i, 2) < min[1]) ? p(i, 2) : min[1];
            min[2] = (p(i, 4) < min[2]) ? p(i, 4) : min[2];

            // if (p(i,0) < min[0]) min[0] = p(i,0);
            // if (p(i,2) < min[1]) min[1] = p(i,2);
            // if (p(i,4) < min[2]) min[2] = p(i,4);
        }
    }

    template <>
    KOKKOS_INLINE_FUNCTION void
    particle_reducer<max_tag>::operator()(const int i, value_type max) const
    {
        if (masks(i)) {
            max[0] = (p(i, 0) > max[0]) ? p(i, 0) : max[0];
            max[1] = (p(i, 2) > max[1]) ? p(i, 2) : max[1];
            max[2] = (p(i, 4) > max[2]) ? p(i, 4) : max[2];

            // if (p(i,0) > max[0]) max[0] = p(i,0);
            // if (p(i,2) > max[1]) max[1] = p(i,2);
            // if (p(i,4) > max[2]) max[2] = p(i,4);
        }
    }

    template <>
    KOKKOS_INLINE_FUNCTION void
    particle_reducer<mom2_tag>::operator()(const int i, value_type sum) const
    {
        if (masks(i)) {
            for (int j = 0; j < 6; ++j) {
                double diff_j = p(i, j) - dev_mean(j);
                for (int k = 0; k <= j; ++k) {
                    double diff_k = p(i, k) - dev_mean(k);
                    sum[j * 6 + k] += diff_j * diff_k;
                }
            }
        }
    }
}

karray1d
Core_diagnostics::calculate_mean(Bunch const& bunch)
{
    using core_diagnostics_impl::mean_tag;
    using core_diagnostics_impl::particle_reducer;

    karray1d mean("mean", 6);

    auto bparts = bunch.get_bunch_particles();
    if (bparts.get_memory_location() == MemoryLocation::Host) {
        std::runtime_error("Bunch particles are active on Host memory space!");
    }

    auto particles = bparts.parts;
    auto masks = bparts.masks;
    const int npart = bunch.size();

    particle_reducer<mean_tag> pr(particles, masks);
    Kokkos::parallel_reduce("cal_mean", npart, pr, mean);
    Kokkos::fence();

    MPI_Allreduce(
        MPI_IN_PLACE, mean.data(), 6, MPI_DOUBLE, MPI_SUM, bunch.get_comm());

    for (int i = 0; i < 6; ++i)
        mean(i) /= bunch.get_total_num();

    return mean;
}

double
Core_diagnostics::calculate_z_mean(Bunch const& bunch)
{
    using core_diagnostics_impl::particle_reducer;
    using core_diagnostics_impl::z_mean_tag;

    double mean = 0;

    auto bparts = bunch.get_bunch_particles();
    if (bparts.get_memory_location() == MemoryLocation::Host) {
        std::runtime_error("Bunch particles are active on Host memory space!");
    }
    auto particles = bparts.parts;
    auto masks = bparts.masks;
    const int npart = bunch.size();

    particle_reducer<z_mean_tag> pr(particles, masks);
    Kokkos::parallel_reduce(npart, pr, mean);
    Kokkos::fence();

    MPI_Allreduce(
        MPI_IN_PLACE, &mean, 1, MPI_DOUBLE, MPI_SUM, bunch.get_comm());

    mean = mean / bunch.get_total_num();

    return mean;
}

karray1d
Core_diagnostics::calculate_abs_mean(Bunch const& bunch)
{
    using core_diagnostics_impl::abs_mean_tag;
    using core_diagnostics_impl::particle_reducer;

    karray1d abs_mean("abs_mean", 6);

    auto bparts = bunch.get_bunch_particles();
    if (bparts.get_memory_location() == MemoryLocation::Host) {
        std::runtime_error("Bunch particles are active on Host memory space!");
    }
    auto particles = bparts.parts;
    auto masks = bparts.masks;
    const int npart = bunch.size();

    particle_reducer<abs_mean_tag> pr(particles, masks);
    Kokkos::parallel_reduce("cal_abs_mean", npart, pr, abs_mean);
    Kokkos::fence();

    MPI_Allreduce(MPI_IN_PLACE,
                  abs_mean.data(),
                  6,
                  MPI_DOUBLE,
                  MPI_SUM,
                  bunch.get_comm());

    for (int i = 0; i < 6; ++i)
        abs_mean(i) /= bunch.get_total_num();

    return abs_mean;
}

karray1d
Core_diagnostics::calculate_std(Bunch const& bunch, karray1d const& mean)
{
    using core_diagnostics_impl::particle_reducer;
    using core_diagnostics_impl::std_tag;

    karray1d std("std", 6);

    auto bparts = bunch.get_bunch_particles();
    if (bparts.get_memory_location() == MemoryLocation::Host) {
        std::runtime_error("Bunch particles are active on Host memory space!");
    }
    auto particles = bparts.parts;
    auto masks = bparts.masks;
    const int npart = bunch.size();

    particle_reducer<std_tag> pr(particles, masks, mean);
    Kokkos::parallel_reduce(npart, pr, std);
    Kokkos::fence();

    MPI_Allreduce(
        MPI_IN_PLACE, std.data(), 6, MPI_DOUBLE, MPI_SUM, bunch.get_comm());
    for (int i = 0; i < 6; ++i)
        std(i) = std::sqrt(std(i) / bunch.get_total_num());

    return std;
}

karray2d_row
Core_diagnostics::calculate_sum2(Bunch const& bunch, karray1d const& mean)
{
    using core_diagnostics_impl::mom2_tag;
    using core_diagnostics_impl::particle_reducer;

    karray2d_row sum2("sum2", 6, 6);

    auto bparts = bunch.get_bunch_particles();
    if (bparts.get_memory_location() == MemoryLocation::Host) {
        std::runtime_error("Bunch particles are active on Host memory space!");
    }
    auto particles = bparts.parts;
    auto masks = bparts.masks;
    auto npart = bunch.size();

    particle_reducer<mom2_tag> pr(particles, masks, mean);
    Kokkos::parallel_reduce(npart, pr, sum2);
    Kokkos::fence();

    for (int i = 0; i < 5; ++i)
        for (int j = i + 1; j < 6; ++j)
            sum2(i, j) = sum2(j, i);

    MPI_Allreduce(
        MPI_IN_PLACE, sum2.data(), 36, MPI_DOUBLE, MPI_SUM, bunch.get_comm());

    return sum2;
}

karray2d_row
Core_diagnostics::calculate_mom2(Bunch const& bunch, karray1d const& mean)
{
    auto sum2 = calculate_sum2(bunch, mean);
    karray2d_row mom2("mom2", 6, 6);

    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            mom2(i, j) = sum2(i, j) / bunch.get_total_num();

    return mom2;
}

karray1d
Core_diagnostics::calculate_min(Bunch const& bunch)
{
    using core_diagnostics_impl::min_tag;
    using core_diagnostics_impl::particle_reducer;

    karray1d min("min", 3);
    min(0) = 1e100;
    min(1) = 1e100;
    min(2) = 1e100;

    auto bparts = bunch.get_bunch_particles();
    if (bparts.get_memory_location() == MemoryLocation::Host) {
        std::runtime_error("Bunch particles are active on Host memory space!");
    }
    auto particles = bparts.parts;
    auto masks = bparts.masks;
    const int npart = bunch.size();

    particle_reducer<min_tag> pr(particles, masks);
    Kokkos::parallel_reduce("cal_min", npart, pr, min);
    Kokkos::fence();

    MPI_Allreduce(
        MPI_IN_PLACE, min.data(), 3, MPI_DOUBLE, MPI_MIN, bunch.get_comm());

    return min;
}

karray1d
Core_diagnostics::calculate_max(Bunch const& bunch)
{
    using core_diagnostics_impl::max_tag;
    using core_diagnostics_impl::particle_reducer;

    karray1d max("max", 3);
    max(0) = -1e100;
    max(1) = -1e100;
    max(2) = -1e100;

    auto bparts = bunch.get_bunch_particles();
    if (bparts.get_memory_location() == MemoryLocation::Host) {
        std::runtime_error("Bunch particles are active on Host memory space!");
    }
    auto particles = bparts.parts;
    auto masks = bparts.masks;
    const int npart = bunch.size();

    particle_reducer<max_tag> pr(particles, masks);
    Kokkos::parallel_reduce("cal_max", npart, pr, max);
    Kokkos::fence();

    MPI_Allreduce(
        MPI_IN_PLACE, max.data(), 3, MPI_DOUBLE, MPI_MAX, bunch.get_comm());

    return max;
}

karray1d
Core_diagnostics::calculate_spatial_mean_stddev(Bunch const& bunch)
{
    using core_diagnostics_impl::particle_reducer;
    using core_diagnostics_impl::spatial_mean_stddev_tag;

    auto bparts = bunch.get_bunch_particles();
    if (bparts.get_memory_location() == MemoryLocation::Host) {
        std::runtime_error("Bunch particles are active on Host memory space!");
    }
    auto particles = bparts.parts;
    auto masks = bparts.masks;
    const int npart = bunch.size();

    const auto total_bunch_particles = bunch.get_total_num();
    const auto local_bunch_capacity = bunch.size();

    karray1d mean_and_stddev("mean_stddev", 6);

    particle_reducer<spatial_mean_stddev_tag> pr(particles, masks);
    Kokkos::parallel_reduce("spatial_mean_stddev", npart, pr, mean_and_stddev);
    Kokkos::fence();

    if (MPI_Allreduce(MPI_IN_PLACE,
                      mean_and_stddev.data(),
                      6,
                      MPI_DOUBLE,
                      MPI_SUM,
                      bunch.get_comm()) != MPI_SUCCESS) {
        std::runtime_error("MPI_Allreduce error");
    }

    /* mean_and_stddev has sum_and_sumsquares on each MPI rank */
    for (int j = 0; j < 3; j++) {
        mean_and_stddev(j) = mean_and_stddev(j) / total_bunch_particles;
        mean_and_stddev(j + 3) =
            std::sqrt((mean_and_stddev(j + 3)) / total_bunch_particles -
                      std::pow(mean_and_stddev(j), 2));
    }

    return mean_and_stddev;
}

double
calculate_median_cdt(Bunch const& bunch)
{
    throw std::runtime_error("bunch cdt median not yet implemented");
}

std::vector<double>
Core_diagnostics::kokkos_view_to_stl_vector(karray1d const& view)
{
    if (!view.span_is_contiguous()) {
        std::runtime_error("non-contigipus view provided!");
    }
    size_t view_size = view.size();
    std::vector<double> retval;
    retval.resize(view_size);
    size_t data_size = sizeof(double) * view_size;
    std::memcpy(retval.data(), view.data(), data_size);
    return retval;
}

std::vector<double>
Core_diagnostics::kokkos_view_to_stl_vector(karray2d_row const& view)
{
    if (!view.span_is_contiguous()) {
        std::runtime_error("non-contigipus view provided!");
    }
    size_t view_size = view.size();
    std::vector<double> retval;
    retval.resize(view_size);
    size_t data_size = sizeof(double) * view_size;
    std::memcpy(retval.data(), view.data(), data_size);
    return retval;
}
