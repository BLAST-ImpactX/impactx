#include <sstream>
#include <stdexcept>

#include "core_diagnostics.h"
#include "populate.h"
#include "populate_host.h"
#include "synergia/utils/floating_point.h"

constexpr auto pi = Kokkos::numbers::pi_v<double>;

namespace {
    bool
    is_symmetric66(const_karray2d_row m)
    {
        bool symmetric = true;
        const double tolerance = 1.0e-14;

        for (int i = 0; i < 6; ++i) {
            for (int j = i + 1; j < 6; ++j) {
                if (!floating_point_equal(m(i, j), m(j, i), tolerance)) {
                    symmetric = false;
                }
            }
        }

        return symmetric;
    }
}

void
adjust_moments(Bunch& bunch,
               const_karray1d means,
               const_karray2d_row covariances)
{
    if (!is_symmetric66(covariances))
        throw std::runtime_error(
            "adjust_moments: covariance matrix must be symmetric");

    // calculate_mean and mom2 are performed on device memory, so we need to
    // copy the particle data from host to device first. checkout is not
    // necessary since the core diagnostics do not change the particle data
    bunch.checkin_particles();

    karray1d bunch_mean = Core_diagnostics::calculate_mean(bunch);
    karray2d_row bunch_mom2 =
        Core_diagnostics::calculate_mom2(bunch, bunch_mean);

    int num_particles = bunch.size();
    int num_particles_slots = bunch.capacity();

    // auto strides = bunch.get_particle_strides();
    // int num_particles_slots = strides[1];

    adjust_moments_host(means.data(),
                        covariances.data(),
                        bunch_mean.data(),
                        bunch_mom2.data(),
                        num_particles,
                        num_particles_slots,
                        bunch.get_host_particles().data());
}

namespace {
    void
    fill_unit_6d(Distribution& dist,
                 HostParticles particles,
                 const_karray2d_row covariances,
                 int start,
                 int end)
    {
        for (int j = 0; j < 6; ++j) {
            const double scale = sqrt(covariances(j, j));

            for (int p = start; p < end; ++p)
                particles(p, j) = dist.get_unit_gaussian() * scale;
        }
    }

    inline bool
    good(ConstHostParticles particles, const_karray1d limits, int index)
    {
        bool retval = true;

        for (int i = 0; i < 6; ++i) {
            double val = particles(index, i);
            double limit = limits[i];

            if ((limit > 0) && ((val > limit) or (val < -limit)))
                retval = false;
        }

        return retval;
    }

    void
    strip_unit_6d(Bunch& bunch,
                  const_karray1d limits,
                  int& total_num,
                  int& local_num)
    {
        auto particles = bunch.get_host_particles();
        local_num = bunch.get_local_num();

        int index = 0;
        while (index < local_num) {
            if (good(particles, limits, index)) {
                ++index;
            } else {
                int last = local_num - 1;
                if (good(particles, limits, last)) {
                    particles(index, 0) = particles(last, 0);
                    particles(index, 1) = particles(last, 1);
                    particles(index, 2) = particles(last, 2);
                    particles(index, 3) = particles(last, 3);
                    particles(index, 4) = particles(last, 4);
                    particles(index, 5) = particles(last, 5);
                }

                --local_num;
            }
        }

        MPI_Allreduce(
            &local_num, &total_num, 1, MPI_INT, MPI_SUM, bunch.get_comm());
    }
}

void
populate_6d(Distribution& dist,
            Bunch& bunch,
            const_karray1d means,
            const_karray2d_row covariances)
{
    if (bunch.size() != bunch.get_local_num()) {
        throw std::runtime_error(
            "populate_6d: "
            "cannot populate bunches that has already lost particles.");
    }

    karray1d limits("limits", 6);
    for (int i = 0; i < 6; ++i)
        limits[i] = 0.0;

    populate_6d_truncated(dist, bunch, means, covariances, limits);
}

void
populate_6d_truncated(Distribution& dist,
                      Bunch& bunch,
                      const_karray1d means,
                      const_karray2d_row covariances,
                      const_karray1d limits)
{
#if 0
    multi_array_assert_size(means, 6, "populate_6d: means");
    multi_array_assert_size(covariances, 6, 6, "populate_6d: covariances");
    multi_array_assert_size(limits, 6, "populate_6d: limits");
#endif

    // deep copy from device to host
    bunch.checkout_particles();

    auto particles = bunch.get_host_particles();

    karray2d_row unit_covariances("unit_covariances", 6, 6);
    karray1d zero_means("zero_means", 6);

    bool truncated(false);

    for (int i = 0; i < 6; ++i) {
        double n = limits[i];

        if (n > 0) {
            truncated = true;

            double cutoff_integral =
                (exp(-n * n / 2.0)) *
                (sqrt(pi) * exp(n * n / 2.0) * erf(n / sqrt(2.0)) -
                 sqrt(2.0) * n) /
                (sqrt(pi));

            unit_covariances(i, i) = 1.0 / (cutoff_integral * cutoff_integral);
        } else {
            unit_covariances(i, i) = 1.0;
        }
    }

    int start = 0;
    int end = bunch.get_local_num();

    fill_unit_6d(dist, particles, unit_covariances, start, end);

    if (truncated) {
        adjust_moments(bunch, zero_means, unit_covariances);

        int iteration = 0;
        int total_num, local_num;

        strip_unit_6d(bunch, limits, total_num, local_num);

        while (total_num < bunch.get_total_num()) {
            ++iteration;
            const int max_iterations = 50;
            if (iteration > max_iterations) {
                throw std::runtime_error(
                    "populate_6d_truncated: "
                    "maximum number of truncation iterations exceeded. "
                    "Algorithm known to fail ~< 2.5 sigma.");
            }

            fill_unit_6d(dist, particles, unit_covariances, local_num, end);
            adjust_moments(bunch, zero_means, unit_covariances);
            strip_unit_6d(bunch, limits, total_num, local_num);
        }
    }

    adjust_moments(bunch, means, covariances);

    // copy to device
    bunch.checkin_particles();

    // check
    bunch.check_pz2_positive();
}

karray2d_row
get_correlation_matrix(const_karray2d_row one_turn_map,
                       double arms,
                       double brms,
                       double crms,
                       double beta,
                       std::array<int, 3> const& rms_index)
{

    for (int idx : rms_index) {
        if (idx < 0 || idx > 5)
            throw std::runtime_error("only valid rms indices (x=0, xp=1, y=2, "
                                     "yp=3, z=4, dpp=5) are allowed");
    }

    karray2d_row correlation_matrix("corr_matrix", 6, 6);

    get_correlation_matrix_host(correlation_matrix.data(),
                                one_turn_map.data(),
                                arms,
                                brms,
                                crms,
                                beta,
                                rms_index);

    return correlation_matrix;
}

void
populate_transverse_gaussian(Distribution& dist,
                             Bunch& bunch,
                             const_karray1d means,
                             const_karray2d_row covariances,
                             double cdt)
{
    // deep copy from device to host
    bunch.checkout_particles();
    auto particles = bunch.get_host_particles();

    for (int i = 0; i < 4; ++i) {
        for (int p = 0; p < bunch.size(); ++p)
            particles(p, i) = dist.get_unit_gaussian();
    }

    // cdt
    for (int p = 0; p < bunch.size(); ++p)
        particles(p, 4) = dist.get_uniform(0.0, 1.0);

    // dpop
    for (int p = 0; p < bunch.size(); ++p)
        particles(p, 5) = dist.get_unit_gaussian();

    // copy of original means and covariances
    karray1d means_modified("means_modified", means.layout());
    karray2d_row covariances_modified("covariances_modified",
                                      covariances.layout());

    Kokkos::deep_copy(means_modified, means);
    Kokkos::deep_copy(covariances_modified, covariances);

    means_modified[4] = 0.0;

    // Symmetry requires no correlations with the cdt coordinate. Make a copy
    // of the covariance matrix and manually set all correlations to zero.
    for (int k = 0; k < 6; ++k)
        covariances_modified(k, 4) = covariances_modified(4, k) = 0.0;

    covariances_modified(4, 4) = cdt * cdt / 12.0;
    adjust_moments(bunch, means_modified, covariances_modified);

    // copy to device
    bunch.checkin_particles();
}
