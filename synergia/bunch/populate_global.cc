#include <sstream>
#include <stdexcept>

#include "synergia/bunch/core_diagnostics.h"
#include "synergia/bunch/populate_global.h"
#include "synergia/utils/floating_point.h"

#include "synergia/utils/pcg/pcg_random.hpp"
#include <random>

constexpr auto pi = Kokkos::numbers::pi_v<double>;

namespace {
    void
    fill_unit_6d(std::vector<pcg64>& rngs,
                 HostParticles parts,
                 HostParticleMasks masks,
                 HostParticleMasks valid,
                 const_karray2d_row covariances)
    {
        const int np = valid.extent(0);
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        for (int j = 0; j < 6; ++j) {
            const double scale = sqrt(covariances(j, j));

            for (int p = 0; p < np; ++p) {
                if (masks(p) && !valid(p)) parts(p, j) = dist(rngs[p]) * scale;
            }
        }
    }

    void
    fill_unit_6d(uint64_t seed,
                 HostParticles parts,
                 HostParticleMasks masks,
                 int np,
                 const_karray2d_row covariances)
    {
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        std::array<double, 6> scales;
        for (int j = 0; j < 6; ++j)
            scales[j] = sqrt(covariances(j, j));

        for (int p = 0; p < np; ++p) {
            if (!masks(p)) continue;

            pcg64 rng(seed, (uint64_t)parts(p, 6));

            for (int j = 0; j < 6; ++j)
                parts(p, j) = dist(rng) * scales[j];
        }
    }

    inline bool
    good(ConstHostParticles parts, int index, const_karray1d limits)
    {
        bool retval = true;

        for (int i = 0; i < 6; ++i) {
            double val = parts(index, i);
            double limit = limits[i];

            if ((limit > 0) && ((val > limit) or (val < -limit)))
                retval = false;
        }

        return retval;
    }

    int
    strip_unit_6d(ConstHostParticles parts,
                  HostParticleMasks masks,
                  HostParticleMasks valid,
                  const_karray1d limits)
    {
        const int np = valid.extent(0);
        int failed_num = 0;

        for (int p = 0; p < np; ++p) {
            if (masks(p) && !valid(p)) {
                if (good(parts, p, limits))
                    valid(p) = 1;
                else
                    ++failed_num;
            }
        }

        return failed_num;
    }
}

void
populate_global_6d(uint64_t seed,
                   Bunch& bunch,
                   const_karray1d means,
                   const_karray2d_row covariances)
{
    karray1d limits("limits", 6);
    for (int i = 0; i < 6; ++i)
        limits[i] = 0.0;

    populate_global_6d_truncated(seed, bunch, means, covariances, limits);
}

void
populate_global_6d_truncated(uint64_t seed,
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
    // not for particles, but for masks
    bunch.checkout_particles();

    const int np = bunch.size();
    auto parts = bunch.get_host_particles();
    auto masks = bunch.get_host_particle_masks();

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

    if (truncated) {
        // a mask indicating whether the generated particle is valid
        // defaults to 0 for all particles
        HostParticleMasks valid("valid", np);

        // random number generators, one for each particle
        std::vector<pcg64> rngs;
        rngs.reserve(np);

        for (int i = 0; i < np; ++i)
            rngs[i] = pcg64(seed, (uint64_t)parts(i, 6));

        // loop to generate all particles
        // total_bad is init to 1 because the MPI_Allreduce has to be
        // called by all ranks in the bunch. If the total_bad is init
        // to np, then its possible that some of the ranks who does not
        // have local particles will not run the allreduce.
        const int max_iters = 50;
        int iter = 0;
        int total_bad = 1;

        while (total_bad && iter < max_iters) {
            fill_unit_6d(rngs, parts, masks, valid, unit_covariances);
            adjust_moments(bunch, zero_means, unit_covariances);
            int bad = strip_unit_6d(parts, masks, valid, limits);

            MPI_Allreduce(
                &bad, &total_bad, 1, MPI_INT, MPI_SUM, bunch.get_comm());

            ++iter;
        }

        if (total_bad) {
            throw std::runtime_error(
                "populate_6d_truncated: "
                "maximum number of truncation iterations exceeded. "
                "Algorithm known to fail ~< 2.5 sigma.");
        }
    } else {
        fill_unit_6d(seed, parts, masks, np, unit_covariances);
    }

    // adjust
    adjust_moments(bunch, means, covariances);

    // copy to device
    bunch.checkin_particles();

    // check
    bunch.check_pz2_positive();
}
