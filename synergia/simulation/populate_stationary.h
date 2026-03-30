#ifndef POPULATE_STATIONARY_H_
#define POPULATE_STATIONARY_H_

#include "synergia/foundation/distribution.h"
#include "synergia/simulation/lattice_simulator.h"

/// Populate a bunch with a shell of particles having fixed constant actions
/// but uniform in phase angle in all three planes in normal form space.
/// Any distribution or statistic constructed from phase space variables
/// will be stationary under propagation through the lattice (up to
/// limited statistics) if there is no other physics.
/// @param dist the distribution generator
/// @param bunch the bunch
/// @param actions std::vector<double> (3) the three mean actions

std::array<std::complex<double>, 3>
get_6d_normal_form_coords(Distribution& dist,
                          std::array<double, 3> const& actions)
{
    std::array<std::complex<double>, 3> part;

    for (int c = 0; c < 3; ++c) {
        double phase =
            dist.get_uniform(0.0, 2.0 * Kokkos::numbers::pi_v<double>);
        double square_root_action =
            sqrt(-actions[c] * log(1.0 - dist.get_uniform(0.0, 1.0)));

        part[c] = std::complex<double>(square_root_action * sin(phase),
                                       -square_root_action * cos(phase));
    }

    return part;
}

template <unsigned int order>
void
populate_6d_stationary_clipped_longitudinal_gaussian(
    Distribution& dist,
    Bunch& bunch,
    std::array<double, 3> const& actions,
    double cdt_min,
    double cdt_max,
    NormalForm<order> const& nf)
{
#ifdef SYNERGIA_ENABLE_CUDA

    // no implementation

#else

    const int max_tries = 100;

    auto parts = bunch.get_host_particles();
    auto np = bunch.size();

    for (int p = 0; p < np; ++p) {
        std::array<double, 6> test_p;
        int curr_try = 0;

        while (curr_try < max_tries) {
            auto nf_p = get_6d_normal_form_coords(dist, actions);
            bool good_particle = true;

            for (int phase = 0; phase < 4; ++phase) {
                test_p = nf.cnvDataFromNormalForm(nf_p);

                if (test_p[4] < cdt_min || test_p[4] > cdt_max) {
                    good_particle = false;
                    break;
                }

                double a2r = nf_p[2].real();
                double a2i = nf_p[2].imag();
                nf_p[2] = std::complex<double>(a2i, -a2r);
            }

            // found a good one
            if (good_particle) break;

            // have another try
            ++curr_try;
        }

        if (curr_try == max_tries)
            throw std::runtime_error(
                "populate stationary: couldnt produce good particle");

        for (int i = 0; i < 6; ++i)
            parts(p, i) = test_p[i];
    }

    bunch.checkin_particles();

#endif // SYNERGIA_ENABLE_CUDA
}

#endif /* POPULATE_STATIONARY_H_ */
