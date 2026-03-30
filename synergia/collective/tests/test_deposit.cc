#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "catch2/matchers/catch_matchers.hpp"
#include "synergia/foundation/physical_constants.h"

#include "synergia/collective/deposit.h"

const double mass = 100.0;
const double total_energy = 125.0;

const int total_num = 100;
const double real_num = 2.0e12;
auto const id = Bunch::id;

TEST_CASE("NoParticles", "[NoParticles]")
{
    Four_momentum fm(mass, total_energy);
    Reference_particle ref(pconstants::proton_charge, fm);

    Bunch bunch(ref, 0, 0, Commxx());
    CHECK(bunch.get_local_num() == 0);

    Rectangular_grid_domain domain(
        {4, 4, 4}, {1e-6, 1e-6, 1e-6}, {0, 0, 0}, false);

    const std::array<int, 3> dims{8, 8, 8};
    karray1d_dev rho_dev("rho_dev", 512);

    deposit_charge_rectangular_3d_kokkos_scatter_view(
        rho_dev, domain, dims, bunch);

    karray1d_hst rho_dev_hst = Kokkos::create_mirror_view(rho_dev);
    Kokkos::deep_copy(rho_dev_hst, rho_dev);
    Kokkos::fence();

    for (int idx = 0; idx < 512; idx++) {
        REQUIRE_THAT(rho_dev_hst(0), Catch::Matchers::WithinAbs(0, 0.01));
    }
}

TEST_CASE("OneParticle", "[OneParticle]")
{
    Four_momentum fm(mass, total_energy);
    Reference_particle ref(pconstants::proton_charge, fm);

    Bunch bunch(ref, 1, 1, Commxx());
    CHECK(bunch.get_local_num() == 1);

    bunch.checkout_particles();
    auto bunch_parts = bunch.get_host_particles();
    bunch_parts.access(0, 0) = 1.55;
    bunch_parts.access(0, 2) = 1.55;
    bunch_parts.access(0, 4) = 1.55;
    bunch.checkin_particles();

    Rectangular_grid_domain domain({6, 6, 6}, {6, 6, 6}, {2, 2, 2}, false);

    auto h = domain.get_cell_size();

    double weight0 = (bunch.get_real_num() / bunch.get_total_num()) *
                     bunch.get_particle_charge() * pconstants::e /
                     (h[0] * h[1] * h[2]);

    const std::array<int, 3> dims{6, 6, 6};
    karray1d_dev rho_dev("rho_dev", dims[0] * dims[1] * dims[2]);

    deposit_charge_rectangular_3d_kokkos_scatter_view(
        rho_dev, domain, dims, bunch);

    karray1d_hst rho_dev_hst = Kokkos::create_mirror_view(rho_dev);
    Kokkos::deep_copy(rho_dev_hst, rho_dev);
    Kokkos::fence();

    // where we expect the particle to be deposited
    std::array<int, 2> deposit_locs{1, 2};

    for (auto x : deposit_locs) {
        auto fact_x = (x == 0) ? 0.95 : 0.0;
        for (auto y : deposit_locs) {
            auto fact_y = (y == 0) ? 0.95 : 0.0;
            for (auto z : deposit_locs) {
                auto fact_z = (z == 0) ? 0.95 : 0.0;
                REQUIRE_THAT(
                    rho_dev_hst(z * dims[0] * dims[1] + y * dims[0] + z),
                    Catch::Matchers::WithinAbs(fact_x * fact_y * fact_z, 0.01));
            }
        }
    }

    auto sums = 0.0;
    for (int idx = 0; idx < dims[0] * dims[1] * dims[2]; idx++) {
        sums += rho_dev_hst(idx) / weight0;
    }

    // one particle is deposited
    REQUIRE_THAT(sums, Catch::Matchers::WithinAbs(1, 0.01));
}
