

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Kokkos_Core.hpp>
#include <Kokkos_NumericTraits.hpp>
#include <Kokkos_Random.hpp>

#include "synergia/bunch/bunch.h"
#include "synergia/bunch/bunch_particles.h"
#include "synergia/foundation/physical_constants.h"

constexpr double mass = pconstants::mp;
constexpr double KE0 = 0.8;
constexpr double total_energy = mass + KE0;
constexpr int total_num = 16;
constexpr double real_num = 2.0e12;
constexpr auto x = Bunch::x;
constexpr auto xp = Bunch::xp;
constexpr auto y = Bunch::y;
constexpr auto yp = Bunch::yp;
constexpr auto cdt = Bunch::cdt;
constexpr auto dpop = Bunch::dpop;
constexpr auto id = Bunch::id;
constexpr int num_parts = 16;

TEST_CASE("Bunch", "[Bunch]")
{
    Four_momentum fm(mass, total_energy);
    Reference_particle ref(pconstants::proton_charge, fm);
    Bunch bunch(ref, num_parts, 1e13, Commxx());

    bunch.checkout_particles();
    auto parts = bunch.get_host_particles();

    auto localnum = bunch.get_local_num();
    for(auto i=0; i<localnum; ++i) {
        for(auto j=0; j<6; ++j) {
            parts(i, j) = 0.0;
        }
        parts(i, Bunch::cdt) = i * 0.05;
        std::cout << "before offset parts(i, 4): " << parts(i, 4) << std::endl;
    }

    bunch.checkin_particles();

    // Offset the bunch cdt
    bunch.offset_bunch_particles_cdt(-0.02);

    // particles changed on the device, so check them out for testing    
    bunch.checkout_particles();

    auto new_parts = bunch.get_host_particles();

    for (auto i=0; i<localnum; ++i) {
        std::cout << "after offset new_parts(i, 4): " << new_parts(i, 4) << std::endl;
        // check that cdt has been offset
        CHECK_THAT(parts(i, Bunch::cdt), Catch::Matchers::WithinAbs(i*0.05-0.02, 1.0e-12));
    }
}
