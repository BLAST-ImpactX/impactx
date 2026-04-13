

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
    }
    // also check transverse momenta
    for (auto i=0; i<localnum; ++i) {
        parts(i, Bunch::xp) = i * 1.0e-4;
        parts(i, Bunch::yp) = (localnum-i) * 1.0e-4;
    }

    constexpr double dE = 0.050; // 50 MeV
    double energy = ref.get_total_energy();
    double old_p = sqrt(energy*energy - mass*mass);
	

    // Add some energy to each particle by changing dp/p

    // what is that in dp/p
    double new_energy = energy + dE;
    double new_p = sqrt(new_energy*new_energy - mass*mass);
    double new_dpop = new_p/old_p - 1.0;

    for (auto i=0; i<localnum; ++i) {
        parts(i, Bunch::dpop) = new_dpop;
    }

    bunch.checkin_particles();

    // adjust the reference energy
    bunch.adjust_bunch_particles_reference_energy(new_energy);

    bunch.checkout_particles();
    auto new_parts = bunch.get_host_particles();

    for (auto i=0; i<localnum; ++i) {
        std::cout << "new_parts(i, 5): " << new_parts(i, 5) << std::endl;
        // check that dp/p is now close to 0
        CHECK_THAT(parts(i, Bunch::dpop), Catch::Matchers::WithinAbs(0.0, 1.0e-12));
        // transverse momenta
		std::cout << "new_parts(" << i << "(, 1): " << new_parts(i, 1) << std::endl;
        CHECK_THAT(parts(i, Bunch::xp), Catch::Matchers::WithinRel(i*1.0e-4*old_p/new_p,  1.0e-12));
		std::cout << "new_parts(" << i << "(, 3): " << new_parts(i, 3) << std::endl;
        CHECK_THAT(parts(i, Bunch::yp), Catch::Matchers::WithinRel((localnum-i)*1.0e-4*old_p/new_p,  1.0e-12));
    }
}
