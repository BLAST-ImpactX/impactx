#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "catch2/matchers/catch_matchers.hpp"
#include "synergia/collective/space_charge_3d_fd.h"
#include "synergia/foundation/physical_constants.h"

TEST_CASE("PointCharge", "[PointCharge]")
{
    PetscErrorCode ierr;

    auto logger = Logger(0, LoggerV::DEBUG);
    auto simlogger = Logger(0, LoggerV::INFO_STEP);

    const double mass = 100.0;
    const double total_energy = 125.0;
    const int total_num = 100;
    const double real_num = 2.0e12;

    /* scope to prevent deallocate after finalize errors for Kokkos views! */
    {
        Four_momentum fm(mass, total_energy);
        Reference_particle ref(pconstants::proton_charge, fm);

        auto bsim =
            Bunch_simulator::create_single_bunch_simulator(ref, 1, 1, Commxx());

        Bunch& bunch = bsim.get_bunch();
        {
            bunch.checkout_particles();
            auto bunch_parts = bunch.get_host_particles();
            bunch_parts.access(0, 0) = 0;
            bunch_parts.access(0, 2) = 0;
            bunch_parts.access(0, 4) = 0;
            // print intital coordinates
            logger << "before kick, particle at (x y z):" << '\n';
            for (int k = 0; k < 1; ++k) {
                logger << k << ": " << bunch_parts(k, 0) << ", "
                       << bunch_parts(k, 2) << ", " << bunch_parts(k, 4)
                       << '\n';
            }
            bunch.checkin_particles();
        }

        // space charge operator options
        auto sc_ops = Space_charge_3d_fd_options(17, 17, 17);
        sc_ops.comm_group_size = 1;

        // set domain
        std::array<double, 3> offset = {0, 0, 0};
        std::array<double, 3> size = {4, 4, 4};
        sc_ops.set_fixed_domain(offset, size);

        // space charge operator
        auto sc = Space_charge_3d_fd(sc_ops);

        // apply space charge operator
        sc.apply(bsim, 1e-6, simlogger);

        {
            bunch.checkout_particles();
            auto bunch_parts = bunch.get_host_particles();
            // print final coordinates
            logger << "after kick, particle at (x y z):" << '\n';
            for (int k = 0; k < 1; ++k) {
                logger << k << ": " << bunch_parts(k, 0) << ", "
                       << bunch_parts(k, 2) << ", " << bunch_parts(k, 4)
                       << '\n';
            }
            logger << "after kick, particle momenta:" << '\n';
            for (int k = 0; k < 1; ++k) {
                logger << k << ": " << bunch_parts(k, 1) << ", "
                       << bunch_parts(k, 3) << ", " << bunch_parts(k, 5)
                       << '\n';
            }

            REQUIRE_THAT(bunch_parts(0, 1), Catch::Matchers::WithinAbs(0, .01));
            REQUIRE_THAT(bunch_parts(0, 3), Catch::Matchers::WithinAbs(0, .01));
            REQUIRE_THAT(bunch_parts(0, 5), Catch::Matchers::WithinAbs(0, .01));
        }
    }
}
