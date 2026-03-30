#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "catch2/matchers/catch_matchers.hpp"
#include "synergia/foundation/physical_constants.h"
#include "synergia/lattice/madx_reader.h"

#include "synergia/simulation/independent_stepper_elements.h"
#include "synergia/simulation/propagator.h"

struct propagator_fixture {
    Logger screen;
    Lattice lattice;
    std::unique_ptr<Propagator> propagator;
    std::unique_ptr<Bunch_simulator> sim;

    propagator_fixture(std::string const& file)
        : screen(0, LoggerV::INFO_TURN)
        , lattice(MadX_reader().get_lattice("machine",
                                            "./lattices/" + file + ".seq"))
        , propagator()
        , sim()
    {
        // Proton 1.5GeV/c
        auto fm = Four_momentum(pconstants::mp);
        fm.set_momentum(1.5);
        auto ref = Reference_particle(1, fm);

        // set refpart to lattice
        lattice.set_reference_particle(ref);

        // propagator
        propagator = std::make_unique<Propagator>(
            lattice, Independent_stepper_elements(1));

        // bunch simulator
        sim = std::make_unique<Bunch_simulator>(
            Bunch_simulator::create_single_bunch_simulator(
                ref, 16, 1e10, Commxx()));
    }

    void
    propagate()
    {
        propagator->propagate(*sim, screen, 1);
    }

    Bunch&
    bunch()
    {
        return sim->get_bunch();
    }

    void
    print_lattice()
    {
        Logger l(0, LoggerV::DEBUG);
        lattice.print(l);
    }
};

karray2d_row
read_madx_output(std::string const& seq)
{
    karray2d_row out(seq, 16, 6);

    std::ifstream file("./lattices/madx_" + seq + ".out");
    std::string line;

    int row = 0;
    while (std::getline(file, line)) {
        if (row >= 26) {
            std::stringstream ss(line);

            double val;
            int col = 0;

            while (ss >> val) {
                if (col >= 2 && col <= 7) out(row - 26, col - 2) = val;
                ++col;
            }
        }
        ++row;
    }

    return out;
}

void
propagate_libff(std::string const& seq, double tolerance = 1e-15)
{
    std::cout << "libff propagate " << seq << "\n";

    propagator_fixture pf(seq);
    pf.print_lattice();

    // init particles
    auto& b = pf.bunch();

    double s2o2 = sqrt(2.0) / 2.0;
    double s3o2 = sqrt(3.0) / 2.0;
    double offset = 1.0e-3;

    b.checkout_particles();
    auto parts = b.get_host_particles();

    parts(0, 0) = offset;
    parts(0, 2) = 0.0;

    parts(1, 0) = offset * s3o2;
    parts(1, 2) = offset * 0.5;

    parts(2, 0) = offset * s2o2;
    parts(2, 2) = offset * s2o2;

    parts(3, 0) = offset * 0.5;
    parts(3, 2) = offset * s3o2;

    parts(4, 0) = 0.0;
    parts(4, 2) = offset;

    parts(5, 0) = -offset * 0.5;
    parts(5, 2) = offset * s3o2;

    parts(6, 0) = -offset * s2o2;
    parts(6, 2) = offset * s2o2;

    parts(7, 0) = -offset * s3o2;
    parts(7, 2) = offset * 0.5;

    parts(8, 0) = -offset;
    parts(8, 2) = 0.0;

    parts(9, 0) = -offset * s3o2;
    parts(9, 2) = -offset * 0.5;

    parts(10, 0) = -offset * s2o2;
    parts(10, 2) = -offset * s2o2;

    parts(11, 0) = -offset * 0.5;
    parts(11, 2) = -offset * s3o2;

    parts(12, 0) = 0.0;
    parts(12, 2) = -offset;

    parts(13, 0) = offset * 0.5;
    parts(13, 2) = -offset * s3o2;

    parts(14, 0) = offset * s2o2;
    parts(14, 2) = -offset * s2o2;

    parts(15, 0) = offset * s3o2;
    parts(15, 2) = -offset * 0.5;

    b.checkin_particles();

    // propagate
    pf.propagate();

    // print
    b.checkout_particles();
    parts = b.get_host_particles();

    // madx output
    auto madx = read_madx_output(seq);

    // check
    for (int p = 0; p < 16; ++p) {
#if 0
        std::cout << "libff p=" << p << ": ";
        for(int i=0; i<5; ++i) std::cout << parts(p, i) << ", ";
        std::cout << "\n";

        std::cout << "madx  p=" << p << ": ";
        for(int i=0; i<5; ++i) std::cout << madx(p, i) << ", ";
        std::cout << "\n";
#endif

        REQUIRE_THAT(parts(p, 0),
                     Catch::Matchers::WithinAbs(madx(p, 0), tolerance));
        REQUIRE_THAT(parts(p, 1),
                     Catch::Matchers::WithinAbs(madx(p, 1), tolerance));
        REQUIRE_THAT(parts(p, 2),
                     Catch::Matchers::WithinAbs(madx(p, 2), tolerance));
        REQUIRE_THAT(parts(p, 3),
                     Catch::Matchers::WithinAbs(madx(p, 3), tolerance));
        REQUIRE_THAT(parts(p, 4),
                     Catch::Matchers::WithinAbs(-madx(p, 4), tolerance));
        REQUIRE_THAT(parts(p, 5),
                     Catch::Matchers::WithinAbs(madx(p, 5), tolerance));
    }
}

TEST_CASE("cfsbend_dipole", "[libFF][Elements]")
{
    propagate_libff("cfsbend_dipole", 1e-7);
}

TEST_CASE("cfsbend_dipole_skew", "[libFF][Elements]")
{
    propagate_libff("cfsbend_dipole_skew", 1e-7);
}

TEST_CASE("cfsbend_sext", "[libFF][Elements]")
{
    propagate_libff("cfsbend_sext", 1e-7);
}

TEST_CASE("cfsbend_sext_skew", "[libFF][Elements]")
{
    propagate_libff("cfsbend_sext_skew", 1e-7);
}

TEST_CASE("mpole_k1", "[libFF][Elements]")
{
    propagate_libff("mpole_k1");
}

TEST_CASE("mpole_k1_tilt", "[libFF][Elements]")
{
    propagate_libff("mpole_k1_tilt");
}

TEST_CASE("mpole_k1s", "[libFF][Elements]")
{
    propagate_libff("mpole_k1s");
}

TEST_CASE("mpole_k2", "[libFF][Elements]")
{
    propagate_libff("mpole_k2");
}

TEST_CASE("mpole_k2_tilt", "[libFF][Elements]")
{
    propagate_libff("mpole_k2_tilt");
}

TEST_CASE("mpole_k2s", "[libFF][Elements]")
{
    propagate_libff("mpole_k2s");
}

TEST_CASE("mpole_k3", "[libFF][Elements]")
{
    propagate_libff("mpole_k3");
}

TEST_CASE("mpole_k3_tilt", "[libFF][Elements]")
{
    propagate_libff("mpole_k3_tilt");
}

TEST_CASE("mpole_k3s", "[libFF][Elements]")
{
    propagate_libff("mpole_k3s");
}

TEST_CASE("mpole_k4", "[libFF][Elements]")
{
    propagate_libff("mpole_k4");
}

TEST_CASE("mpole_k4_tilt", "[libFF][Elements]")
{
    propagate_libff("mpole_k4_tilt");
}

TEST_CASE("mpole_k4s", "[libFF][Elements]")
{
    propagate_libff("mpole_k4s");
}

TEST_CASE("mpole_k5", "[libFF][Elements]")
{
    propagate_libff("mpole_k5");
}

TEST_CASE("mpole_k5_tilt", "[libFF][Elements]")
{
    propagate_libff("mpole_k5_tilt");
}

TEST_CASE("mpole_k5s", "[libFF][Elements]")
{
    propagate_libff("mpole_k5s");
}
