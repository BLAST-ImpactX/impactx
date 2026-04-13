#include <catch2/catch_test_macros.hpp>

#include "synergia/lattice/madx_reader.h"

#include "synergia/simulation/independent_stepper_elements.h"
#include "synergia/simulation/propagator.h"

struct propagator_fixture {
    Logger screen;
    Lattice lattice;
    Propagator propagator;
    std::unique_ptr<Bunch_simulator> sim;

    propagator_fixture(std::string const& seq)
        : screen(0, LoggerV::INFO_TURN)
        , lattice(MadX_reader().get_lattice(seq, "fodo.madx"))
        , propagator(lattice, Independent_stepper_elements(1))
        , sim()
    {
        auto ref = lattice.get_reference_particle();
        auto fm = ref.get_four_momentum();
        // fm.set_momentum(fm.get_momentum()*0.25);
        fm.set_momentum(fm.get_momentum() * 0.95);
        // fm.set_momentum(3.0);
        ref.set_four_momentum(fm);

        sim = std::make_unique<Bunch_simulator>(
            Bunch_simulator::create_single_bunch_simulator(
                ref, 1, 1e09, Commxx()));

        // sim = Bunch_simulator::create_single_bunch_simulator(
        //         ref, 1, 1e09, Commxx());

        // propagate options (start turn, num turns)
        // sim.set_turns(0, 1);
    }

    void
    propagate()
    {
        propagator.propagate(*sim, screen, 1);
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

void
propagate_libff(std::string const& seq)
{
    std::cout << "libff propagate " << seq << "\n";

    propagator_fixture pf(seq);
    pf.print_lattice();
    auto& b = pf.bunch();

    b.checkout_particles();
    auto parts = b.get_host_particles();
    for (int i = 0; i < 6; ++i)
        parts(0, i) = 0.1;
    b.checkin_particles();

    pf.propagate();

    b.checkout_particles();
    parts = b.get_host_particles();
    std::cout << std::setprecision(16);
    for (int i = 0; i < 6; ++i)
        std::cout << parts(0, i) << "\n";
    std::cout << "\n";
}

TEST_CASE("sbend", "[libFF][Elements]")
{
    CHECK(true);
    propagate_libff("seq_sbend");
}

TEST_CASE("drift", "[libFF][Elements]")
{
    CHECK(true);
    propagate_libff("seq_drift");
}

TEST_CASE("quadrupole", "[libFF][Elements]")
{
    CHECK(true);
    propagate_libff("seq_quadrupole");
}

TEST_CASE("sbend cf", "[libFF][Elements]")
{
    CHECK(true);
    propagate_libff("seq_cfsbend");
}

TEST_CASE("sbend cf 2", "[libFF][Elements]")
{
    CHECK(true);
    propagate_libff("seq_cfsbend2");
}

TEST_CASE("octupole", "[libFF][Elements]")
{
    CHECK(true);
    propagate_libff("seq_octupole");
}

TEST_CASE("octupole with tilt", "[libFF][Elements]")
{
    CHECK(true);
    propagate_libff("seq_octupole2");
}

TEST_CASE("sextupole", "[libFF][Elements]")
{
    CHECK(true);
    propagate_libff("seq_sextupole");
}

TEST_CASE("sextupole with tilt", "[libFF][Elements]")
{
    CHECK(true);
    propagate_libff("seq_sextupole2");
}

TEST_CASE("rfcavity", "[libFF][Elements]")
{
    CHECK(true);
    propagate_libff("seq_rfc");
}

TEST_CASE("thin multipole 1 -- mad8 format", "[libFF][Elements]")
{
    CHECK(true);
    propagate_libff("seq_mp1");
}

TEST_CASE("thin multipole 2 -- mad8 format with tilts", "[libFF][Elements]")
{
    CHECK(true);
    propagate_libff("seq_mp2");
}

TEST_CASE("thin multipole 3 -- madx format", "[libFF][Elements]")
{
    CHECK(true);
    propagate_libff("seq_mp3");
}

TEST_CASE("thin multipole 4 -- madx format with tilts", "[libFF][Elements]")
{
    CHECK(true);
    propagate_libff("seq_mp4");
}

TEST_CASE("hkicker", "[libFF][Elements]")
{
    CHECK(true);
    propagate_libff("seq_hkicker");
}

TEST_CASE("vkicker", "[libFF][Elements]")
{
    CHECK(true);
    propagate_libff("seq_vkicker");
}

TEST_CASE("kicker", "[libFF][Elements]")
{
    CHECK(true);
    propagate_libff("seq_kicker");
}

TEST_CASE("long hkicker", "[libFF][Elements]")
{
    CHECK(true);
    propagate_libff("seq_long_hkicker");
}

TEST_CASE("long hkicker simple", "[libFF][Elements]")
{
    CHECK(true);
    propagate_libff("seq_long_hkicker_simple");
}

TEST_CASE("long vkicker", "[libFF][Elements]")
{
    CHECK(true);
    propagate_libff("seq_long_vkicker");
}

TEST_CASE("long vkicker simple", "[libFF][Elements]")
{
    CHECK(true);
    propagate_libff("seq_long_vkicker_simple");
}

TEST_CASE("long kicker", "[libFF][Elements]")
{
    CHECK(true);
    propagate_libff("seq_long_kicker");
}

TEST_CASE("long kicker simple", "[libFF][Elements]")
{
    CHECK(true);
    propagate_libff("seq_long_kicker_simple");
}

#if 0
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);

    int result = Catch::Session().run(argc, argv);

    Kokkos::finalize();
    MPI_Finalize();

    return result;
}

#endif
