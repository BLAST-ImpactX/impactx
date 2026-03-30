#include <catch2/catch_test_macros.hpp>

#include "synergia/lattice/madx_reader.h"

#include "synergia/simulation/independent_stepper_elements.h"
#include "synergia/simulation/propagator.h"

#include "synergia/bunch/populate.h"

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

#include <synergia/foundation/physical_constants.h>
#include <synergia/libFF/ff_element.h>

TEST_CASE("foil")
{
    auto mass = pconstants::mp;
    auto charge = pconstants::proton_charge;
    auto energy = 1.0;

    Reference_particle ref(charge, mass, energy);
    std::cout << "pref = " << ref.get_momentum()
              << ", beta = " << ref.get_beta() << "\n";

    // comm world
    Commxx comm;

    // trigon bunch to get the one-turn-map
    bunch_t<double> tb(ref, 1e6, 1.0e13, comm);

#if 0
    auto tparts = tb.get_host_particles();

    // init value
    tparts(0, 0) = 0.0;
    tparts(0, 1) = 0.0;
    tparts(0, 2) = 0.0;
    tparts(0, 3) = 0.0;
    tparts(0, 4) = 0.0;
    tparts(0, 5) = 0.0;

    // check in
    tb.checkin_particles();
#endif

#if 1
    karray1d means("means", 6);
    for (int i = 0; i < 6; ++i)
        means(i) = 0.0;

    karray2d_row covariances("covariances", 6, 6);
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            covariances(i, j) = 0.0;

    covariances(0, 0) = 1e-6;
    covariances(1, 1) = 1e-6;
    covariances(2, 2) = 1e-6;
    covariances(3, 3) = 1e-6;
    covariances(4, 4) = 1e-6;
    covariances(5, 5) = 1e-6;

    Random_distribution dist(5, Commxx::world_rank());
    populate_6d(dist, tb, means, covariances);
#endif

    Lattice_element ele("foil", "f");

    ele.set_double_attribute("xmin", -0.05);
    ele.set_double_attribute("xmax", 0.05);
    ele.set_double_attribute("ymin", -0.05);
    ele.set_double_attribute("ymax", 0.05);
    ele.set_double_attribute("thick", 600);
    ele.set_double_attribute("simple", 1);

    FF_element::apply(ele, tb);

    Lattice_element quad("quadrupole", "q");
    ele.set_double_attribute("l", 1.0);
    ele.set_double_attribute("k1", 0.3);

    FF_element::apply(quad, tb);

    // checkout particles
    tb.checkout_particles();

    // print timer
    Logger screen(0, LoggerV::DEBUG);
    simple_timer_print(screen);
}
