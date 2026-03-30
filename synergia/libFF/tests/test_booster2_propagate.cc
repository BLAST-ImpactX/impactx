#include <catch2/catch_test_macros.hpp>

#include "synergia/lattice/madx_reader.h"

#include "synergia/simulation/independent_stepper_elements.h"
#include "synergia/simulation/propagator.h"

// Fixture that takes a string definition of an element,
// stuffs it into a lattice and creates a simulator/bunch ready to
// propagate, and a propagator to do the propagation.

struct propagator_fixture {
    Logger screen;
    Lattice lattice;
    Propagator propagator;
    std::unique_ptr<Bunch_simulator> sim;

    propagator_fixture()
        : lattice(propagator_fixture::construct_test_lattice())
        , screen(0, LoggerV::INFO_TURN)
        , propagator(Propagator(lattice, Independent_stepper_elements(1)))
        , sim()
    {

        const Reference_particle& refpart = lattice.get_reference_particle();

        sim = std::make_unique<Bunch_simulator>(
            Bunch_simulator::create_single_bunch_simulator(
                refpart, 8, 1e09, Commxx()));
    }

    Lattice
    construct_test_lattice()
    {
        std::string booster_madx(R"foo(
ncells=24;
turn_voltage=1.0; ! 1 MV /turn
beam, particle=proton,energy=pmass+0.8;

b: sbend, l=2.0, angle=(pi/(2*ncells));
f: quadrupole, l=2.0, k1=1/16.2;
d: quadrupole, l=2.0, k1=-1/16.7;
rfc: rfcavity, l=0.0, volt=turn_voltage/ncells, harmon=96, lag=(1/120.0);

cell: sequence, l=20.0, refer=centre;
fodo_1: f, at=1.5;
fodo_1a: b, at=5.0;
fodo_2: d, at=8.5;
fodo_3: d, at=11.5;
fodo_3a: b, at=15.0;
fodo_4: f, at=18.5;
fodo_5: rfc, at=20.0;
endsequence;

booster: sequence, l=480.0, refer=entry;
cell, at=0.0;
cell, at=20.0;
cell, at=40.0;
cell, at=60.0;
cell, at=80.0;
cell, at=100.0;
cell, at=120.0;
cell, at=140.0;
cell, at=160.0;
cell, at=180.0;
cell, at=200.0;
cell, at=220.0;
cell, at=240.0;
cell, at=260.0;
cell, at=280.0;
cell, at=300.0;
cell, at=320.0;
cell, at=340.0;
cell, at=360.0;
cell, at=380.0;
cell, at=400.0;
cell, at=420.0;
cell, at=440.0;
cell, at=460.0;
endsequence;)foo");

        // construct the lattice from the given element by sandwiching
        // it between the lattice_head and lattice_tail
        MadX_reader reader;
        reader.parse(booster_madx);
        Lattice lattice(reader.get_lattice("booster"));
        // Turn off RF cavities the same way it is done in
        // calculate_closed_orbit
        for (auto& ele : lattice.get_elements())
            if (ele.get_type() == element_type::rfcavity) {
                ele.set_double_attribute("volt", 0.0);
                ele.set_double_attribute("lag", 0.0);
                ele.set_double_attribute("freq", 0.0);
            }

        return lattice;
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
propagate_test_elem()
{
    propagator_fixture pf;

    auto& b = pf.bunch();

    b.checkout_particles();
    auto parts = b.get_host_particles();
    for (int i = 0; i < 6; ++i)
        parts(0, i) = 0.0;
    b.checkin_particles();

    pf.propagate();

    b.checkout_particles();
    parts = b.get_host_particles();
    std::cout << std::scientific << std::setprecision(16);
    for (int i = 0; i < 6; ++i)
        std::cout << parts(0, i) << "\n";
    std::cout << "\n";

    // For propagating particle at 0, all the transverse elements should
    // remain at 0

    // on Ryzen 7, x hits 7.82e-17
    for (int i = 0; i < 4; ++i) {
        CHECK(std::abs(parts(0, i)) < 1.0e-16);
    }
}

TEST_CASE("booster_noCF")
{
    propagate_test_elem();
}
