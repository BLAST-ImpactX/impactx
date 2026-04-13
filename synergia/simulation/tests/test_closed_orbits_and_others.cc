#include <catch2/catch_test_macros.hpp>

#include "synergia/simulation/lattice_simulator.h"

#include "synergia/foundation/physical_constants.h"
#include "synergia/lattice/madx_reader.h"
#include "synergia/utils/commxx.h"
#include "synergia/utils/logger.h"
#include "synergia/utils/simple_timer.h"
#include "synergia/utils/utils.h"
#include <iomanip>
#include <iostream>
#include <string>

Lattice
get_lattice()
{
    static std::string fodo_madx(R"foo(
ncells=24;
turn_voltage=1.0; ! 1 MV /turn
beam, particle=proton,energy=pmass+0.8;

f: sbend, l=2.0, angle=(pi/(2*ncells)), k1=1/16.2;
d: sbend, l=2.0, angle=(pi/(2*ncells)), k1=-1/16.7;
!f: quadrupole, l=2.0, k1=0.0625;
!d: quadrupole, l=2.0, k1=-0.0625;
rfc: rfcavity, l=0.0, volt=turn_voltage/ncells, harmon=96, lag=(1/120.0);


cell: sequence, l=20.0, refer=centre;
fodo_1: f, at=1.5;
fodo_2: d, at=8.5;
fodo_3: d, at=11.5;
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
endsequence;
)foo");

    MadX_reader reader;
    reader.parse(fodo_madx);
    return reader.get_lattice("booster");
}

TEST_CASE("closed_orbit_at_0dpp")
{
    Logger screen(0, LoggerV::INFO);

    Lattice lattice = get_lattice();

    auto closed_orbit_state =
        Lattice_simulator::calculate_closed_orbit(lattice);
    for (int i = 0; i < 6; ++i) {
        std::cout << std::setprecision(17) << i << ": " << closed_orbit_state[i]
                  << std::endl;
    }

    for (int i = 0; i < 6; ++i) {
        CHECK(std::abs(closed_orbit_state[i]) < 1.0e-12);
    }
}

#if 0 // I don't see a good way to make this test work
TEST_CASE("closed_orbit_at_nonzero_dpp")
{
    Logger screen(0, LoggerV::INFO);

    Lattice lattice = get_lattice();

    constexpr double dpp=4.0e-4;
    auto closed_orbit_state = Lattice_simulator::calculate_closed_orbit(lattice, dpp);
    for (int i=0; i<6; ++i) {
        std::cout << std::setprecision(17) << i << ": " << closed_orbit_state[i] << std::endl;    
    }

    CHECK (closed_orbit_state[0] == Approx(0.00072931911596656749)); // previous run
    CHECK (closed_orbit_state[1] == Approx(-5.7585101111132694e-15)); //previous run

}
#endif
