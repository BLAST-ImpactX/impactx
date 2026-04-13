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

TEST_CASE("get_tunes_dpp_offset")
{
    Logger screen(0, LoggerV::INFO);

    Lattice lattice = get_lattice();
    Reference_particle refpart(lattice.get_reference_particle());
    double beta = refpart.get_beta();

    constexpr double dpp = 4.0e-4;
    auto tunes(Lattice_simulator::calculate_tune_and_cdt(lattice, dpp));

    screen << "x tune: " << std::setprecision(17) << tunes[0] << std::endl;
    screen << "y tune: " << std::setprecision(17) << tunes[1] << std::endl;

    // from outputs of test_get_tunes2
    constexpr double chromx = -60.349667735076068;
    constexpr double chromy = -21.774885376635297;
    constexpr double tunex_dpp0 = 0.16533644789863577;
    constexpr double tuney_dpp0 = 0.4360932921548214;
    // CHECK (tunes[0] == Approx(tunex_dpp0 + chromx*dpp).epsilon(1.0e-8));
    // CHECK (tunes[1] == Approx(tuney_dpp0 + chromy*dpp).epsilon(1.0e-8));
}
