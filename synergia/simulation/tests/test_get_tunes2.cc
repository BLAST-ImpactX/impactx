#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "catch2/matchers/catch_matchers.hpp"
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

TEST_CASE("get_tunes")
{
    /*
        tunes =
       synergia.simulation.Lattice_simulator.calculate_tune_and_cdt(lattice_fixture)
        print('Tunes: x: ', tunes[0], ', y: ', tunes[1])
        assert tunes[0] == pytest.approx(0.16533417, rel=1.0e-4)
        assert tunes[1] == pytest.approx(0.436092762, rel=1.0e-4)
    */
    Logger screen(0, LoggerV::INFO);

    Lattice lattice = get_lattice();

    screen << "Read lattice, length: " << lattice.get_length() << ", "
           << lattice.get_elements().size() << " elements" << std::endl;

    auto tunes(Lattice_simulator::calculate_tune_and_cdt(lattice));

    screen << "x tune: " << std::setprecision(17) << tunes[0] << " should be "
           << 0.16533417 << std::endl;
    screen << "y tune: " << std::setprecision(17) << tunes[1] << " should be "
           << 0.436092762 << std::endl;

    REQUIRE_THAT(tunes[0], Catch::Matchers::WithinRel(0.16533417, 1.0e-4));
    REQUIRE_THAT(tunes[1], Catch::Matchers::WithinRel(0.436092762, 1.0e-4));
}

TEST_CASE("get_chromaticities")
{
    Logger screen(0, LoggerV::INFO);

    Lattice lattice = get_lattice();
    Reference_particle refpart(lattice.get_reference_particle());
    double beta = refpart.get_beta();

    auto chroms(Lattice_simulator::get_chromaticities(lattice));

    screen << "H chrom: " << std::setprecision(17)
           << chroms.horizontal_chromaticity << " should be "
           << -65.68223681 * beta << std::endl;
    screen << "V chromaticity: " << std::setprecision(17)
           << chroms.vertical_chromaticity << " should be "
           << -25.5636675 * beta << std::endl;

    // CHECK( chroms.horizontal_chromaticity ==
    // Approx(-65.68223681*beta).epsilon(1.0e-4) ); CHECK(
    // chroms.vertical_chromaticity == Approx(-25.5636675*beta).epsilon(1.0e-4)
    // );
}
