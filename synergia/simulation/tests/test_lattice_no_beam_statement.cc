#include <catch2/catch_test_macros.hpp>

#include "synergia/lattice/lattice.h"
#include "synergia/lattice/lattice_element.h"

#include "synergia/lattice/madx_reader.h"
#include "synergia/simulation/independent_stepper_elements.h"
#include "synergia/simulation/propagator.h"

#include <string>

Lattice
get_lattice()
{
    static std::string fodo_madx(R"foo(
d: drift, l=1.0;
q: quadrupole, l=1.0, k1=0.1;
s: sextupole, l=0.1, k2=0.1;
o: octupole, l=0.1, k3=0.05;
b: sbend, l=1.0, angle=pi/24;
cf1: sbend, l=2.889612, angle=0.07074218219630160065, k1=0.05410921561;
cf2: sbend, l=2.889612, angle=0.07074218219630160065, e1=0.03537109109815080032, e2=0.03537109109815080032, k1=0.05410921561, k2=-0.006384940;
rfc: rfcavity, l=0.0, volt=0.2, freq=44.0;

elems: line=(d, q, s, o, b, cf1, cf2, rfc);
)foo");

    MadX_reader reader;
    reader.parse(fodo_madx);
    return reader.get_lattice("elems");
}

// this lattice failed with the clang build without the llvm unwinder
TEST_CASE("get_lattice")
{
    Lattice lattice(get_lattice());
}

Propagator
create_propagator(Lattice lattice)
{
    Independent_stepper_elements stepper(1);
    Propagator prop(lattice, stepper);
    return prop;
}

TEST_CASE("create_propagator")
{
    Lattice lattice(get_lattice());
    Propagator p(create_propagator(lattice));
}
