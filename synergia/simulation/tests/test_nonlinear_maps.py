#!/usr/bin/env python

import pytest
import numpy
import synergia

# This tests that the 3rd order maps for a set of common element types are not
# totally wacky.

# elements_to_test is just being used as a container for a set of lattice elements

elements_to_test = """
! define beam statement so test passes on ubuntu-clang
beam, particle=proton, energy=pmass+0.4;

d: drift, l=1.0;
q: quadrupole, l=1.0, k1=0.1;
s: sextupole, l=0.1, k2=0.1;
o: octupole, l=0.1, k3=0.05;
b: sbend, l=1.0, angle=pi/24;
cf1: sbend, l=2.889612, angle=0.07074218219630160065, k1=0.05410921561;
cf2: sbend, l=2.889612, angle=0.07074218219630160065, e1=0.03537109109815080032, e2=0.03537109109815080032, k1=0.05410921561, k2=-0.006384940;
rfc: rfcavity, l=0.0, volt=0.2, freq=44.0;

elems: line=(d, q, s, o, b, cf1, cf2, rfc);
"""

# cf2: sbend, angle=0.07074218219630160065, e1=0.03537109109815080032, e2=0.03537109109815080032, k1=0.05410921561, k2=-0.006384940;


def test_element_mapping():
    reader = synergia.lattice.MadX_reader()
    reader.parse(elements_to_test)
    testlattice = reader.get_lattice("elems")

    elems = testlattice.get_elements()
    # peel off each element, create a lattice with the element sandwiched between 0 length drifts and try to get maps
    for e in elems:
        refpart = synergia.foundation.Reference_particle(
            1,
            synergia.foundation.pconstants.mp,
            0.8 + synergia.foundation.pconstants.mp,
        )
        lattice = synergia.lattice.Lattice("foo")
        lattice.set_reference_particle(refpart)
        d = synergia.lattice.Lattice_element("drift", "dr1")
        d.set_double_attribute("l", 0.0)
        lattice.append(d)
        lattice.append(e)
        lattice.append(d)

        print("testing element: ", e.get_type_name(), ": ", e.get_name())
        mapping = synergia.simulation.Lattice_simulator.get_one_turn_map_o3(lattice)

        # or iterate through the components and fields
        for comp in range(6):
            trigon = mapping.component(comp)

            for pwr in range(trigon.power() + 1):
                for idx in range(trigon.count(pwr)):
                    # a non-wacko term will have a magnitude < 1e10
                    # this failed before PR#137 fixing the trigon log calculation
                    assert abs(trigon.get_term(pwr, idx)) < 1.0e10

        del mapping
        del lattice
        del d
        del refpart


if __name__ == "__main__":
    test_element_mapping()
