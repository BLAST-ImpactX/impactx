#!/usr/bin/env python

from synergia.lattice import Lattice, Lattice_element, MadX_reader
from synergia.foundation import Reference_particle
import pytest

lattice_str = """
    start = 0.35;

    x = 1.0;
    y = {1, 2, 3, 4};
    z = "abc";

    o: drift, l=0.2;

    a: quadrupole, l=0.0, k1=x+1.0;
    b: quadrupole, l=o->l;
    c: quadrupole, l=0.4;
    d: quadrupole, l=o->l*4, k1=o->l*5;

    seq: sequence, l=3.0;
    a, at=0;
    b, at=0;
    c, at=start+0.1;
    d, at=0.6, from=c;
    endsequence;
"""


def test_dynamic_lattice():
    reader = MadX_reader()
    reader.parse(lattice_str)

    lattice = reader.get_dynamic_lattice("seq")
    print(lattice)

    elms = lattice.get_elements()

    # a->k1 == 2.0
    assert elms[0].get_double_attribute("k1") == pytest.approx(2.0, 1e-12)

    # set x
    lattice.set_variable("x", 3.0)

    # a->k1 == 4.0 now
    assert elms[0].get_double_attribute("k1") == pytest.approx(4.0, 1e-12)

    # d->l = o->l*4 == 0.8
    assert elms[4].get_double_attribute("l") == pytest.approx(0.8, 1e-12)

    # d->k1 = o->l*5  == 1.0
    assert elms[4].get_double_attribute("k1") == pytest.approx(1.0, 1e-12)

    # set o->l to 0.3
    lattice.get_lattice_tree().set_element_attribute("o", "l", 0.3)

    # d->l = o->l*4 == 1.2
    assert elms[4].get_double_attribute("l") == pytest.approx(1.2, 1e-12)

    # d->k1 = o->l*5  == 1.5
    assert elms[4].get_double_attribute("k1") == pytest.approx(1.5, 1e-12)

    # parsing error
    with pytest.raises(RuntimeError, match="parse"):
        elms[3].set_double_attribute("k1", "o->l*3+")

    # no throw
    elms[3].set_double_attribute("k1", "o->l*3")

    # k1 = 0.3*3 = 0.9
    assert elms[3].get_double_attribute("k1") == pytest.approx(0.9, 1e-12)
