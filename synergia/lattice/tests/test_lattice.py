#!/usr/bin/env python
import pytest
import synergia
from synergia.lattice import Lattice, Lattice_element
from synergia.foundation import Reference_particle

name = "foo"
charge = -1
mass = 100.0
total_energy = 125.0
tolerance = 1.0e-12


def test_construct():
    lattice = Lattice(name)


def test_get_name():
    lattice = Lattice(name)
    assert lattice.get_name() == name


def test_set_reference_particle():
    lattice = Lattice(name)
    reference_particle = Reference_particle(charge, mass, total_energy)
    lattice.set_reference_particle(reference_particle)


def test_get_reference_particle():
    lattice = Lattice(name)
    reference_particle = Reference_particle(charge, mass, total_energy)
    lattice.set_reference_particle(reference_particle)
    assert lattice.get_reference_particle().equal(reference_particle, tolerance)


def test_throws_without_reference_particle():
    lattice = Lattice(name)
    with pytest.raises(RuntimeError):
        lattice.get_reference_particle()


def test_set_reference_particle_energy():
    lattice = Lattice(name)
    reference_particle = Reference_particle(charge, mass, total_energy)
    lattice.set_reference_particle(reference_particle)
    lattice.get_reference_particle().set_total_energy(total_energy + 1)
    assert lattice.get_reference_particle().get_total_energy() == total_energy + 1.0

def test_copy_constructor():
    lattice_txt = """
q1: quadrupole, l=0.5, k1=0.2;
b1: sbend, angle=0.05, l=2.0;
machine: sequence, refer=entry, l=6.0;
    q1, at=1.0;
    b1, at=2.5;
    endsequence;
    beam, particle=proton, energy=1.5;
"""
    
    reader = synergia.lattice.MadX_reader()
    reader.parse(lattice_txt)
    lattice1 = reader.get_lattice('machine')
    assert lattice1.get_length() == 6.0

    lattice2 = lattice1
    assert lattice2.get_length() == 6.0

    # Are lattice1 and lattice2 the same object?
    assert id(lattice1) == id(lattice2)
    # Change to an element in lattice2 is reflected in lattice1
    elem = lattice2.get_elements()[1]
    elem.set_double_attribute('xyzzy', 3.0)
    assert lattice1.get_elements()[1].get_double_attribute('xyzzy') == 3.0

    # lattice3 should be a seperate object with copies of the elements
    lattice3 = synergia.lattice.Lattice(lattice1)
    assert lattice3.get_length() == 6.0

    assert id(lattice1) != id(lattice3)

    # change parameters in an element of lattice3 to check that it doesn't
    # affect lattice1
    elem = lattice3.get_elements()[0]
    elem.set_double_attribute('foo', 4)
    assert not lattice1.get_elements()[0].has_double_attribute('foo')
