#!/usr/bin/env python

import pytest

import synergia


@pytest.fixture
def prop_fixture():
    fodo_madx = """
beam, particle=proton,pc=3.0;

o: drift, l=8.0;
f: quadrupole, l=2.0, k1=0.071428571428571425, a1=1.0;
d: quadrupole, l=2.0, k1=-0.071428571428571425;

fodo: sequence, l=20.0, refer=entry;
fodo_1: f, at=0.0;
fodo_2: o, at=2.0;
fodo_3: d, at=10.0;
fodo_4: o, at=12.0;
endsequence;
"""
    reader = synergia.lattice.MadX_reader()
    reader.parse(fodo_madx)
    lattice = reader.get_lattice("fodo")
    lattice.set_all_string_attribute("extractor_type", "libff")
    stepper = synergia.simulation.Independent_stepper_elements(1)
    propagator = synergia.simulation.Propagator(lattice, stepper)

    return {"lattice": lattice, "propagator": propagator}


def test_create_prop(prop_fixture):
    assert True


def test_prop_get_lattice_elements(prop_fixture):
    nlattice_elem = len(prop_fixture["lattice"].get_elements())
    nprop_latt_elem = len(prop_fixture["propagator"].get_lattice_elements())
    assert nlattice_elem == nprop_latt_elem
    assert True


def test_prop_get_lattice(prop_fixture):
    nlattice_elem = len(prop_fixture["lattice"].get_elements())
    nprop_latt_elem = len(prop_fixture["propagator"].get_lattice().get_elements())
    assert nlattice_elem == nprop_latt_elem
    assert True


def test_modify_lattice_elements(prop_fixture):
    # print(prop_fixture['propagator'].get_lattice().get_elements()[0])
    orig_a1 = (
        prop_fixture["propagator"]
        .get_lattice()
        .get_elements()[0]
        .get_double_attribute("a1")
    )
    new_a1 = orig_a1 + 100.0
    prop_fixture["propagator"].get_lattice().get_elements()[0].set_double_attribute(
        "a1", new_a1
    )
    assert (
        prop_fixture["propagator"]
        .get_lattice()
        .get_elements()[0]
        .get_double_attribute("a1")
        == new_a1
    )

    # There is a problem with the iterating off the end of a the pybind11
    # wrapped slices structure on ubuntu-clang so we can't convert it into
    # a python list. Doing so generates an uncatchable SEGV.
    # slices = list(prop_fixture['propagator'].get_lattice_element_slices())
    # assert slices[0].get_lattice_element().get_double_attribute('a1') == new_a1

    # iterate over slices to pick out the first one which is the
    # one that has been modified
    c = 0
    for s in prop_fixture["propagator"].get_lattice_element_slices():
        print("testing slice ", c)
        if c == 0:
            assert s.get_lattice_element().get_double_attribute("a1") == new_a1
            break
        # actually should never get here
        c = c + 1

    assert True


def test_modify_reference_particle_energy(prop_fixture):
    orig_energy = (
        prop_fixture["propagator"]
        .get_lattice()
        .get_reference_particle()
        .get_total_energy()
    )
    new_energy = orig_energy + 2.0
    prop_fixture["propagator"].get_lattice().get_reference_particle().set_total_energy(
        new_energy
    )
    assert (
        prop_fixture["propagator"]
        .get_lattice()
        .get_reference_particle()
        .get_total_energy()
        == new_energy
    )


def test_set_get_checkpoint(prop_fixture):
    init_period = prop_fixture["propagator"].get_checkpoint_period()
    cp_period = init_period + 1
    prop_fixture["propagator"].set_checkpoint_period(cp_period)
    assert cp_period == prop_fixture["propagator"].get_checkpoint_period()


def test_set_get_final_checkpoint(prop_fixture):
    init_fcp = prop_fixture["propagator"].get_final_checkpoint()
    prop_fixture["propagator"].set_final_checkpoint(not init_fcp)
    assert prop_fixture["propagator"].get_final_checkpoint() is not init_fcp
