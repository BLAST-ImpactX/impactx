#!/usr/bin/env python

import pytest
import numpy
import synergia

mass = 1.0
betagamma = 3/4
gamma = 5/4
energy = mass*gamma

@pytest.fixture
def bunch_fixture():
    refpart = synergia.foundation.Reference_particle(1, mass, energy)
    macroparticles=10
    realparticles=5e10
    bunch = synergia.bunch.Bunch(refpart, macroparticles, realparticles, synergia.utils.Commxx())
    return bunch

def test_bunch_refpart(bunch_fixture):
    assert bunch_fixture.get_reference_particle().get_total_energy() == energy
    new_energy = energy+1.0
    bunch_fixture.get_reference_particle().set_total_energy(new_energy)
    assert bunch_fixture.get_reference_particle().get_total_energy() == new_energy

def test_bunch_design_refpart(bunch_fixture):
    assert bunch_fixture.get_design_reference_particle().get_total_energy() == energy
    new_energy = energy+1.0
    bunch_fixture.get_design_reference_particle().set_total_energy(new_energy)
    # the energy should not have changed because this is not a binding-by-reference
    assert bunch_fixture.get_design_reference_particle().get_total_energy() == new_energy

