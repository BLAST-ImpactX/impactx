#!/usr/bin/env python

import pytest
import numpy as np
import synergia

mass = 1.0
betagamma = 3/4
gamma = 5/4
energy = mass*gamma

stdz = 1.2
emit = 12.0e-6
alphax = -0.1
betax = 22.0
alphay = 0.05
betay = 5.0
betaz = 900

sxx = emit*betax
sxxp = -alphax/betax*sxx
sxpxp = sxx*(1 + alphax**2)/betax

syy = emit*betay
syyp = -alphay/betay*syy
sypyp = syy * (1 + alphay**2)/betay

szz = stdz**2
szpzp = stdz**2/betaz**2

bunch_covariances = np.array([ [sxx,  sxxp, 0, 0, 0, 0 ],
                         [sxxp, sxpxp, 0, 0, 0, 0 ],
                         [0, 0, syy, syyp, 0, 0 ],
                         [0, 0, syyp, sypyp, 0, 0 ],
                         [0, 0, 0, 0, szz, 0 ],
                         [ 0, 0, 0, 0, 0, szpzp]])

@pytest.fixture
def bunch_fixture():
    refpart = synergia.foundation.Reference_particle(1, mass, energy)
    macroparticles=10
    realparticles=5e10
    bunch = synergia.bunch.Bunch(refpart, macroparticles, realparticles, synergia.utils.Commxx())
    return bunch

@pytest.fixture
def bunch_fixture2():
    refpart = synergia.foundation.Reference_particle(1, mass, energy)
    macroparticles=10000
    realparticles=5e10
    bunch = synergia.bunch.Bunch(refpart, macroparticles, realparticles, synergia.utils.Commxx())

    return bunch

                       
def test_bunch_1(bunch_fixture):
    assert bunch_fixture.get_total_num() == 10
    assert bunch_fixture.get_real_num() == pytest.approx(5e10)
    bunch_fixture.set_real_num(6e10)
    assert bunch_fixture.get_real_num() == pytest.approx(6e10)

def test_bunch_core_diagnostics(bunch_fixture2):
    dist = synergia.foundation.PCG_random_distribution(1234567, synergia.utils.Commxx.World.rank())

    bunch_means = np.zeros(6)
    synergia.bunch.populate_6d( dist, 
        bunch_fixture2, 
        bunch_means,
        bunch_covariances)

    mean = synergia.bunch.Core_diagnostics.calculate_mean(bunch_fixture2)
    std = synergia.bunch.Core_diagnostics.calculate_std(bunch_fixture2, mean)
    mom2 = synergia.bunch.Core_diagnostics.calculate_mom2(bunch_fixture2, mean)

    for i in range(6):
        assert mean[i] == pytest.approx(0.0)

    assert std[0] == pytest.approx(np.sqrt(sxx))
    assert std[1] == pytest.approx(np.sqrt(sxpxp))
    assert std[2] == pytest.approx(np.sqrt(syy))
    assert std[3] == pytest.approx(np.sqrt(sypyp))
    assert std[4] == pytest.approx(np.sqrt(szz))
    assert std[5] == pytest.approx(np.sqrt(szpzp))

    for i in range(6):
        for j in range(6):
            assert mom2[i, j] == pytest.approx(bunch_covariances[i, j])
