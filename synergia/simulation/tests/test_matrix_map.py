#!/usr/bin/env python
import sys, os
import numpy as np
import synergia
import pytest

macroparticles = 16
realparticles = 4.0e10

# fake twiss parameters
beta1 = 33
alpha1 = 0.01
mu1 = 2 * np.pi * 0.63

m11 = np.cos(mu1) + alpha1 * np.sin(mu1)
m12 = beta1 * np.sin(mu1)
m21 = -np.sin(mu1) * (1 + alpha1**2) / beta1
m22 = np.cos(mu1) - alpha1 * np.sin(mu1)

beta2 = 12
alpha2 = -0.4
mu2 = 2 * np.pi * 0.74

m33 = np.cos(mu2) + alpha1 * np.sin(mu2)
m34 = beta2 * np.sin(mu2)
m43 = -np.sin(mu2) * (1 + alpha2**2) / beta2
m44 = np.cos(mu2) - alpha2 * np.sin(mu2)

mus = 2 * np.pi * 0.003
betas = 300
m55 = np.cos(mus)
m56 = betas * np.sin(mus)
m65 = -np.sin(mus) / betas
m66 = np.cos(mus)

nturns = 100


# prop_fixture is a propagator
@pytest.fixture
def prop_fixture():
    # booster-like lattice
    channel_madx = """
beam, particle=proton,energy=pmass+0.8;

! reproduce in madx file

beta1 = 33;
alpha1 = 0.01;
mu1 = 2*pi*0.63;

m11 = cos(mu1) + alpha1*sin(mu1);
m12 = beta1 * sin(mu1);
m21 = -sin(mu1) * (1 + alpha1^2)/beta1;
m22 = cos(mu1) - alpha1*sin(mu1);

beta2 = 12;
alpha2 = -.4;
mu2 = 2*pi*0.74;

m33 = cos(mu2) + alpha1*sin(mu2);
m34 = beta2 * sin(mu2);
m43 = -sin(mu2) * (1 + alpha2^2)/beta2;
m44 = cos(mu2) - alpha2*sin(mu2);

mus = 2*pi*0.003;
betas = 300;
m55 = cos(mus);
m56 = betas*sin(mus);
m65 = -sin(mus)/betas;
m66 = cos(mus);

m: matrix, l=0, rm11=m11, rm12=m12, rm21=m21, rm22=m22,
           rm33=m33, rm34=m34, rm43=m43, rm44=m44,
           rm55=m55, rm56=m56, rm65=m65, rm66=m66;
channel: sequence, l=0.0;
m, at=0.0;
endsequence;

"""

    reader = synergia.lattice.MadX_reader()
    reader.parse(channel_madx)
    lattice = reader.get_lattice("channel")
    print(lattice)
    lattice.set_all_string_attribute("extractor_type", "libff")
    stepper = synergia.simulation.Independent_stepper_elements(1)
    propagator = synergia.simulation.Propagator(lattice, stepper)

    return propagator


def test_lattice_map(prop_fixture):
    lattice = prop_fixture.get_lattice()
    map = synergia.simulation.Lattice_simulator.get_linear_one_turn_map(lattice)
    # print(map)
    assert map[0, 0] == pytest.approx(m11)
    assert map[0, 1] == pytest.approx(m12)
    assert map[1, 0] == pytest.approx(m21)
    assert map[1, 1] == pytest.approx(m22)
    assert map[2, 2] == pytest.approx(m33)
    assert map[2, 3] == pytest.approx(m34)
    assert map[3, 2] == pytest.approx(m43)
    assert map[3, 3] == pytest.approx(m44)
    assert map[4, 4] == pytest.approx(m55)
    assert map[4, 5] == pytest.approx(m56)
    assert map[5, 4] == pytest.approx(m65)
    assert map[5, 5] == pytest.approx(m66)

    # all the others should be 0
    for i in range(6):
        for j in range(6):
            if (i == 0 or i == 1) and (j == 0 or j == 1):
                continue
            if (i == 2 or i == 3) and (j == 2 or j == 3):
                continue
            if (i == 4 or i == 5) and (j == 4 or j == 5):
                continue
            assert map[i, j] == pytest.approx(0.0)


def create_simulator(ref_part):
    sim = synergia.simulation.Bunch_simulator.create_single_bunch_simulator(
        ref_part, macroparticles, realparticles
    )
    bunch = sim.get_bunch()
    bunch.checkout_particles()
    lp = bunch.get_particles_numpy()
    lp[:, 0:6] = 0.0
    bunch.checkin_particles()
    return sim


def main():
    pf = prop_fixture()
    test_lattice_map(pf)


if __name__ == "__main__":
    main()
