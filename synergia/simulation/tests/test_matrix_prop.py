#!/usr/bin/env python
import sys, os
import numpy as np
import synergia
import pytest

macroparticles = 40  # 8 particles in transverse ring x 5 longitudinal positions
realparticles = 4.0e10

# extracted map from 400 MeV Booster simulation
booster_map = np.array(
    [
        [
            -2.36091394e-01,
            -3.26760289e01,
            0.00000000e00,
            0.00000000e00,
            1.50784719e-04,
            3.94071286e00,
        ],
        [
            2.87151957e-02,
            -2.61249729e-01,
            0.00000000e00,
            0.00000000e00,
            6.41040556e-07,
            -9.02291366e-02,
        ],
        [
            0.00000000e00,
            0.00000000e00,
            3.62462342e-01,
            -4.88366765e00,
            0.00000000e00,
            0.00000000e00,
        ],
        [
            0.00000000e00,
            0.00000000e00,
            1.77021951e-01,
            3.73786762e-01,
            0.00000000e00,
            0.00000000e00,
        ],
        [
            -1.30316509e-01,
            5.56865107e00,
            0.00000000e00,
            0.00000000e00,
            9.96790543e-01,
            3.03415425e02,
        ],
        [
            2.12179267e-07,
            2.34390926e-04,
            0.00000000e00,
            0.00000000e00,
            -4.66014699e-05,
            9.89011964e-01,
        ],
    ]
)


# prop_fixture is a propagator
@pytest.fixture
def pf():
    # booster-like lattice
    channel_madx = """
beam, particle=proton,energy=pmass+0.4;

m: matrix, rm11=-0.2360913935508699, 
rm12=-32.676028880719926, 
rm15=0.00015078471867412094, 
rm16=3.940712863966235, 
rm21=0.02871519571777345, 
rm22=-0.26124972905688154, 
rm25=6.410405557641018e-07, 
rm26=-0.09022913657323245, 
rm33=0.36246234212348744, 
rm34=-4.883667648489127, 
rm43=0.17702195094897383, 
rm44=0.3737867616380023, 
rm51=-0.1303165086379494, 
rm52=5.56865106574654, 
rm55=0.9967905427138485, 
rm56=303.41542490997386, 
rm61=2.1217926701948382e-07, 
rm62=0.00023439092605032628, 
rm65=-4.660146987991736e-05, 
rm66=0.9890119640563689;

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


def create_simulator(ref_part):
    sim = synergia.simulation.Bunch_simulator.create_single_bunch_simulator(
        ref_part, macroparticles, realparticles
    )
    bunch = sim.get_bunch()
    s2o2 = np.sqrt(2.0) / 2
    bunch.checkout_particles()
    lp = bunch.get_particles_numpy()

    dr = 0.001

    # populate longitudinal coordinates first
    k = 0
    for iz in range(-2, 3):
        # loop over longitudinal position
        cdt = 0.5 * iz
        for j in range(8):
            lp[k, 0:6] = 0.0
            lp[k, 4] = cdt
            k = k + 1

    # populate transverse coordinates
    k = 0
    for j in range(5):
        lp[k, 0] = dr
        k = k + 1

        lp[k, 0] = s2o2 * dr
        lp[k, 2] = s2o2 * dr
        k = k + 1

        lp[k, 2] = dr
        k = k + 1

        lp[k, 0] = -s2o2 * dr
        lp[k, 2] = s2o2 * dr
        k = k + 1

        lp[k, 0] = -dr
        k = k + 1

        lp[k, 0] = -s2o2 * dr
        lp[k, 2] = -s2o2 * dr
        k = k + 1

        lp[k, 2] = -dr
        k = k + 1

        lp[k, 0] = s2o2 * dr
        lp[k, 2] = -s2o2 * dr
        k = k + 1

    bunch.checkin_particles()
    return sim


def test_lattice_prop(pf):
    refpart = pf.get_lattice().get_reference_particle()
    sim = create_simulator(refpart)

    bunch = sim.get_bunch(0, 0)
    bunch.checkout_particles()
    lp = bunch.get_particles_numpy()

    # copy original particles
    op = np.array(lp)[:, :6]

    simlog = synergia.utils.Logger()
    # propagate 1 turn
    pf.propagate(sim, simlog, 1)

    bunch.checkout_particles()
    prop_particles = bunch.get_particles_numpy()

    # multiply the map with the original particles
    # particle array is Nx6 so to propagate it with the map it would have
    # to be transposed to 6xN, We can instead transpose the map.
    op_by_map = np.dot(op, booster_map.transpose())

    for i in range(40):
        for j in range(6):
            # print(i, j)
            assert op_by_map[i, j] == pytest.approx(prop_particles[i, j])


def main():
    pf = pf()
    test_lattice_prop(pf)


if __name__ == "__main__":
    main()
