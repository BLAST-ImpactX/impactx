#!/usr/bin/env python3
#
# Copyright 2022-2025 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# This set of tests is used for performance benchmarking of individual
# elements pushes (as micro-benchmarks). We use this file to rapidly
# evaluate performance changes when tuning beamline element performance
# on CPUs and GPUs.
#
# -*- coding: utf-8 -*-

from functools import partial

import pytest

from impactx import ImpactX, Map6x6, distribution, elements, twiss

# benchmark config
rounds = 5
npart = 1_000_000  # increase this to >10M or even 100M to avoid L1/L2/L3 cache effects on some hardware.

# element config
nslice = 1
mapsteps = 4


@pytest.fixture(scope="function")
def sim():
    class SimContextManager:
        def __enter__(self):
            self.sim = ImpactX()

            # set numerical parameters and IO control
            self.sim.particle_shape = 2  # B-spline order
            self.sim.space_charge = False
            self.sim.diagnostics = False  # benchmarking
            self.sim.slice_step_diagnostics = False

            self.sim.init_grids()

            beam = self.sim.particle_container()
            beam.clear_particles()

            # load a 1 GeV electron beam with an initial
            # unnormalized rms emittance of 1 nm
            kin_energy_MeV = 1.0e3  # reference energy
            bunch_charge_C = 1.0e-9  # used with space charge

            #   reference particle
            ref = self.sim.particle_container().ref_particle()
            ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_kin_energy_MeV(
                kin_energy_MeV
            )

            #   particle bunch
            distr = distribution.Waterbag(
                **twiss(
                    beta_x=1.0,
                    beta_y=1.0,
                    beta_t=1.0,
                    emitt_x=1e-09,
                    emitt_y=1e-09,
                    emitt_t=1e-06,
                    alpha_x=0.0,
                    alpha_y=0.0,
                    alpha_t=0.0,
                )
            )
            self.sim.add_particles(bunch_charge_C, distr, npart)
            assert self.sim.particle_container().total_number_of_particles() == npart

            self.sim.backup_beam = self.sim.particle_container().make_alike()
            assert self.sim.backup_beam.total_number_of_particles() == 0

            self.sim.backup_beam.add_particles(
                self.sim.particle_container(), local=True
            )
            assert self.sim.backup_beam.total_number_of_particles() == npart

            return self.sim

        def __exit__(self, exc_type, exc_value, traceback):
            self.sim.finalize()
            del self.sim

    with SimContextManager() as sim:
        yield sim


def pc_setup(sim):
    """Fresh PC for each call of benchmark"""

    assert sim.backup_beam.total_number_of_particles() == npart

    # instead of drawing from the distribution again, create a 2nd
    # particle container and copy the same initial particles again.
    pc = sim.particle_container()
    pc.clear_particles()
    pc.add_particles(sim.backup_beam, local=True)

    assert pc.total_number_of_particles() == npart

    return (pc,), {}


def test_Aperture(benchmark, sim):
    el = elements.Aperture(
        name="collimator", aperture_x=4.0e-5, aperture_y=4.0e-5, shape="rectangular"
    )
    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


def test_Buncher(benchmark, sim):
    el = elements.Buncher(name="shortrf1", V=0.01, k=15.0)
    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


def test_CFbend(benchmark, sim):
    el = elements.CFbend(ds=0.5, rc=7.613657587094493, k=-7.057403, nslice=nslice)
    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


def test_ChrDrift(benchmark, sim):
    chrdrift = elements.ChrDrift(name="drift1", ds=0.25, nslice=nslice)
    benchmark.pedantic(chrdrift.push, setup=partial(pc_setup, sim), rounds=rounds)


def test_ChrPlasmaLens(benchmark, sim):
    el = elements.ChrPlasmaLens(
        name="q1", ds=0.331817852986604588, k=2.98636067687944129, unit=0, nslice=nslice
    )
    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


def test_ChrQuad(benchmark, sim):
    el = elements.ChrQuad(name="quad1", ds=1.0, k=1.0, nslice=nslice)
    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


def test_ChrAcc(benchmark, sim):
    el = elements.ChrAcc(
        name="acc", ds=1.8, ez=10871.950994502130424, bz=1.0e-12, nslice=nslice
    )
    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


def test_ConstF(benchmark, sim):
    el = elements.ConstF(name="constf1", ds=2.0, kx=1.0, ky=1.0, kt=1.0, nslice=nslice)
    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


def test_DipEdge(benchmark, sim):
    rc = 10.3462283686195526  # bend radius (meters)
    psi = 0.048345620280243  # pole face rotation angle (radians)
    el = elements.DipEdge(name="dipedge1", psi=-psi, rc=-rc, g=0.0, K2=0.0)
    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


def test_Drift(benchmark, sim):
    el = elements.Drift(name="drift1", ds=0.25, nslice=nslice)
    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


# def test_Empty(benchmark, sim):
#    el = elements.Empty()
#    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


def test_ExactDrift(benchmark, sim):
    el = elements.ExactDrift(name="drift1", ds=0.25, nslice=nslice)
    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


def test_ExactMultipole(benchmark, sim):
    el = elements.ExactMultipole(
        name="quad1",
        ds=1.0,
        k_normal=[0.0, 1.0],
        k_skew=[0.0, 0.0],
        unit=0,
        mapsteps=mapsteps,
        nslice=nslice,
    )
    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


def test_ExactQuad(benchmark, sim):
    el = elements.ExactQuad(
        name="quad1",
        ds=1.0,
        k=1.0,
        int_order=4,
        nslice=nslice,
        mapsteps=mapsteps,
    )
    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


def test_ExactSbend(benchmark, sim):
    el = elements.ExactSbend(name="bend", ds=1.0, phi=10.0, B=0.0, nslice=nslice)
    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


def test_Kicker(benchmark, sim):
    el = elements.Kicker(name="kick1", xkick=2.0e-3, ykick=0.0, unit="dimensionless")
    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


def test_LinearMap(benchmark, sim):
    R1 = Map6x6.identity()

    ds1 = 0.25
    R1[1, 2] = ds1
    R1[3, 4] = ds1
    R1[5, 6] = ds1 / 16.6464  # ds / (beta*gamma^2)
    el = elements.LinearMap(name="drift1", R=R1, ds=ds1)
    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


# def test_Marker(benchmark, sim):
#    el = elements.Marker(name="marker")
#    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


def test_Multipole(benchmark, sim):
    el = elements.Multipole(
        name="thin_octupole", multipole=4, K_normal=65.0, K_skew=6.0
    )
    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


def test_NonlinearLens(benchmark, sim):
    el = elements.NonlinearLens(name="nllens", knll=4.0e-6, cnll=0.01)
    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


def test_PlaneXYRot(benchmark, sim):
    el = elements.PlaneXYRot(name="rotation1", angle=90.0)
    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


# def test_Programmable(benchmark, sim):
#    el = elements.Programmable(...)
#    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


def test_PRot(benchmark, sim):
    el = elements.PRot(name="rotation1", phi_in=0.0, phi_out=-5.0)
    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


def test_Quad(benchmark, sim):
    el = elements.Quad(name="quad1", ds=1.0, k=1.0, nslice=nslice)
    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


# TODO QuadEdge example missing https://github.com/BLAST-ImpactX/impactx/issues/986
# def test_QuadEdge(benchmark, sim):
#     el = elements.QuadEdge(
#         name="quadedge", ...
#     )
#     benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


def test_RFCavity(benchmark, sim):
    el = elements.RFCavity(
        name="rf",
        ds=1.31879807,
        escale=62.0,
        freq=1.3e9,
        phase=85.5,
        cos_coefficients=[
            0.1644024074311037,
            -0.1324009958969339,
            4.3443060026047219e-002,
            8.5602654094946495e-002,
            -0.2433578169042885,
            0.5297150596779437,
            0.7164884680963959,
            -5.2579522442877296e-003,
            -5.5025369142193678e-002,
            4.6845673335028933e-002,
            -2.3279346335638568e-002,
            4.0800777539657775e-003,
            4.1378326533752169e-003,
            -2.5040533340490805e-003,
            -4.0654981400000964e-003,
            9.6630592067498289e-003,
            -8.5275895985990214e-003,
            -5.8078747006425020e-002,
            -2.4044337836660403e-002,
            1.0968240064697212e-002,
            -3.4461179858301418e-003,
            -8.1201564869443749e-004,
            2.1438992904959380e-003,
            -1.4997753525697276e-003,
            1.8685171825676386e-004,
        ],
        sin_coefficients=[
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ],
        mapsteps=mapsteps,
        nslice=nslice,
    )
    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


def test_Sbend(benchmark, sim):
    el = elements.Sbend(name="sbend1", ds=0.5, rc=-10.346, nslice=nslice)
    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


def test_ShortRF(benchmark, sim):
    el = elements.ShortRF(name="shortrf1", V=1000.0, freq=1.3e9, phase=-89.5)
    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


def test_SoftQuadrupole(benchmark, sim):
    el = elements.SoftQuadrupole(
        name="quad1",
        ds=1.0,
        gscale=1.0,
        cos_coefficients=[2],
        sin_coefficients=[0],
        mapsteps=mapsteps,
        nslice=nslice,
    )
    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


def test_SoftSolenoid(benchmark, sim):
    el = elements.SoftSolenoid(
        name="sol1",
        ds=6.0,
        bscale=1.233482899483985,
        cos_coefficients=[
            0.350807812299706,
            0.323554693720069,
            0.260320578919415,
            0.182848575294969,
            0.106921016050403,
            4.409581845710694e-002,
            -9.416427163897508e-006,
            -2.459452716865687e-002,
            -3.272762575737291e-002,
            -2.936414401076162e-002,
            -1.995780078926890e-002,
            -9.102893342953847e-003,
            -2.456410658713271e-006,
            5.788233017324325e-003,
            8.040408292420691e-003,
            7.480064552867431e-003,
            5.230254569468851e-003,
            2.447614547094685e-003,
            -1.095525090532255e-006,
            -1.614586867387170e-003,
            -2.281365457438345e-003,
            -2.148709081338292e-003,
            -1.522541739363011e-003,
            -7.185505862719508e-004,
            -6.171194824600157e-007,
            4.842109305036943e-004,
            6.874508102002901e-004,
            6.535550288205728e-004,
            4.648795813759210e-004,
            2.216564722797528e-004,
            -4.100982995210341e-007,
            -1.499332112463395e-004,
            -2.151538438342482e-004,
            -2.044590946652016e-004,
            -1.468242784844341e-004,
        ],
        sin_coefficients=[
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ],
        mapsteps=mapsteps,
        nslice=nslice,
    )
    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


def test_Sol(benchmark, sim):
    el = elements.Sol(name="sol1", ds=3.820395, ks=0.8223219329893234)
    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


# def test_Source(benchmark, sim):
#     el = elements.Source(
#         "openPMD", "../solenoid.py/diags/openPMD/monitor.h5"
#     )
#     benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


def test_TaperedPL(benchmark, sim):
    focal_length = 0.5  # focal length in m
    dtaper = 11.488289081903567  # 1/(horizontal dispersion in m)
    dk = 1.0 / focal_length
    el = elements.TaperedPL(name="pl", k=dk, taper=dtaper, unit=0)
    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)


def test_ThinDipole(benchmark, sim):
    el = elements.ThinDipole(name="kick", theta=0.45, rc=1.0)
    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)
