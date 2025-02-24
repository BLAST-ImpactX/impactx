"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from . import setup_server

server, state, ctrl = setup_server()

from impactx import Config, ImpactX

from .Analyze.plotsMain import generate_phase_space
from .Input.distributionParameters.distributionMain import distribution_parameters
from .Input.latticeConfiguration.latticeMain import lattice_elements

# Call MPI_Init and MPI_Finalize only once:
if Config.have_mpi:
    from mpi4py import MPI  # noqa


def run_simulation():
    """
    This tests runs a simulation on ImpactX
    based on user inputs from the dashboard.
    """
    sim = ImpactX()

    # space charge selections
    if state.space_charge:
        sim.max_level = state.max_level
        sim.n_cell = state.n_cell
        sim.particle_shape = state.particle_shape
        sim.poisson_solver = state.poisson_solver
        sim.space_charge = state.space_charge
        sim.dynamic_size = state.dynamic_size
        sim.prob_relative = state.prob_relative
        sim.blocking_factor_x = [state.blocking_factor_x]
        sim.blocking_factor_y = [state.blocking_factor_y]
        sim.blocking_factor_z = [state.blocking_factor_z]

        if state.poisson_solver == "multigrid":
            sim.mlmg_relative_tolerance = state.mlmg_relative_tolerance
            sim.mlmg_absolute_tolerance = state.mlmg_absolute_tolerance
            sim.mlmg_max_iters = state.mlmg_max_iters
            sim.mlmg_verbosity = state.mlmg_verbosity
    # csr
    if state.csr:
        sim.csr = state.csr
        sim.csr_bins = state.csr_bins

    sim.particle_shape = state.particle_shape
    sim.slice_step_diagnostics = True
    sim.init_grids()

    # init particle beam
    kin_energy_MeV = state.kin_energy_MeV
    bunch_charge_C = state.bunch_charge_C
    npart = state.npart

    #   reference particle
    pc = sim.particle_container()
    ref = pc.ref_particle()
    ref.set_charge_qe(state.charge_qe).set_mass_MeV(state.mass_MeV).set_kin_energy_MeV(
        kin_energy_MeV
    )

    distribution = distribution_parameters()
    sim.add_particles(bunch_charge_C, distribution, npart)

    lattice_configuration = lattice_elements()

    sim.lattice.extend(lattice_configuration)

    # simulate
    sim.evolve()

    generate_phase_space(pc)

    sim.finalize()
