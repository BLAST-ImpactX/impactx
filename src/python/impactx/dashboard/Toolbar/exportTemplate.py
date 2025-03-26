"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from .. import setup_server
from ..Input.distributionParameters.distributionFunctions import DistributionFunctions
from ..Input.latticeConfiguration.latticeMain import parameter_input_checker_for_lattice

server, state, ctrl = setup_server()

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------


def build_distribution_list():
    """
    Generates an instance of distribution inputs
    as a string for exporting purposes.
    """
    distribution_name = state.distribution
    parameters = DistributionFunctions.convert_distribution_parameters_to_valid_type()

    indentation = " " * (8 if state.distribution_type == "Twiss" else 4)
    distribution_parameters = ",\n".join(
        f"{indentation}{key}={value}" for key, value in parameters.items()
    )

    if state.distribution_type == "Twiss":
        return (
            f"distr = distribution.{distribution_name}(\n"
            f"    **twiss(\n"
            f"{distribution_parameters},\n"
            f"    )\n"
            f")"
        )
    else:
        return (
            f"distr = distribution.{distribution_name}(\n{distribution_parameters},\n)"
        )


def build_lattice_list():
    """
    Generates a string representation of lattice element
    inputs for export purposes.
    """

    lattice_elements = ",\n    ".join(
        f"elements.{element['name']}("
        + ", ".join(
            f"{key}={value}"
            for key, value in parameter_input_checker_for_lattice(element).items()
        )
        + ")"
        for element in state.selected_lattice_list
    )

    return f"lattice_configuration = [\n    {lattice_elements}\n]"


def build_space_charge_or_csr():
    """
    Generates simulation content for space charge
    and csr.
    """
    content = ""

    if state.space_charge:
        content += f"""# Space Charge
sim.space_charge = {state.space_charge}
sim.dynamic_size = {state.dynamic_size}
sim.poisson_solver = '{state.poisson_solver}'
sim.particle_shape = {state.particle_shape}
sim.max_level = {state.max_level}
sim.n_cell = {state.n_cell}
sim.blocking_factor_x = {state.blocking_factor_x}
sim.blocking_factor_y = {state.blocking_factor_y}
sim.blocking_factor_z = {state.blocking_factor_z}
sim.prob_relative = {state.prob_relative}
"""
        if state.poisson_solver == "multigrid":
            content += f"""
# Space Charge - Multigrid-Specific Numerical Options
sim.mlmg_relative_tolerance = {state.mlmg_relative_tolerance}
sim.mlmg_absolute_tolerance = {state.mlmg_absolute_tolerance}
sim.mlmg_max_iters = {state.mlmg_max_iters}
sim.mlmg_verbosity = {state.mlmg_verbosity}
"""
    if state.csr:
        content += f"""# Coherent Synchrotron Radiation
sim.csr = {state.csr}
sim.csr_bins = {state.csr_bins}
"""
        if not state.space_charge:
            content += f"""
sim.particle_shape = {state.particle_shape}
"""
    if not content:
        content = f"""
sim.particle_shape = {state.particle_shape}
"""

    return content


# -----------------------------------------------------------------------------
# Trame setup
# -----------------------------------------------------------------------------


def input_file():
    """
    This function creates the template to export
    dashboard user inputs into a python script.
    """
    script = f"""
from impactx import ImpactX, distribution, elements, twiss

sim = ImpactX()

{build_space_charge_or_csr()}
sim.slice_step_diagnostics = True

sim.init_grids()

# Initialize particle beam
kin_energy_MeV = {state.kin_energy_MeV}
bunch_charge_C = {state.bunch_charge_C}
npart = {state.npart}

# Reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe({state.charge_qe}).set_mass_MeV({state.mass_MeV}).set_kin_energy_MeV(kin_energy_MeV)

{build_distribution_list()}
sim.add_particles(bunch_charge_C, distr, npart)

{build_lattice_list()}
sim.lattice.extend(lattice_configuration)

# Simulate
sim.track_particles()

# Clean Shutdown
sim.finalize()
"""

    return script
