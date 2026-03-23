#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Eric G. Stern
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-


import mpi4py.MPI as MPI
import numpy as np
import synergia
from scipy import constants
from syn2_to_impactx import syn2_to_impactx, unroll_impactx_lattice

from impactx import ImpactX, elements

c = constants.c
e = constants.e
mp = synergia.foundation.pconstants.mp

sim = ImpactX()

# Run parameters

harmonic_number = 84  # harmonic number
enable_rf = True
rf_volt = 200.0e-6  # 200 KV

turns = 100

DEBUG = True

total_booster_charge = 6.7e12
filled_buckets = 81
bunch_charge_C = e * total_booster_charge / filled_buckets

myrank = MPI.COMM_WORLD.rank

lattice_line = "booster"
lattice_file = "sbbooster-cooked.madx"


# ========================================================================

# set the voltage and tune the lattice
# voltage is total voltage in GV


def set_rf(lattice, voltage, harmno, bunch_phase_offset, phase, above_transition=False):

    # above transition, phase needs to be (pi - phase) for longitudinal stability
    if above_transition:
        phase_set = np.pi - phase
    else:
        phase_set = phase

    # Offset phase by how far the bunch has shifted
    phase_set = phase_set + bunch_phase_offset

    if DEBUG and myrank == 0:
        print("setrf: lattice: ", id(lattice))
        print(
            "set_rf: voltage=",
            voltage,
            ", harmno = ",
            harmno,
            ", phase = ",
            phase_set,
            end="",
        )
    # count RF cavities
    cavities = 0
    for elem in lattice.get_elements():
        if elem.get_type() == synergia.lattice.element_type.rfcavity:
            cavities = cavities + 1
    if DEBUG and myrank == 0:
        print(" for ", cavities, " cavities")

    # Set the RF cavity voltage
    for elem in lattice.get_elements():
        if elem.get_type() == synergia.lattice.element_type.rfcavity:
            elem.set_double_attribute("volt", 1000 * voltage / cavities)
            elem.set_double_attribute("lag", phase_set / (2 * np.pi))
            elem.set_double_attribute("harmon", harmno)

    synergia.simulation.Lattice_simulator.tune_circular_lattice(lattice)

    for elem in lattice.get_elements():
        if (
            DEBUG
            and myrank == 0
            and elem.get_type() == synergia.lattice.element_type.rfcavity
        ):
            print("set_rf: ", elem)
            break

    return lattice


# ========================================================================


def get_lattice():
    # read the lattice in from a MadX sequence file
    lattice_raw = synergia.lattice.MadX_reader().get_lattice(lattice_line, lattice_file)
    lattice_raw.set_all_string_attribute("extractor_type", "libff")

    if enable_rf:
        lattice_with_rf = set_rf(
            lattice_raw,
            voltage=rf_volt,
            harmno=harmonic_number,
            bunch_phase_offset=0,
            phase=0.0,
            above_transition=False,
        )

    else:
        # tune the lattice (set frequency)
        synergia.simulation.Lattice_simulator.tune_circular_lattice(lattice_raw)

    return lattice_raw


# ========================================================================


# ========================================================================


def main():

    # Read the lattice
    lattice = get_lattice()
    print(
        "Read lattice, length = {}, {} elements".format(
            lattice.get_length(), len(lattice.get_elements())
        )
    )

    # Assume the lattice sets the reference particle
    refpart = lattice.get_reference_particle()

    energy = refpart.get_total_energy()
    momentum = refpart.get_momentum()
    gamma = refpart.get_gamma()
    beta = refpart.get_beta()

    if myrank == 0:
        print("Beam parameters")
        print("energy: ", energy, "GeV")
        print("momentum: ", momentum, "GeV/c")
        print("gamma: ", gamma)
        print("beta: ", beta)
        print()

    # Get original tunes and chromaticities
    rf_freq = synergia.simulation.Lattice_simulator.get_rf_frequency(lattice)
    (orig_xtune, orig_ytune, orig_orbit_cdt) = (
        synergia.simulation.Lattice_simulator.calculate_tune_and_cdt(lattice)
    )

    chrom = synergia.simulation.Lattice_simulator.get_chromaticities(lattice)
    orig_hchrom = chrom.horizontal_chromaticity
    orig_vchrom = chrom.vertical_chromaticity
    if myrank == 0:
        print("Orbit length: ", beta * orig_orbit_cdt, "m")
        print("RF frequency: ", rf_freq, "Hz")
        print("original horizontal tune: ", orig_xtune)
        print("original vertical tune: ", orig_ytune)
        print("original horizontal chromaticity: ", orig_hchrom)
        print("original vertical chromaticity: ", orig_vchrom)
        print()

    # calculate lattice functions to calculate initial bunch population
    synergia.simulation.Lattice_simulator.CourantSnyderLatticeFunctions(lattice)
    synergia.simulation.Lattice_simulator.calc_dispersions(lattice)

    # lattice functions at the last element gives the
    # values for the overall ring
    lf = lattice.get_elements()[-1].lf
    beta_x = lf.beta.hor
    alpha_x = lf.alpha.hor
    beta_y = lf.beta.ver
    alpha_y = lf.alpha.ver
    psi_x = lf.psi.hor
    psi_y = lf.psi.ver
    disp_x = lf.dispersion.hor
    dprime_x = lf.dPrime.hor

    if myrank == 0:
        print("lattice functions after adjustment")
        print("beta_x: ", beta_x)
        print("alpha_x", alpha_x)
        print("beta_y: ", beta_y)
        print("alpha_y: ", alpha_y)
        print("dispersion_x: ", disp_x)
        print("dispersion prime_x: ", dprime_x)
        print()

    # The rest of this is standard ImpactX

    # set numerical parameters and IO control
    sim.particle_shape = 2  # B-spline order
    sim.space_charge = False
    # sim.diagnostics = False  # benchmarking
    sim.diagnostics = True
    sim.slice_step_diagnostics = True

    # domain decomposition & space charge mesh
    sim.init_grids()

    # Create the ImpactX reference particle
    ref = sim.particle_container().ref_particle()
    ref.set_charge_qe(refpart.get_charge())
    ref.set_mass_MeV(refpart.get_mass() * 1000.0)
    kin_energy = (refpart.get_total_energy() - refpart.get_mass()) * 1000.0
    ref.set_kin_energy_MeV(kin_energy)

    # charge to mass ratio
    qm_eev = 1.0 / (1.0e-9 * mp)

    # insert the converted MAD-X->Synergia lattice
    init_monitor = False
    final_monitor = True
    ix_lattice = syn2_to_impactx(lattice, init_monitor, final_monitor)

    # save the lattice

    if myrank == 0:
        # first the synergia lattice
        with open("booster_synergia_lattice.txt", "w") as f:
            print(lattice, file=f)

        # impactx lattice
        with open("booster_impactx_lattice.txt", "w") as f:
            print(unroll_impactx_lattice(ix_lattice), file=f)

        # save in madx format too
        lattice.export_madx_file("booster_lattice.madx", sanitize=True)

    # the monitor is the last element in the list. Run a 1 turn simulation
    # get the initial distribution with a Source element
    source = elements.Source("openPMD", "booster_particles.h5")
    sim.lattice.append(source)
    # add the monitor element to get the initial distribution
    sim.lattice.append(ix_lattice[-1])

    print("begin lattice")
    for e in sim.lattice:
        print(e)
    print("end lattice")
    sim.periods = 1
    sim.track_particles()

    print("running full lattice")
    # now run the full lattice
    sim.lattice.clear()
    sim.lattice.extend(ix_lattice)
    print("full lattice has ", len(sim.lattice), "elements")

    # run simulation
    if turns > 0:
        sim.periods = turns
        sim.track_particles()

    # clean shutdown
    sim.finalize()
    pass


if __name__ == "__main__":
    main()
    pass
