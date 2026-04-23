#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Eric G. Stern
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-


import mpi4py.MPI as MPI
import numpy as np
from scipy import constants
from impactx import ImpactX, distribution, elements, synmadx, twiss
from syn2_to_impactx import syn2_to_impactx, unroll_impactx_lattice

pi = constants.pi
c = constants.c
eV = constants.eV
mp = 1.0e-9 * constants.m_p * c**2/eV # proton mass in GeV

# 
# these lattice functions are calculated with Synergia3
# from the sbbooster-cooked.madx file.
alpha_x = -1.298673960026007664e-02
beta_x = 3.373645362843065243e+01
alpha_y = 6.089861210659328755e-03
beta_y = 5.252517912567207681e+00

disp_x = 3.187407765856291153e+00
disp_px = 1.136005067625678322e-03

# emittances from PIP-II CDR v0.3
emit_x = 16.0e-6 # normalized 95% emit
emit_y = 16.0e-6
emit_eV_s = 0.1 # longitudinal emittance 97% eV-s

sim = ImpactX()

# Run parameters

harmonic_number = 84  # harmonic number
enable_rf = True
rf_volt = 200.0e-6  # 200 KV

# number of particles
npart = 10000

turns = 1

DEBUG = True

total_booster_charge = 6.7e12
filled_buckets = 81
bunch_charge_C = eV * total_booster_charge / filled_buckets

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
        if elem.get_type() == synmadx.element_type.rfcavity:
            cavities = cavities + 1
    if DEBUG and myrank == 0:
        print(" for ", cavities, " cavities")

    # set RF frequency assuming the closed orbit length matches the lattice lengh
    # First I need the velocity
    beta = lattice.get_reference_particle().get_beta()

    lattice_length = lattice.get_length()
    rf_freq = harmno * beta * c/lattice_length

    # Set the RF cavity voltage distributing over all RF cavities

    for elem in lattice.get_elements():
        if elem.get_type() == synmadx.element_type.rfcavity:
            elem.set_double_attribute("volt", 1000 * voltage / cavities)
            elem.set_double_attribute("lag", phase_set / (2 * np.pi))
            elem.set_double_attribute("harmon", harmno)
            elem.set_double_attribute("freq", rf_freq*1.0e-6) # MAD-X convention frequency in MHz

    for elem in lattice.get_elements():
        if (
            DEBUG
            and myrank == 0
            and elem.get_type() == synmadx.element_type.rfcavity
        ):
            print("set_rf: ", elem)
            break

    return lattice


# ========================================================================


def get_lattice():
    # read the lattice in from a MadX sequence file
    lattice_raw = synmadx.MadX_reader().get_lattice(lattice_line, lattice_file)

    if enable_rf:
        lattice_with_rf = set_rf(
            lattice_raw,
            voltage=rf_volt,
            harmno=harmonic_number,
            bunch_phase_offset=0,
            phase=0.0,
            above_transition=False,
        )

    return lattice_raw

# ========================================================================

def main():

    # Read the lattice
    lattice = get_lattice()

    if myrank == 0:
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

    # Create distribution and add particles
    distr = distribution.Gaussian(
        **twiss(
            beta_x = beta_x,
            alpha_x = alpha_x,
            emitt_x = emit_x/(ref.beta_gamma*6), # normalized 95% emit -> geometric
            beta_y = beta_y,
            alpha_y = alpha_y,
            emitt_y = emit_y/(ref.beta_gamma*6), # normalized 95% emit -> geometric
            emitt_t = emit_eV_s * 1.0e-6 * c/(ref.mass_MeV * ref.beta_gamma * 6 * pi), # 97% emit eV-s -> RMS emit
            beta_t = 1258.0, # seems to go from 1158 close to the center to
                            # 1258 at about 1.25m
            # dispersion from Synergia so it needs conversion from dp/p to dE/p
            dispersion_x = disp_x/ref.beta,
            dispersion_px = disp_px/ref.beta
            )
        )

    sim.add_particles(bunch_charge_C, distr, npart)

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

    # add the monitor element to get the initial distribution
    sim.lattice.append(ix_lattice[-1])
    sim.periods = 1
    sim.track_particles()

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

if __name__ == "__main__":
    main()

