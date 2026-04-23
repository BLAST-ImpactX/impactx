#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Eric G. Stern, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import argparse
import matplotlib.pyplot as plt

from booster_impactx_lattice import get_lattice

from impactx import ImpactX, elements

# options to run this script
parser = argparse.ArgumentParser(
    description="Plot the Fermilab Booster lattice survey."
)
parser.add_argument(
    "--save-png", action="store_true", help="non-interactive run: save to PNGs"
)
args = parser.parse_args()

# reference particle (800 MeV proton, as in run_simple_booster.py)
sim = ImpactX()
sim.init_grids()
ref = sim.particle_container().ref_particle()
ref.set_species("proton").set_kin_energy_MeV(800.0)

# load the Booster lattice
lattice = elements.KnownElementsList()
lattice.extend(get_lattice())

# survey plot
fig = plt.figure(figsize=(12, 4.8))
lattice.plot_survey(ref=ref)
plt.tight_layout()
if args.save_png:
    plt.savefig("simple_booster_survey.png")
else:
    plt.show()

# clean shutdown
sim.finalize()
