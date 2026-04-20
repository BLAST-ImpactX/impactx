#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Eric G. Stern, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import math

from booster_impactx_lattice import get_lattice

from impactx import elements

# load the Booster lattice
lattice = elements.KnownElementsList()
lattice.extend(get_lattice())

# total lattice length (nominal Fermilab Booster circumference: 474.20 m)
total_s = sum(element.ds for element in lattice)
expected_s = 474.20
# tolerance loose enough to hold in single precision accumulation across the ring
assert math.isclose(total_s, expected_s, rel_tol=1.0e-4), (
    f"lattice length {total_s} m does not match expected {expected_s} m"
)

# initial/final beam (TODO)
# series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
# last_step = list(series.iterations)[-1]
# initial = series.iterations[1].particles["beam"].to_df()
# final = series.iterations[last_step].particles["beam"].to_df()

# compare number of particles (TODO)
# num_particles = 10000
# assert num_particles == len(initial)
# assert num_particles == len(final)

# TODO for Eric: add more tests as needed, e.g., checking beam moments
