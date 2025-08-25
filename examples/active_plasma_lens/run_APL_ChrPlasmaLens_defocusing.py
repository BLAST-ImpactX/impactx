#!/usr/bin/env python3
#
# Copyright 2022-2025 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell, Kyrre Sjobak
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

from run_APL_ChrPlasmaLens import run_APL_ChrPlasmaLens

# Run the ChrPlasmaLens APL test in defocusing mode
# (rigiditiy is negative. Gradient given in [T/m])
run_APL_ChrPlasmaLens(+1000, 1.0e-3, 100e-6)
