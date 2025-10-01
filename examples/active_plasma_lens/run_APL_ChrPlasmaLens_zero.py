#!/usr/bin/env python3
#
# Copyright 2022-2025 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell, Kyrre Sjobak
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

from run_APL import run_APL_tracking

# Run the ChrPlasmaLens APL test in no-field mode
run_APL_tracking(0.0, 1e-3, 10e-6, lensType="ChrPlasmaLens")

# Gives the same output -- for creating the analysis file
# run_APL_tracking(0.0, 1e-3, 10e-6, lensType='ChrDrift')
