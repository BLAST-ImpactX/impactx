#!/usr/bin/env python3
#
# Copyright 2025 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-


from impactx import my_run

verbose = False
# mode = "forward"
mode = "backward"
inputs_file_beam = "examples/fodo_space_charge/input_fodo_envelope_sc.in"


eps = 1e-12

mode = "gradient-free"
print(f"{mode}:")
values = my_run(-85.0, -eps, mode, inputs_file_beam, verbose=verbose)
error_zero_lo = values["error"]
values = my_run(-85.0, eps, mode, inputs_file_beam, verbose=verbose)
error_zero_hi = values["error"]
print(f"  error_zero_lo={error_zero_lo}, error_zero_hi={error_zero_hi}")
print(f"  FD gradient: derror_dq2_k={(error_zero_hi - error_zero_lo) / (2 * eps)}")

values = my_run(-85.0, 0.0, mode, inputs_file_beam, verbose=verbose)
print(f"  error_forward={values['error']}")

mode = "forward"
print(f"\n{mode}:")
values = my_run(-85.0, 0.0, mode, inputs_file_beam, verbose=verbose)
print(f"  error_forward={values['error']} derror_dq2_k={values['derror_dq2_k']}")

mode = "backward"
print(f"\n{mode}:")
values = my_run(-85.0, 0.0, mode, inputs_file_beam, verbose=verbose)
print(f"  error_forward={values['error']} derror_dq2_k={values['derror_dq2_k']}")
