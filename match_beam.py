#!/usr/bin/env python3
#
# Copyright 2025 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import time

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

from impactx import my_run

verbose = False
# mode = "forward"
mode = "backward"
inputs_file_beam = "examples/fodo_space_charge/input_fodo_envelope_sc.in"

history = {}
runtime = {}
optimizer = None
grads = False


class ProcessTimer:
    def __init__(self):
        self.elapsed_time = 0

    def __enter__(self):
        self.start_time = time.process_time_ns()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        end_time = time.process_time_ns()
        self.elapsed_time = end_time - self.start_time


def objective(parameters: tuple) -> float:
    """
    A function that is evaluated by the optimizer.

    Parameters
    ----------
    parameters: tuple
      quadrupole strengths k of quad 1 and quad 2.

    Returns
    -------
    The L2 norm of alpha and beta of the beam at the end of the
    simulation.
    """
    global grads

    if verbose:
        print(f"Run objective with parameters={parameters}...")

    if not grads:
        use_mode = "gradient-free"
    else:
        use_mode = mode

    q1_k, q2_k = parameters

    with ProcessTimer() as time_ns:
        values = my_run(q1_k, q2_k, use_mode, inputs_file_beam, verbose)
    # print(f"optimizer (grads={grads}): {time_ns.elapsed_time}ns")
    error = values["error"]
    # print(q1_k, q2_k, error, values)

    history[optimizer].append([q1_k, q2_k])
    runtime[optimizer] = min(time_ns.elapsed_time, runtime[optimizer])

    if verbose:
        print(f"error={error}, q1_k={q1_k}, q2_k={q2_k}")

    if np.isnan(error):
        error = 1.0e99

    if grads:
        return (error, [values["derror_dq1_k"], values["derror_dq2_k"]])
    else:
        return error


def optimize_and_plot(opti, jac=None, plot=True):
    global grads, history, optimizer, runtime
    grads = jac == True  # noqa
    optimizer = opti
    history[optimizer] = []
    runtime[optimizer] = 9e99
    res = minimize(
        objective,
        initial_quad_strengths,
        method=optimizer,
        jac=jac,
        tol=1.0e-8,
        options=options,
        bounds=bounds,
    )
    ln = len(history[optimizer])
    print(optimizer, ln)
    q1_k, q2_k = zip(*history[optimizer])
    if plot:
        opt_str = optimizer
        if optimizer == "CG":
            opt_str = "Conjugate Gradient"
        plt.plot(q1_k, q2_k, label=f"{opt_str}: {ln}x")
    print("  Optimal parameters for k:", res.x)
    print("  L2 norm of alpha & beta at the optimum:", res.fun)
    print(f"  Min. runtime: {runtime[optimizer] / ln}ns")


if __name__ == "__main__":
    # Initial guess for the quadrople strengths
    initial_quad_strengths = np.array([-85, 80])

    # Bounds for values to test: (min, max)
    positive = (0, None)
    negative = (None, 0)
    bounds = [negative, positive]

    # optimizer specific values
    #   https://docs.scipy.org/doc/scipy/reference/optimize.minimize-neldermead.html
    #   https://docs.scipy.org/doc/scipy/reference/optimize.minimize-lbfgsb.html
    options = {
        "maxiter": 1000,
    }

    fig = plt.figure(figsize=(4.5, 2.2))

    # Call the optimizer
    optimize_and_plot("Nelder-Mead", False)
    optimize_and_plot("CG", "2-point", False)

    # gradient-based: CG, BFGS, Newton-CG, L-BFGS-B, TNC, SLSQP, dogleg, trust-ncg, trust-krylov, trust-exact and trust-constr
    optimize_and_plot("CG", True)

    # analytical result:
    #   q1_k = -103.12574100336
    #   q2_k = -q1_k

    fig.gca().annotate(
        "initial quadrupole strengths",  # Annotation text
        xy=(-85, 80),  # Point to annotate (arrow points here)
        xytext=(-95.8, 62),  # Text position (arrow starts here)
        arrowprops=dict(
            facecolor="black", shrink=0.1, width=1.5, headwidth=4, headlength=6
        ),
    )
    fig.gca().annotate(
        "solution",  # Annotation text
        xy=(-103.12574100336, 103.12574100336),  # Point to annotate (arrow points here)
        xytext=(-106, 113),  # Text position (arrow starts here)
        arrowprops=dict(
            facecolor="black", shrink=0.1, width=1.5, headwidth=4, headlength=6
        ),
    )

    plt.gca().set_xlim(-108, -80)
    plt.gca().set_ylim(60, 120)
    plt.xlabel(r"$k_1$  [1/m$^2$]")
    plt.ylabel(r"$k_2$  [1/m$^2$]")
    plt.legend(loc="upper right")
    plt.tight_layout()
    plt.show()
