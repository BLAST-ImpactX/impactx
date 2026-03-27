import numpy as np
from scipy.optimize import minimize_scalar

from impactx import push


def objective(parameter, ref, element):
    """
    A function that is evaluated by the optimizer.

    Parameters
    ----------
    parameter:
      rf cavity phase

    Returns
    -------
    Negative of the RF cavity energy gain in MeV (to minimize).
    """

    # adjust the RF cavity phase
    phase_opt = parameter
    element.phase = phase_opt

    # store the incoming energy and copy the reference particle
    KE_in = ref.kin_energy_MeV
    ref_copy = ref.copy()

    # push the copy of the reference particle
    push(ref_copy, element)
    KE_fin = ref_copy.kin_energy_MeV

    # evaluate the objective
    loss = KE_in - KE_fin

    if np.isnan(loss):
        loss = 1.0e99

    return loss


def optimize(ref, element):

    # optimizer specific options
    options = {"maxiter": 2000, "disp": 1}

    # Call the optimizer
    res = minimize_scalar(
        objective,
        method="bounded",
        args=(ref, element),
        tol=1.0e-8,
        options=options,
        bounds=(-180, 180),
    )

    # Optimization result
    phase_opt = res.x
    e_gain = -1.0 * res.fun

    return phase_opt, e_gain
