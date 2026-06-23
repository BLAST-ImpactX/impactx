"""
This file is part of ImpactX

Copyright 2023 ImpactX contributors
Authors: Axel Huebl
License: BSD-3-Clause-LBNL
"""

from __future__ import annotations

__all__: list[str] = [
    "ix_beam_moments_history",
    "ix_pc_add_n_particles",
    "ix_pc_plot_mpl_phasespace",
    "register_ImpactXParticleContainer_extension",
]

def _as_real_device_vector(arr):
    """
    Convert an input into the AMReX real ``DeviceVector`` that
    ``AddNParticles`` expects, copying the data as needed.

    Accepts ``None`` (passed through), a NumPy or CuPy array (or array-like),
    or any pyAMReX ``PODVector``. A ``DeviceVector_real`` is passed through
    unchanged; any other allocator is copied via ``to_device`` (which picks
    the AMReX copy direction across memory spaces).
    """

def ix_beam_moments_history(self):
    """
    Return the history of the beam as calculated by the reduced beam characteristics on every step.
    """

def ix_pc_add_n_particles(
    self, x, y, t, px, py, pt, qm, bunch_charge=None, w=None, sx=None, sy=None, sz=None
):
    """
    Add new particles to the container for fixed s.

    The coordinate and weight arguments accept NumPy or CuPy arrays (or
    array-likes), as well as pyAMReX ``PODVector`` objects. Inputs are copied
    into device-compatible PODVectors as needed.

    Either the total charge (``bunch_charge``) or the weight of each particle
    (``w``) must be provided.

    Note: This can only be used *after* the grids have been created, i.e. after
    ``ImpactX.init_grids`` has been called.

    Parameters
    ----------
    x, y, t, px, py, pt : array_like
        Particle positions (x, y, time-of-flight c*t) and momenta.
    qm : float
        Charge over mass in 1/eV.
    bunch_charge : float, optional
        Total charge within a bunch in C.
    w : array_like, optional
        Weight of each particle: how many real particles to represent.
    sx, sy, sz : array_like, optional
        Spin components in x, y, z.
    """

def ix_pc_plot_mpl_phasespace(self, num_bins=50, root_rank=0):
    """
    Plot the longitudinal and transverse phase space projections with matplotlib.

    Parameters
    ----------
    self : ImpactXParticleContainer_*
        The particle container class in ImpactX
    num_bins : int, default=50
        The number of bins for spatial and momentum directions per plot axis.
    root_rank : int, default=0
        MPI root rank to reduce to in parallel runs.

    Returns
    -------
    A matplotlib figure with containing the plot.
    For MPI-parallel ranks, the figure is only created on the root_rank.
    """

def register_ImpactXParticleContainer_extension(ixpc):
    """
    ImpactXParticleContainer helper methods
    """
