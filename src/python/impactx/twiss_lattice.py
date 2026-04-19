#!/usr/bin/env python3
#
# Copyright 2026 ImpactX contributors
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-
"""
Linear-optics (Twiss) diagnostics for an ImpactX lattice.

This module computes the standard linear-optics observables along a
beam-line or ring:

    * beta and alpha functions (envelopes of the eigen-motion),
    * phase advances mu,
    * fractional tunes nu = mu / (2 pi) per one turn (for rings),
    * closed-orbit dispersion (D_x, D_px, D_y, D_py).

The inputs are the 6x6 linear transport maps M(s) accumulated along the
lattice, as provided by ``sim.lattice.map_trace(ref)`` from the C++ side.
The longitudinal coordinate is s (path length along the design orbit).

Parameterization
----------------
We use the coupled 6D parameterization of Wolski [Wolski2006]_. Concretely,
the 6x6 one-turn map M is diagonalized by a real symplectic matrix A,

    A^{-1} M A = block_diag(R_1, R_2, R_3),   R_k = [[cos 2 pi nu_k,  sin 2 pi nu_k],
                                                     [-sin 2 pi nu_k, cos 2 pi nu_k]],

where each mode k corresponds to a conjugate eigenvalue pair of M. The
column structure of A encodes the coupled Mais-Ripken beta functions
[Mais1982]_ [Willeke1988]_:
beta_{k,x} = A[0, 2k]^2 + A[0, 2k+1]^2 measures how much eigen-mode k
contributes to the horizontal position variance (and similarly for y, t).
In the uncoupled limit, these reduce to the classical Courant-Snyder
parameters [CourantSnyder1958]_.

Advantages over the classical Edwards-Teng parameterization
[EdwardsTeng1973]_:

    * remains well-defined and smooth through points where the
      Edwards-Teng block decoupling is singular (so-called "mode flip"
      near tune degeneracy),
    * extends uniformly to 6D, including synchro-betatron coupling from
      RF cavities,
    * produces the same beta_{k,j} functions used by, e.g., Xsuite and Bmad.

Phase-space ordering
--------------------
All 6x6 matrices in this module use the ordering

    z = (x, p_x, y, p_y, t, p_t)

with integer indices 0..5 respectively. Here x, y are transverse
positions in meters; p_x, p_y are the transverse canonical momenta
normalized by the reference momentum; t is arrival time of a particle
relative to the reference multiplied by the speed of light (meters);
p_t is the momentum deviation relative to the reference (dimensionless).

Canonical symplectic form
-------------------------
The canonical symplectic form for this ordering is the block-diagonal
matrix

    J = diag(omega, omega, omega),   omega = [[0, 1], [-1, 0]],

so that A is symplectic if A^T J A = J. Courant-Snyder / Twiss functions
are reconstructed from A via the Wolski formulas above.

Scope and limitations
---------------------
The current implementation does **not** solve the closed-orbit equation:
it assumes the reference particle passed to ``map_trace`` is already on
the closed orbit for rings (this is the common case for a nominally
matched design lattice). It also does not include synchrotron-radiation
damping, amplitude-dependent detuning, or other non-linear effects. For
those, a truncated-power-series (TPSA) or normal-form approach is needed,
which is out of scope here.

References
----------
.. [Wolski2006] A. Wolski, "Alternative approach to general coupled linear
    optics", Phys. Rev. ST Accel. Beams 9, 024001 (2006).
    doi:10.1103/PhysRevSTAB.9.024001.
.. [CourantSnyder1958] E. D. Courant and H. S. Snyder, "Theory of the
    Alternating-Gradient Synchrotron", Annals of Physics 3, 1-48 (1958).
    doi:10.1016/0003-4916(58)90012-5.
.. [Mais1982] H. Mais and G. Ripken, "Theory of coupled synchro-betatron
    oscillations", DESY Report M-82-05 (1982).
.. [Willeke1988] F. Willeke and G. Ripken, "Methods of beam optics",
    DESY Report 88-114 (1988).
.. [EdwardsTeng1973] D. A. Edwards and L. C. Teng, "Parametrization of
    Linear Coupled Motion in Periodic Systems", IEEE Trans. Nucl. Sci.
    20, 885-888 (1973). doi:10.1109/TNS.1973.4327279.
"""

from __future__ import annotations

import warnings
from collections.abc import Mapping
from typing import Any

import numpy as np


def _as_numpy_6x6(M: Any) -> np.ndarray:
    """Coerce a 6x6 matrix to a C-ordered ``float64`` numpy array.

    Accepts two input forms:

      * The pybind11-exposed ``Map6x6`` type from amrex/ImpactX (a thin
        wrapper around ``amrex::SmallMatrix<ParticleReal, 6, 6, F, 1>``).
        This is a Fortran-ordered, 1-indexed matrix type; calling its
        ``.to_numpy()`` method produces a standard 0-indexed, C-ordered
        copy. We go through that method explicitly because a plain
        ``np.asarray`` on ``Map6x6`` would surface the underlying
        Fortran layout (effectively transposing the matrix).
      * A numpy array. Returned as a contiguous ``float64`` copy.
    """
    to_numpy = getattr(M, "to_numpy", None)
    if callable(to_numpy):
        return np.ascontiguousarray(to_numpy(), dtype=np.float64)
    return np.ascontiguousarray(np.asarray(M), dtype=np.float64)


def _symplectic_form() -> np.ndarray:
    """Return the 6x6 canonical symplectic form J for ordering
    (x, p_x, y, p_y, t, p_t).

    J is block-diagonal with three 2x2 blocks ``[[0, 1], [-1, 0]]``. A
    matrix ``A`` is symplectic iff ``A^T J A = J``; the Poisson bracket
    of two phase-space functions f, g satisfies
    ``{f, g} = (grad f)^T J (grad g)``.
    """
    J = np.zeros((6, 6))
    for i in range(3):
        J[2 * i, 2 * i + 1] = 1.0
        J[2 * i + 1, 2 * i] = -1.0
    return J


_J = _symplectic_form()


def _dominant_subspace(v: np.ndarray) -> int:
    """Return the index (0, 1, or 2) of the canonical plane in which the
    length-6 vector ``v`` has the largest squared projection:

      * 0 for the (x, p_x) plane,
      * 1 for the (y, p_y) plane,
      * 2 for the (t, p_t) plane.

    Used to classify each eigen-mode of the one-turn map as predominantly
    horizontal, vertical, or longitudinal, so that the matched beta
    functions carry their usual physical labels.
    """
    projection = np.array(
        [
            abs(v[0]) ** 2 + abs(v[1]) ** 2,
            abs(v[2]) ** 2 + abs(v[3]) ** 2,
            abs(v[4]) ** 2 + abs(v[5]) ** 2,
        ]
    )
    return int(np.argmax(projection))


def wolski_a_from_m(M: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Compute the Wolski parameterization matrix ``A`` of a one-turn map.

    For a stable symplectic one-turn map ``M``, this returns the real
    6x6 matrix ``A`` such that

        A^T J A = J    (A is symplectic),
        A^{-1} M A  =  block_diag(R_1, R_2, R_3),

    where each ``R_k`` is a 2x2 rotation by the angle ``2 pi nu_k`` and
    ``nu_k`` is the fractional tune of mode k. The columns of A contain
    the real and imaginary parts of the properly normalized symplectic
    eigenvectors of M, and they furnish the coupled Mais-Ripken beta /
    alpha functions via :func:`twiss_from_a`.

    Algorithm
    ---------
      1. Diagonalize M numerically with ``numpy.linalg.eig``.
      2. Identify pairs of complex-conjugate eigenvalues; each such pair
         corresponds to one oscillatory (rotation) normal mode of the
         linearized motion. Purely real eigenvalues are grouped separately
         (see "parabolic modes" below).
      3. For each complex pair, pick the eigenvector with positive
         symplectic norm: the quantity
         ``sigma = -i * v^H J v`` is real for eigenvectors of a real
         symplectic matrix, and of the two conjugate partners one has
         ``sigma > 0``. The vector with ``sigma > 0`` is chosen, and
         then rescaled so that ``sigma = 1``. This is the Wolski
         normalization convention.
      4. Assign each normalized mode to the canonical plane (x, y, or t)
         in which its eigenvector has the largest projection, so that
         ``A[:, 0:2]``, ``A[:, 2:4]``, ``A[:, 4:6]`` correspond
         respectively to the horizontal, vertical, and longitudinal
         eigen-modes in the standard uncoupled limit.

    Parabolic modes
    ---------------
    A physical ring without radio-frequency (RF) cavities has a
    longitudinal one-turn block of the form ``[[1, c], [0, 1]]`` for some
    non-zero ``c`` (drift in arrival time). Both eigenvalues are equal to
    one, the block is defective (its Jordan form is *not* diagonal), and
    so no Wolski-type A-block with a finite beta function exists for this
    mode. This routine therefore substitutes the 2x2 identity block for it
    in A, assigns tune zero, and emits a ``RuntimeWarning``. Callers of
    :func:`compute_twiss` see the corresponding ``beta``, ``alpha`` and
    ``mu`` columns replaced by NaN along the whole lattice (see the
    "parabolic-mode masking" section of that function), because
    propagating the defective 2x2 shear would otherwise produce numbers
    that *look* like Twiss functions but have no physical meaning.
    Transverse optics (``beta_1x``, ``beta_2y``, etc.) and the closed-orbit
    dispersion remain valid in this case.

    Parameters
    ----------
    M : (6, 6) array_like
        Real symplectic 6x6 matrix. If M is not close enough to symplectic
        (for example, due to a lattice with gain / instability), the
        routine raises ``ValueError``.

    Returns
    -------
    A : (6, 6) numpy.ndarray
        Real Wolski parameterization matrix (symplectic, block-diagonalizes M).
    tunes : (3,) numpy.ndarray
        Fractional tunes ``nu_1, nu_2, nu_3`` of the three eigen-modes, in
        units of turns (i.e. ``mu_k = 2 pi nu_k``). The ordering corresponds
        to the modes identified with x, y, and t respectively.

    References
    ----------
    Primary reference for the 6D parameterization and the positive-norm
    normalization used here: [Wolski2006]_.
    Historical context of the coupled beta functions: [Mais1982]_,
    [Willeke1988]_. The uncoupled limit reduces to the classical
    Courant-Snyder parameters [CourantSnyder1958]_. A popular alternative
    that is *not* used here, because it is singular at mode-flip points,
    is the Edwards-Teng decoupling [EdwardsTeng1973]_.
    See the module-level docstring for the full bibliography.
    """
    M = _as_numpy_6x6(M)
    if M.shape != (6, 6):
        raise ValueError(f"Expected a 6x6 matrix, got shape {M.shape}.")

    eigvals, eigvecs = np.linalg.eig(M)

    # Identify complex-conjugate eigenvalue pairs (oscillatory modes). For a
    # real symplectic matrix, stable modes come as such pairs on the unit
    # circle.
    used = [False] * 6
    imag_tol = 1e-9
    complex_pairs: list[tuple[int, int]] = []
    for i in range(6):
        if used[i] or abs(eigvals[i].imag) < imag_tol:
            continue
        best_j = -1
        best_err = np.inf
        for j in range(6):
            if used[j] or j == i or abs(eigvals[j].imag) < imag_tol:
                continue
            err = abs(eigvals[i] - np.conj(eigvals[j]))
            if err < best_err:
                best_err = err
                best_j = j
        if best_j >= 0 and best_err < 1e-6 * max(1.0, abs(eigvals[i])):
            complex_pairs.append((i, best_j))
            used[i] = True
            used[best_j] = True

    # For each complex-conjugate pair, select the eigenvector with positive
    # symplectic norm and normalize it to sigma = -i v^H J v = 1. This is
    # the Wolski-convention normalization.
    normalized_modes: list[tuple[complex, np.ndarray]] = []
    for i, j in complex_pairs:
        v_i = eigvecs[:, i]
        sigma_i = (-1j * np.vdot(v_i, _J @ v_i)).real
        if sigma_i > 0:
            lam, v, sigma = eigvals[i], v_i, sigma_i
        else:
            lam, v, sigma = eigvals[j], eigvecs[:, j], -sigma_i
        if sigma <= 0.0:
            # Both partners have non-positive symplectic norm -- numerically
            # inconsistent with a symplectic matrix; skip this mode so the
            # final count check below surfaces the issue cleanly.
            continue
        v = v / np.sqrt(sigma)
        normalized_modes.append((lam, v))

    # Any eigenvalues not assigned to a complex-conjugate pair must be real.
    # For a symplectic matrix they come in pairs (lambda, 1/lambda); the
    # generic instance encountered here is the defective pair (1, 1) of a
    # no-RF longitudinal drift. We cannot build a Wolski A-block for such a
    # mode and substitute identity below.
    remaining_real = [i for i in range(6) if not used[i]]
    n_parabolic_pairs = len(remaining_real) // 2
    if n_parabolic_pairs + len(normalized_modes) != 3:
        raise ValueError(
            f"Expected 3 normal modes total; got {len(normalized_modes)} "
            f"complex-conjugate pairs and {n_parabolic_pairs} real pairs. "
            f"The input matrix may not be sufficiently symplectic."
        )

    if n_parabolic_pairs > 0:
        warnings.warn(
            f"{n_parabolic_pairs} of 3 normal modes of the one-turn map are "
            "not rotational (typically a longitudinal drift without an RF "
            "cavity, for which both eigenvalues equal one and the 2x2 Jordan "
            "block is defective). For these modes, the Wolski A-block is "
            "replaced by the 2x2 identity and the reported tune is zero; "
            "compute_twiss() returns NaN for the corresponding beta, alpha, "
            "and mu along the whole lattice. Transverse beta, alpha, mu, and "
            "dispersion remain valid.",
            RuntimeWarning,
            stacklevel=2,
        )

    # Assign rotational modes to slots (0 -> x, 1 -> y, 2 -> t) by the plane
    # of strongest eigenvector projection. If two modes would prefer the
    # same slot (possible under strong linear coupling), the one with the
    # larger projection wins the slot and the other falls into the first
    # free slot in order.
    slots: list[tuple[complex, np.ndarray] | None] = [None, None, None]
    ranked = sorted(
        normalized_modes,
        key=lambda item: max(
            abs(item[1][0]) ** 2 + abs(item[1][1]) ** 2,
            abs(item[1][2]) ** 2 + abs(item[1][3]) ** 2,
            abs(item[1][4]) ** 2 + abs(item[1][5]) ** 2,
        ),
        reverse=True,
    )
    for lam, v in ranked:
        preferred = _dominant_subspace(v)
        if slots[preferred] is None:
            slots[preferred] = (lam, v)
        else:
            for kk in range(3):
                if slots[kk] is None:
                    slots[kk] = (lam, v)
                    break

    A = np.zeros((6, 6))
    tunes = np.zeros(3)
    for k, entry in enumerate(slots):
        if entry is None:
            # Parabolic (non-rotational) mode: identity A-block, tune = 0.
            A[2 * k, 2 * k] = 1.0
            A[2 * k + 1, 2 * k + 1] = 1.0
        else:
            lam, v = entry
            tunes[k] = np.angle(lam) / (2.0 * np.pi)
            # Column pair for mode k is (sqrt(2) Re v, sqrt(2) Im v); this
            # factor makes A^T J A = J given the v^H (iJ) v = 1 normalization.
            A[:, 2 * k] = np.sqrt(2.0) * v.real
            A[:, 2 * k + 1] = np.sqrt(2.0) * v.imag

    return A, tunes


def a_from_courant_snyder(
    beta_x: float,
    alpha_x: float,
    beta_y: float,
    alpha_y: float,
    beta_t: float = 1.0,
    alpha_t: float = 0.0,
) -> np.ndarray:
    """Construct a Wolski A matrix from uncoupled Courant-Snyder parameters.

    This is intended for transfer-line and linac use cases, where the user
    specifies initial Twiss parameters at the entrance and the Twiss
    functions are then transported along the lattice.

    For a fully uncoupled lattice, each 2x2 block of A reduces to the
    classical Courant-Snyder form [CourantSnyder1958]_

        A_k = [[ sqrt(beta_k),          0            ],
               [ -alpha_k / sqrt(beta_k), 1 / sqrt(beta_k) ]],

    for k in {x, y, t}, which satisfies ``A_k^T omega A_k = omega``.
    The 6x6 matrix is the block-diagonal assembly of the three.

    Parameters
    ----------
    beta_x, beta_y, beta_t : float
        Beta functions (meters). Must be strictly positive. ``beta_t``
        defaults to 1 m because the longitudinal Twiss is often irrelevant
        for transfer-line analysis and a non-zero value is required only to
        keep the 6x6 matrix non-singular.
    alpha_x, alpha_y, alpha_t : float
        Alpha functions (dimensionless). Default to zero.

    Returns
    -------
    A : (6, 6) numpy.ndarray
        Block-diagonal Wolski parameterization matrix.
    """
    A = np.zeros((6, 6))
    for k, (beta, alpha) in enumerate(
        [(beta_x, alpha_x), (beta_y, alpha_y), (beta_t, alpha_t)]
    ):
        if beta <= 0.0:
            raise ValueError(f"beta must be strictly positive; got {beta}.")
        sq_beta = np.sqrt(beta)
        row = 2 * k
        A[row, 2 * k] = sq_beta
        A[row, 2 * k + 1] = 0.0
        A[row + 1, 2 * k] = -alpha / sq_beta
        A[row + 1, 2 * k + 1] = 1.0 / sq_beta
    return A


def twiss_from_a(A: np.ndarray) -> dict[str, float]:
    """Extract the coupled Mais-Ripken Twiss parameters from ``A`` at one s.

    Given the Wolski matrix ``A(s)`` at some longitudinal position, this
    returns a flat dictionary of the nine coupled beta functions and nine
    coupled alpha functions [Mais1982]_ [Willeke1988]_:

        beta_{k, j}(s)  = A[2j, 2k-2]^2 + A[2j, 2k-1]^2,
        alpha_{k, j}(s) = -(A[2j, 2k-2] * A[2j+1, 2k-2]
                            + A[2j, 2k-1] * A[2j+1, 2k-1]),

    for eigen-mode index k in {1, 2, 3} and canonical plane index j mapped
    through {x -> 0, y -> 1, t -> 2}. The diagonal entries
    (beta_{1,x}, beta_{2,y}, beta_{3,t}) are the traditional (uncoupled)
    Twiss beta functions. The off-diagonal entries (beta_{1,y},
    beta_{2,x}, ...) vanish identically in an uncoupled lattice and
    measure the leakage of mode k into the other canonical planes when
    coupling is present.

    Parameters
    ----------
    A : (6, 6) array_like
        Wolski parameterization matrix at a given s.

    Returns
    -------
    twiss : dict[str, float]
        Keys are of the form ``beta_<k><plane>`` and ``alpha_<k><plane>``
        with ``<k>`` in ``{1, 2, 3}`` and ``<plane>`` in ``{x, y, t}``.
    """
    A = np.asarray(A)
    out: dict[str, float] = {}
    planes = [(0, "x"), (2, "y"), (4, "t")]
    for mode in (1, 2, 3):
        col_real = 2 * (mode - 1)
        col_imag = col_real + 1
        for pos_row, plane in planes:
            mom_row = pos_row + 1
            beta = float(A[pos_row, col_real] ** 2 + A[pos_row, col_imag] ** 2)
            alpha = -float(
                A[pos_row, col_real] * A[mom_row, col_real]
                + A[pos_row, col_imag] * A[mom_row, col_imag]
            )
            out[f"beta_{mode}{plane}"] = beta
            out[f"alpha_{mode}{plane}"] = alpha
    return out


def tunes_from_m(M: np.ndarray) -> np.ndarray:
    """Return the three fractional tunes of a one-turn map.

    The fractional tune of mode k is ``nu_k = arg(lambda_k) / (2 pi)``,
    where ``lambda_k`` is the eigenvalue of the one-turn map associated
    with the positive-norm eigenvector of that mode. Fractional tunes are
    defined modulo one; the value returned here lies in ``[-1/2, +1/2]``.

    The ordering corresponds to the assignment of eigen-modes to the
    horizontal, vertical, and longitudinal canonical planes performed by
    :func:`wolski_a_from_m`.

    Parameters
    ----------
    M : (6, 6) array_like
        Real symplectic one-turn matrix.

    Returns
    -------
    tunes : (3,) numpy.ndarray
        Fractional tunes ``(nu_1, nu_2, nu_3)``.
    """
    _A, tunes = wolski_a_from_m(M)
    return tunes


def dispersion_from_m(M: np.ndarray) -> np.ndarray:
    """Closed-orbit dispersion of a periodic ring at the starting s.

    The closed-orbit dispersion is the change in the periodic closed orbit
    induced by a small momentum offset delta = p_t. It satisfies the
    first-order periodic equation

        (I - M_4x4) D  =  M[0:4, 5],

    where ``M_4x4`` is the transverse 4x4 block of the one-turn matrix,
    ``M[0:4, 5]`` is the transverse-to-longitudinal coupling column, and
    ``D = (D_x, D_{p_x}, D_y, D_{p_y})`` is the periodic dispersion four-vector.
    This assumes the ring is "longitudinally matched", i.e. a particle
    with constant p_t returns to its own starting phase after one turn;
    this holds for rings without synchrotron-radiation damping and
    without net energy gain over a turn.

    Parameters
    ----------
    M : (6, 6) array_like
        One-turn matrix.

    Returns
    -------
    D : (4,) numpy.ndarray
        Periodic dispersion ``(D_x, D_{p_x}, D_y, D_{p_y})``.

    Raises
    ------
    numpy.linalg.LinAlgError
        If ``(I - M_4x4)`` is singular, e.g. for an integer transverse
        tune or an unstable lattice. The closed-orbit dispersion is not
        well-defined in that case.
    """
    M = _as_numpy_6x6(M)
    A_mat = np.eye(4) - M[0:4, 0:4]
    rhs = M[0:4, 5]
    try:
        return np.linalg.solve(A_mat, rhs)
    except np.linalg.LinAlgError as exc:
        raise np.linalg.LinAlgError(
            "Closed-orbit dispersion requires (I - M_4x4) to be non-singular. "
            "The current lattice has an integer transverse tune or is "
            "unstable, so the periodic dispersion is not well-defined."
        ) from exc


def propagate_dispersion(D0: np.ndarray, M_cum: np.ndarray) -> np.ndarray:
    """Transport the dispersion four-vector from s=0 to position s.

    Given the cumulative linear transport matrix ``M_cum`` from the start
    of the lattice to some later position s, the dispersion four-vector
    ``D(s) = (D_x, D_{p_x}, D_y, D_{p_y})`` satisfies

        D(s) = M_cum[0:4, 0:4] * D(0) + M_cum[0:4, 5].

    The first term transports the incoming dispersion through the 4x4
    transverse block; the second term is the dispersion *generated* by the
    transverse-longitudinal coupling column (non-zero in dipoles and other
    elements with an energy-dependent effective length).

    Parameters
    ----------
    D0 : (4,) array_like
        Dispersion at s = 0.
    M_cum : (6, 6) array_like
        Cumulative transport map from s = 0 to the position of interest.

    Returns
    -------
    D : (4,) numpy.ndarray
        Dispersion at the position described by ``M_cum``.
    """
    M_cum = _as_numpy_6x6(M_cum)
    D0 = np.asarray(D0, dtype=np.float64)
    return M_cum[0:4, 0:4] @ D0 + M_cum[0:4, 5]


def _unwrap_phase(prev: float, curr: float) -> float:
    """Continuously unwrap a phase angle against the previous sample.

    ``numpy.arctan2`` returns values in the principal interval (-pi, pi],
    but a physically continuous phase advance mu(s) may grow far beyond
    that range along a lattice. This helper adds or subtracts integer
    multiples of 2 pi from ``curr`` so that it lies within pi of
    ``prev``, producing a continuous sequence when applied step by step
    along the lattice. It works in both signs: for a mode with positive
    tune mu increases with s, while for negative tune it decreases.

    Parameters
    ----------
    prev : float
        Previously unwrapped phase sample (radians).
    curr : float
        New phase sample from ``arctan2`` (radians).

    Returns
    -------
    float
        ``curr`` adjusted by an integer multiple of 2 pi.
    """
    two_pi = 2.0 * np.pi
    while curr < prev - np.pi:
        curr += two_pi
    while curr > prev + np.pi:
        curr -= two_pi
    return curr


def compute_twiss(
    sim: Any,
    init: Mapping[str, float] | None = None,
) -> dict[str, np.ndarray]:
    """Compute the Twiss (linear-optics) table along a lattice.

    Two lattice categories are supported:

    * **Periodic ring** (``init=None``). The one-turn map
      ``M`` is obtained via ``sim.map_trace()``, its Wolski eigendecomposition
      provides the matched initial ``A_0``, and the periodic closed-orbit
      dispersion is solved from ``(I - M_4x4) D_0 = M[0:4, 5]``. Fractional
      tunes are extracted from the eigenvalue phases. The reference particle
      stored in ``sim.beam.ref`` is assumed to be on the closed orbit.

      *Parabolic-mode masking.* If ``wolski_a_from_m`` reports one or more
      modes as parabolic (tune identically zero, identity A-block; typically
      a longitudinal drift in a ring without RF), the corresponding
      ``beta_k*``, ``alpha_k*`` and ``mu_k`` columns are returned as NaN
      along the whole lattice instead of the meaningless shear-propagated
      values that a naive ``M(s) @ A_0`` would produce. Dispersion and
      transverse Twiss remain valid.
    * **Transfer line / linac** (``init`` given). The user provides initial
      uncoupled Courant-Snyder parameters (``beta_x``, ``alpha_x``,
      ``beta_y``, ``alpha_y``, optionally ``beta_t`` and ``alpha_t``) and an
      initial dispersion (``disp_x``, ``disp_px``, ``disp_y``, ``disp_py``,
      all defaulting to zero). Twiss is then transported from the entrance
      through the lattice; no matching is enforced.

    The computation walks the cumulative linear map ``M(s)`` returned by
    ``sim.map_trace()`` and evaluates

        A(s)  = M(s) A_0,
        D(s)  = M(s)[0:4, 0:4] D_0 + M(s)[0:4, 5],

    then extracts beta, alpha from ``A(s)`` via :func:`twiss_from_a` and
    phase advance ``mu_k(s)`` as the continuous arctan2 of the first-row
    columns of A(s).

    Parameters
    ----------
    sim : impactx.ImpactX
        Simulation object with ``sim.lattice`` populated and
        ``sim.beam.ref`` set to the reference orbit at the starting s.
    init : Mapping or None
        For line mode, a dict-like containing ``beta_x``, ``beta_y``,
        and optionally ``alpha_x``, ``alpha_y``, ``beta_t``, ``alpha_t``,
        ``disp_x``, ``disp_px``, ``disp_y``, ``disp_py``. If ``None``,
        the lattice is treated as a periodic ring and Twiss is matched
        to the one-turn map.

    Returns
    -------
    dict of numpy.ndarray
        One entry per element boundary of the lattice (plus an initial
        entry at the start). Keys:

          * ``s`` -- integrated path length in meters (shape ``(n,)``).
          * ``name``, ``type`` -- element identifiers (shape ``(n,)``,
            dtype object).
          * ``beta_1x``, ``beta_1y``, ..., ``beta_3t`` -- coupled
            Mais-Ripken beta functions, units of meters. NaN for a
            parabolic mode in ring mode (see above).
          * ``alpha_1x``, ..., ``alpha_3t`` -- coupled alpha functions,
            dimensionless. NaN for a parabolic mode in ring mode.
          * ``mu_1``, ``mu_2``, ``mu_3`` -- cumulative phase advance per
            mode, in radians, continuously unwrapped and referenced to
            zero at s = 0. NaN for a parabolic mode in ring mode.
          * ``disp_x``, ``disp_px``, ``disp_y``, ``disp_py`` --
            dispersion four-vector along the lattice.
          * ``beta_x = beta_1x`` etc. -- convenience aliases matching the
            uncoupled MAD-X naming convention.
          * For periodic rings only, ``tunes`` -- a length-3 array of the
            fractional tunes ``(nu_1, nu_2, nu_3)``. Note that this entry
            has shape ``(3,)``, unlike the per-position arrays.
    """
    # Route through sim.lattice.map_trace(ref) so the reference particle is
    # copied by pybind11 on the way in (see KnownElementsList.map_trace).
    # This preserves the caller's ``sim.beam.ref`` in case they decide to run
    # a particle / envelope tracking simulation afterward.
    trace = sim.lattice.map_trace(sim.beam.ref)
    if not trace:
        raise ValueError(
            "Lattice is empty -- sim.lattice.map_trace() returned no entries."
        )

    # The end-to-end (one-turn, for a ring) map is the last trace entry.
    M_oneturn = _as_numpy_6x6(trace[-1]["M"])

    # Modes whose per-position beta / alpha / mu should be masked to NaN
    # after propagation. Populated only in ring (periodic) mode, when
    # wolski_a_from_m reports a parabolic mode (tune exactly zero, A-block
    # replaced by the identity). Masking is not applied in line mode
    # because there the user-supplied A_0 is by construction a valid
    # Courant-Snyder envelope and its propagation is physically meaningful
    # (even the longitudinal envelope growing through a drift).
    parabolic_modes: set[int] = set()

    if init is None:
        # Periodic ring: Wolski-match to the one-turn map.
        A0, tunes = wolski_a_from_m(M_oneturn)
        D0 = dispersion_from_m(M_oneturn)
        is_periodic = True
        parabolic_modes = {k for k, nu in enumerate(tunes) if nu == 0.0}
    else:
        # Transfer line / linac: use user-supplied uncoupled Courant-Snyder.
        A0 = a_from_courant_snyder(
            beta_x=float(init["beta_x"]),
            alpha_x=float(init.get("alpha_x", 0.0)),
            beta_y=float(init["beta_y"]),
            alpha_y=float(init.get("alpha_y", 0.0)),
            beta_t=float(init.get("beta_t", 1.0)),
            alpha_t=float(init.get("alpha_t", 0.0)),
        )
        D0 = np.array(
            [
                float(init.get("disp_x", 0.0)),
                float(init.get("disp_px", 0.0)),
                float(init.get("disp_y", 0.0)),
                float(init.get("disp_py", 0.0)),
            ]
        )
        tunes = None
        is_periodic = False

    # Sweep along the trace, transporting A and D, and accumulating the
    # (continuously unwrapped) phase advance of each mode.
    n = len(trace)
    s_arr = np.zeros(n)
    names: list[str] = []
    types: list[str] = []
    twiss_rows: list[dict[str, float]] = []
    disp_rows = np.zeros((n, 4))
    mu = np.zeros((n, 3))

    prev_mu = np.zeros(3)
    for i, entry in enumerate(trace):
        s_arr[i] = entry["s"]
        names.append(entry["name"])
        types.append(entry["type"])

        M_cum = _as_numpy_6x6(entry["M"])
        A_here = M_cum @ A0

        twiss_rows.append(twiss_from_a(A_here))
        disp_rows[i] = propagate_dispersion(D0, M_cum)

        # Phase advance of mode k: take the dominant-plane first row of A,
        # which rotates by angle mu_k as A is transported, and extract its
        # instantaneous phase via arctan2. The result is unwrapped against
        # the previous sample so the sequence is continuous along s.
        for mode_k in range(3):
            pos_row = 2 * mode_k
            col_real = 2 * mode_k
            col_imag = col_real + 1
            phi = float(
                np.arctan2(A_here[pos_row, col_imag], A_here[pos_row, col_real])
            )
            if i == 0:
                mu[0, mode_k] = phi
                prev_mu[mode_k] = phi
            else:
                phi_unwrapped = _unwrap_phase(prev_mu[mode_k], phi)
                mu[i, mode_k] = phi_unwrapped
                prev_mu[mode_k] = phi_unwrapped

    # Shift mu so that it starts exactly at zero at s = 0 (the raw arctan2
    # of A_0 depends on an arbitrary overall phase of the eigenvector
    # returned by numpy.linalg.eig).
    mu = mu - mu[0:1, :]

    out: dict[str, np.ndarray] = {
        "s": s_arr,
        "name": np.array(names, dtype=object),
        "type": np.array(types, dtype=object),
        "mu_1": mu[:, 0],
        "mu_2": mu[:, 1],
        "mu_3": mu[:, 2],
        "disp_x": disp_rows[:, 0],
        "disp_px": disp_rows[:, 1],
        "disp_y": disp_rows[:, 2],
        "disp_py": disp_rows[:, 3],
    }
    # One array per beta / alpha key, indexed the same as s_arr above.
    for key in twiss_rows[0].keys():
        out[key] = np.array([row[key] for row in twiss_rows])

    if is_periodic:
        # Single 3-vector of fractional tunes, not a per-position array.
        out["tunes"] = tunes

    # Parabolic-mode masking. In ring mode, any mode flagged as
    # non-rotational by wolski_a_from_m (tune == 0, A-block replaced by
    # identity) has no physically meaningful beta / alpha / mu: propagating
    # the identity A-block through a cumulative 2x2 shear produces
    # monotonically growing numbers that a naive reader could mistake for a
    # Courant-Snyder envelope. Replace those columns with NaN so downstream
    # plots and tables signal "undefined" unambiguously.
    for mode_k in parabolic_modes:
        mode_label = mode_k + 1  # keys use 1-based mode indexing
        nan_column = np.full(n, np.nan)
        for plane in ("x", "y", "t"):
            out[f"beta_{mode_label}{plane}"] = nan_column.copy()
            out[f"alpha_{mode_label}{plane}"] = nan_column.copy()
        out[f"mu_{mode_label}"] = nan_column.copy()

    # Convenience aliases matching the uncoupled MAD-X naming. In an
    # uncoupled lattice these are identical to the diagonal Mais-Ripken
    # entries; in a coupled lattice they report only the diagonal
    # contribution and the off-diagonal coupling shows up in beta_1y etc.
    # After parabolic-mode masking, e.g. beta_t picks up the NaN masking
    # on beta_3t for a no-RF ring.
    out["beta_x"] = out["beta_1x"]
    out["alpha_x"] = out["alpha_1x"]
    out["beta_y"] = out["beta_2y"]
    out["alpha_y"] = out["alpha_2y"]
    out["beta_t"] = out["beta_3t"]
    out["alpha_t"] = out["alpha_3t"]

    return out
