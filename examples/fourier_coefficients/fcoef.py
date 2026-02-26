import math

import numpy as np


def read_data(filename):
    """Read data from file."""
    data = np.loadtxt(filename)
    return data[:, 0], data[:, 1]


def calculate_coefficients(z, field_or_gradient, ncoef):
    """Calculate Fourier coefficients of on-axis field data.

    Uses the trapezoidal rule with linear interpolation to compute
    cosine and sine Fourier coefficients of the field profile,
    centered about the midpoint of the data range.

    Parameters
    ----------
    z : numpy.ndarray
        Longitudinal positions in meters, shape (N,).
    field_or_gradient : numpy.ndarray
        On-axis field or field gradient values in arbitrary units,
        shape (N,).  Typically scaled so that the peak absolute value
        is 1; physical units are set later via element parameters such
        as ``bscale`` or ``escale``.
    ncoef : int
        Number of Fourier coefficients to compute.

    Returns
    -------
    cos_coeffs : numpy.ndarray
        Cosine Fourier coefficients, length *ncoef*.
    sin_coeffs : numpy.ndarray
        Sine Fourier coefficients, length *ncoef*.
    """
    ndatareal = len(z)

    zlen = z[-1] - z[0]
    zmid = (z[-1] + z[0]) / 2
    zhalf = zlen / 2.0
    h = zlen / (ndatareal - 1)

    j = np.arange(ncoef)  # (ncoef,)

    # Endpoint correction (trapezoidal rule)
    zz0 = z[0] - zmid
    zz1 = z[-1] - zmid
    angle0 = j * 2 * np.pi * zz0 / zlen  # (ncoef,)
    angle1 = j * 2 * np.pi * zz1 / zlen  # (ncoef,)

    cos_coeffs = (-0.5 * field_or_gradient[0] * np.cos(angle0) * h) / zhalf
    sin_coeffs = (-0.5 * field_or_gradient[0] * np.sin(angle0) * h) / zhalf
    cos_coeffs -= (0.5 * field_or_gradient[-1] * np.cos(angle1) * h) / zhalf
    sin_coeffs -= (0.5 * field_or_gradient[-1] * np.sin(angle1) * h) / zhalf

    # Interior points: interpolate field onto uniform grid, then integrate
    zz_uniform = np.arange(ndatareal) * h + z[0]  # (ndatareal,)
    ez1 = np.interp(zz_uniform, z, field_or_gradient)  # (ndatareal,)

    zz_centered = zz_uniform - zmid  # (ndatareal,)
    # Outer product: angles[i, j] = j * 2 * pi * zz_centered[i] / zlen
    angles = np.outer(zz_centered, j * 2 * np.pi / zlen)  # (ndatareal, ncoef)

    cos_coeffs += np.sum(ez1[:, np.newaxis] * np.cos(angles) * h, axis=0) / zhalf
    sin_coeffs += np.sum(ez1[:, np.newaxis] * np.sin(angles) * h, axis=0) / zhalf

    return cos_coeffs, sin_coeffs


def write_coefficients(cos_coeffs, sin_coeffs, filename):
    """Write coefficients to file."""
    with open(filename, "w") as f:
        f.write("j  cos_coef[j]  sin_coef[j] \n")
        for j, (coef, coef2) in enumerate(zip(cos_coeffs, sin_coeffs)):
            f.write(f"{j} {coef} {coef2}\n")


def write_data(cos_coeffs, sin_coeffs, z, filename):
    """Write data to file."""
    zlen = z[-1] - z[0]
    zmid = (z[-1] + z[0]) / 2

    with open("onaxis_data.out", "w") as f:
        for i, zi in enumerate(z):
            zz = zi - zmid
            tmpsum = 0.5 * cos_coeffs[0]
            tmpsump = 0.0
            tmpsumpp = 0.0
            for j in range(1, len(cos_coeffs)):
                tmpsum += cos_coeffs[j] * math.cos(
                    (j) * 2 * math.pi * zz / zlen
                ) + sin_coeffs[j] * math.sin((j) * 2 * math.pi * zz / zlen)
                tmpsump += (
                    -(j)
                    * 2
                    * math.pi
                    * cos_coeffs[j]
                    * math.sin((j) * 2 * math.pi * zz / zlen)
                    / zlen
                    + (j)
                    * 2
                    * math.pi
                    * sin_coeffs[j]
                    * math.cos((j) * 2 * math.pi * zz / zlen)
                    / zlen
                )
                tmpsumpp += -(((j) * 2 * math.pi / zlen) ** 2) * cos_coeffs[
                    j
                ] * math.cos((j) * 2 * math.pi * zz / zlen) - (
                    (j) * 2 * math.pi / zlen
                ) ** 2 * sin_coeffs[j] * math.sin((j) * 2 * math.pi * zz / zlen)
            f.write(f"{zi} {tmpsum} {tmpsump} {tmpsumpp}\n")


def main():
    z, field_or_gradient = read_data("onaxis_data.in")
    ncoef = int(input("How many Fourier coefficients do you want? "))
    cos_coeffs, sin_coeffs = calculate_coefficients(z, field_or_gradient, ncoef)
    write_coefficients(cos_coeffs, sin_coeffs, "fcoef.out")
    write_data(cos_coeffs, sin_coeffs, z, "onaxis_datax")


if __name__ == "__main__":
    main()
