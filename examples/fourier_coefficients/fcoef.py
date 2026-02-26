import math

import numpy as np


def read_data(filename):
    """Read data from file."""
    data = np.loadtxt(filename)
    return data


def calculate_coefficients(data, ncoef):
    """Calculate Fourier coefficients of on-axis field data.

    Uses the trapezoidal rule with linear interpolation to compute
    cosine and sine Fourier coefficients of the field profile,
    centered about the midpoint of the data range.

    Parameters
    ----------
    data : numpy.ndarray
        Two-column array with shape (N, 2).  Column 0 is the
        longitudinal position z in meters.  Column 1 is the field or
        field gradient in arbitrary units (typically scaled so that the
        peak absolute value is 1; physical units are set later via
        element parameters such as ``bscale`` or ``escale``).
    ncoef : int
        Number of Fourier coefficients to compute.

    Returns
    -------
    cos_coeffs : numpy.ndarray
        Cosine Fourier coefficients, length *ncoef*.
    sin_coeffs : numpy.ndarray
        Sine Fourier coefficients, length *ncoef*.
    """
    zdata = data[:, 0]
    edata = data[:, 1]
    ndatareal = len(zdata)

    zlen = zdata[-1] - zdata[0]
    zmid = (zdata[-1] + zdata[0]) / 2
    zhalf = zlen / 2.0
    h = zlen / (ndatareal - 1)

    j = np.arange(ncoef)  # (ncoef,)

    # Endpoint correction (trapezoidal rule)
    zz0 = zdata[0] - zmid
    zz1 = zdata[-1] - zmid
    angle0 = j * 2 * np.pi * zz0 / zlen  # (ncoef,)
    angle1 = j * 2 * np.pi * zz1 / zlen  # (ncoef,)

    cos_coeffs = (-0.5 * edata[0] * np.cos(angle0) * h) / zhalf
    sin_coeffs = (-0.5 * edata[0] * np.sin(angle0) * h) / zhalf
    cos_coeffs -= (0.5 * edata[-1] * np.cos(angle1) * h) / zhalf
    sin_coeffs -= (0.5 * edata[-1] * np.sin(angle1) * h) / zhalf

    # Interior points: interpolate field onto uniform grid, then integrate
    zz_uniform = np.arange(ndatareal) * h + zdata[0]  # (ndatareal,)
    ez1 = np.interp(zz_uniform, zdata, edata)  # (ndatareal,)

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


def write_data(cos_coeffs, sin_coeffs, zdata, filename):
    """Write data to file."""
    zlen = zdata[-1] - zdata[0]
    zmid = (zdata[-1] + zdata[0]) / 2

    with open("onaxis_data.out", "w") as f:
        for i, z in enumerate(zdata):
            zz = z - zmid
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
            f.write(f"{z} {tmpsum} {tmpsump} {tmpsumpp}\n")


def main():
    data = read_data("onaxis_data.in")
    ncoef = int(input("How many Fourier coefficients do you want? "))
    cos_coeffs, sin_coeffs = calculate_coefficients(data, ncoef)
    write_coefficients(cos_coeffs, sin_coeffs, "fcoef.out")
    write_data(cos_coeffs, sin_coeffs, data[:, 0], "onaxis_datax")


if __name__ == "__main__":
    main()
