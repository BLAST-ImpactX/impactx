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
    emax : float
        Peak absolute field value from the input data (arbitrary units).
    """
    zdata = data[:, 0]
    edata = data[:, 1]
    emax = np.max(np.abs(edata))
    ndatareal = len(zdata)

    zlen = zdata[-1] - zdata[0]
    zmid = (zdata[-1] + zdata[0]) / 2
    zhalf = zlen / 2.0
    h = zlen / (ndatareal - 1)
    pi = math.pi

    cos_coeffs = np.zeros(ncoef)
    sin_coeffs = np.zeros(ncoef)

    """Contribution from the endpoints."""
    for j in range(ncoef):
        zz = zdata[0] - zmid
        cos_coeffs[j] += (
            -0.5 * edata[0] * math.cos((j) * 2 * pi * zz / zlen) * h
        ) / zhalf
        sin_coeffs[j] += (
            -0.5 * edata[0] * math.sin((j) * 2 * pi * zz / zlen) * h
        ) / zhalf
        zz = zdata[-1] - zmid
        cos_coeffs[j] -= (
            0.5 * edata[-1] * math.cos((j) * 2 * pi * zz / zlen) * h
        ) / zhalf
        sin_coeffs[j] -= (
            0.5 * edata[-1] * math.sin((j) * 2 * pi * zz / zlen) * h
        ) / zhalf

    """Contribution from the interior points."""
    for i in range(ndatareal):
        zz = i * h + zdata[0]
        klo = 0
        khi = ndatareal - 1
        while khi - klo > 1:
            k = (khi + klo) // 2
            if zdata[k] > zz:
                khi = k
            else:
                klo = k
        hstep = zdata[khi] - zdata[klo]
        slope = (edata[khi] - edata[klo]) / hstep
        ez1 = edata[klo] + slope * (zz - zdata[klo])

        zz = zdata[0] + i * h - zmid
        for j in range(ncoef):
            cos_coeffs[j] += (ez1 * math.cos((j) * 2 * pi * zz / zlen) * h) / zhalf
            sin_coeffs[j] += (ez1 * math.sin((j) * 2 * pi * zz / zlen) * h) / zhalf

    return cos_coeffs, sin_coeffs, emax


def write_coefficients(cos_coeffs, sin_coeffs, filename):
    """Write coefficients to file."""
    with open(filename, "w") as f:
        f.write("j  cos_coef[j]  sin_coef[j] \n")
        for j, (coef, coef2) in enumerate(zip(cos_coeffs, sin_coeffs)):
            f.write(f"{j} {coef} {coef2}\n")


def write_data(cos_coeffs, sin_coeffs, zdata, emax, filename):
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
    cos_coeffs, sin_coeffs, emax = calculate_coefficients(data, ncoef)
    write_coefficients(cos_coeffs, sin_coeffs, "fcoef.out")
    write_data(cos_coeffs, sin_coeffs, data[:, 0], emax, "onaxis_datax")


if __name__ == "__main__":
    main()
