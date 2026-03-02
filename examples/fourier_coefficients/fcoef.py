import math

import numpy as np


def read_data(filename):
    """Read data from file."""
    data = np.loadtxt(filename)
    return data[:, 0], data[:, 1]


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
    from impactx import fourier_coefficients

    z, field_or_gradient = read_data("onaxis_data.in")
    ncoef = int(input("How many Fourier coefficients do you want? "))
    cos_coeffs, sin_coeffs = fourier_coefficients(z, field_or_gradient, ncoef)
    write_coefficients(cos_coeffs, sin_coeffs, "fcoef.out")
    write_data(cos_coeffs, sin_coeffs, z, "onaxis_datax")


if __name__ == "__main__":
    main()
