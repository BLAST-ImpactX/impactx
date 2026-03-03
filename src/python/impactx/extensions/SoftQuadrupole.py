"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Chad Mitchell, Axel Huebl
License: BSD-3-Clause-LBNL
"""


def register_SoftQuadrupole_extension(cls):
    """Extend SoftQuadrupole with an alternative constructor that accepts
    raw on-axis field/gradient data and computes Fourier coefficients
    internally.

    Parameters
    ----------
    cls : type
        The pybind11 ``SoftQuadrupole`` class to extend.
    """
    _original_init = cls.__init__

    def _new_init(
        self,
        *,
        ds,
        gscale,
        cos_coefficients=None,
        sin_coefficients=None,
        z=None,
        field_or_gradient=None,
        ncoef=None,
        dx=0,
        dy=0,
        rotation=0,
        aperture_x=0,
        aperture_y=0,
        mapsteps=1,
        nslice=1,
        name=None,
    ):
        has_coefficients = cos_coefficients is not None or sin_coefficients is not None
        has_field_data = (
            z is not None or field_or_gradient is not None or ncoef is not None
        )

        if has_coefficients and has_field_data:
            raise ValueError(
                "SoftQuadrupole: provide either (cos_coefficients, sin_coefficients) "
                "or (z, field_or_gradient, ncoef), not both."
            )

        if has_field_data:
            if z is None or field_or_gradient is None or ncoef is None:
                raise ValueError(
                    "SoftQuadrupole: when using field data, all three parameters "
                    "'z', 'field_or_gradient', and 'ncoef' must be provided."
                )
            from ..fourier import fourier_coefficients

            cos_coefficients, sin_coefficients = fourier_coefficients(
                z, field_or_gradient, ncoef
            )

        kwargs = dict(
            ds=ds,
            gscale=gscale,
            dx=dx,
            dy=dy,
            rotation=rotation,
            aperture_x=aperture_x,
            aperture_y=aperture_y,
            mapsteps=mapsteps,
            nslice=nslice,
        )
        if name is not None:
            kwargs["name"] = name
        if cos_coefficients is not None:
            kwargs["cos_coefficients"] = cos_coefficients
        if sin_coefficients is not None:
            kwargs["sin_coefficients"] = sin_coefficients

        _original_init(self, **kwargs)

    _new_init.__doc__ = _original_init.__doc__
    cls.__init__ = _new_init
