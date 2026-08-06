"""
This file is part of ImpactX

Copyright 2026 ImpactX contributors
Authors: Axel Huebl, Chad Mitchell
License: BSD-3-Clause-LBNL
"""


def _as_kick_component(value, num_bins, key, context):
    """Convert one CSR kick component into an AMReX real ``DeviceVector``.

    Accepts NumPy arrays and array-likes (force-cast to the build precision),
    any pyAMReX ``PODVector``, and objects implementing the CUDA array
    interface (e.g., CuPy arrays), converted through the pyAMReX
    ``DeviceVector`` helpers.

    Returns the vector and whether the CUDA array interface was used.
    """
    import numpy as np

    from ..impactx_pybind import Config
    from .ImpactXParticleContainer import _as_real_device_vector

    def error_context():
        return (
            f" for kick component '{key}' in element type "
            f"'{context.element_type}' (name: '{context.element_name}', "
            f"slice {context.slice}); expected {num_bins} values"
        )

    import amrex.space3d as amr

    # pyAMReX PODVectors also export the CUDA array interface on GPU builds:
    # only flag foreign device arrays (e.g., CuPy), which route through CuPy
    is_podvector = type(value).__name__.startswith("PODVector_")
    used_cuda_interface = not is_podvector and hasattr(
        value, "__cuda_array_interface__"
    )
    if not used_cuda_interface and not is_podvector:
        # host arrays and array-likes: force-cast to the build precision
        dtype = np.float64 if Config.precision == "DOUBLE" else np.float32
        try:
            value = np.ascontiguousarray(np.asarray(value, dtype=dtype))
        except (TypeError, ValueError) as e:
            raise RuntimeError(
                "The CSR kick model returned an object that cannot be "
                "converted to a 1D array" + error_context()
            ) from e
        if value.ndim != 1:
            raise RuntimeError(
                "The CSR kick model returned an array that is not 1D" + error_context()
            )

    if isinstance(value, amr.DeviceVector_real):
        # fresh copy: the CSRKick construction consumes (moves from) its
        # inputs, which must not empty a vector owned by the user model
        vec = amr.DeviceVector_real.from_xp(value.to_xp())
    else:
        vec = _as_real_device_vector(value)
    if len(vec) != num_bins:
        raise RuntimeError(
            "The CSR kick model returned an array of wrong length" + error_context()
        )
    return vec, used_cuda_interface


def _wrap_csr_kick_model(model):
    """Wrap a user CSR kick model for the C++ side.

    Validates the model result (a dict with the required key "pt" and the
    optional keys "px" and "py", or a bare array meaning {"pt": array}),
    converts every component into an AMReX real ``DeviceVector`` via pyAMReX,
    and hands them to ImpactX as an ``impactx.CSRKick``.
    """

    def wrapped(profile, context):
        result = model(profile, context)

        if isinstance(result, dict) or hasattr(result, "keys"):
            components = dict(result)
            for key in components:
                if key not in ("pt", "px", "py"):
                    raise RuntimeError(
                        f"The CSR kick model returned an unknown key '{key}' "
                        "(allowed keys: 'pt', 'px', 'py') in element type "
                        f"'{context.element_type}' "
                        f"(name: '{context.element_name}')"
                    )
            if "pt" not in components:
                raise RuntimeError(
                    "The CSR kick model result is missing the required key "
                    f"'pt' in element type '{context.element_type}' "
                    f"(name: '{context.element_name}')"
                )
        else:
            # a bare array is shorthand for {"pt": array}
            components = {"pt": result}

        used_cuda_interface = False
        converted = {}
        for key, value in components.items():
            converted[key], used = _as_kick_component(
                value, profile.num_bins, key, context
            )
            used_cuda_interface = used_cuda_interface or used

        if used_cuda_interface:
            # order device writes issued through CuPy before ImpactX reads
            # the vectors on its own GPU stream
            import cupy

            cupy.cuda.get_current_stream().synchronize()

        from ..impactx_pybind import CSRKick

        return CSRKick(
            pt=converted["pt"],
            px=converted.get("px"),
            py=converted.get("py"),
        )

    return wrapped


def ix_csr_kick_model_get(self):
    """The user-provided CSR kick model callable, or None."""
    return getattr(self, "_csr_kick_model_callable", None)


def ix_csr_kick_model_set(self, model):
    """Install or clear the user-provided CSR kick model."""
    if model is None:
        self._csr_kick_model_callable = None
        self._set_csr_kick_model(None)
        return
    if not callable(model):
        raise TypeError("csr_kick_model must be callable or None")
    # keep the original Python callable alive and introspectable
    self._csr_kick_model_callable = model
    self._set_csr_kick_model(_wrap_csr_kick_model(model))


def register_ImpactX_extension(ix):
    """ImpactX helper methods"""
    ix.csr_kick_model = property(
        ix_csr_kick_model_get,
        ix_csr_kick_model_set,
        doc="User-provided CSR kick model (e.g., an ML surrogate), replacing "
        "the built-in analytic CSR wake model. Called per tracking slice in "
        "CSR-active bend elements as model(profile: CSRProfile, context: "
        "CSRElementContext) and must return per-bin kick forces in Newtons: "
        "a dict with the required key 'pt' and the optional keys 'px'/'py', "
        "or a bare array (= {'pt': array}), each of length profile.num_bins. "
        "Requires csr = True to be applied. Set to None (default) to restore "
        "the built-in analytic model. Under MPI, set it on every rank. It is "
        "invoked on the I/O rank only.",
    )
