import os

# Python 3.8+ on Windows: DLL search paths for dependent
# shared libraries
# Refs.:
# - https://github.com/python/cpython/issues/80266
# - https://docs.python.org/3.8/library/os.html#os.add_dll_directory
if os.name == "nt":
    # add anything in the current directory
    pwd = __file__.rsplit(os.sep, 1)[0] + os.sep
    os.add_dll_directory(pwd)
    # add anything in PATH
    paths = os.environ.get("PATH", "")
    for p in paths.split(";"):
        p_abs = os.path.abspath(os.path.expanduser(os.path.expandvars(p)))
        if os.path.exists(p_abs):
            os.add_dll_directory(p_abs)

# import core bindings to C++
from . import impactx_pybind
from .impactx_pybind import *  # noqa
from .madx_to_impactx import read_beam  # noqa

__version__ = impactx_pybind.__version__
__doc__ = impactx_pybind.__doc__
__license__ = impactx_pybind.__license__
__author__ = impactx_pybind.__author__

from .distribution_input_helpers import twiss  # noqa
from .fourier import fourier_coefficients  # noqa
from .twiss_lattice import compute_twiss as _compute_twiss_lattice
from .extensions.KnownElementsList import (
    FilteredElementsList,
    register_KnownElementsList_extension,
)
from .extensions.ImpactXParticleContainer import (
    register_ImpactXParticleContainer_extension,
)
from .extensions.RFCavity import (
    register_RFCavity_extension,
)
from .extensions.SoftQuadrupole import (
    register_SoftQuadrupole_extension,
)
from .extensions.SoftSolenoid import (
    register_SoftSolenoid_extension,
)

# at this place we can enhance Python classes with additional methods written
# in pure Python or add some other Python logic

# MAD-X file reader for beamline lattice elements
register_KnownElementsList_extension(impactx_pybind.elements.KnownElementsList)

# Public alias on the elements submodule (same class object as in extensions)
impactx_pybind.elements.FilteredElementsList = FilteredElementsList
FilteredElementsList.__module__ = impactx_pybind.elements.__name__

# MAD-X file reader for reference particle
RefPart.load_file = read_beam  # noqa

# Pure Python extensions to ImpactX types
register_ImpactXParticleContainer_extension(impactx_pybind.ImpactXParticleContainer)

# Alternative constructors from raw field data
register_RFCavity_extension(impactx_pybind.elements.RFCavity)
register_SoftQuadrupole_extension(impactx_pybind.elements.SoftQuadrupole)
register_SoftSolenoid_extension(impactx_pybind.elements.SoftSolenoid)

# Expose the Wolski-6D lattice Twiss analysis as a convenience method on the
# ImpactX simulation object. Internally this routes through
# ``sim.lattice.map_trace(sim.beam.ref)`` (the reference particle is copied
# by pybind11 on the way in, so the caller's ``sim.beam.ref`` is preserved),
# and the linear-optics analysis lives in ``impactx.twiss_lattice``.
#
#   sim.twiss()                               -> periodic ring: matched Twiss
#                                                from the Wolski eigendecomposition
#                                                of the one-turn map
#   sim.twiss(init={"beta_x": ..., ...})      -> transfer line / linac: propagate
#                                                uncoupled Courant-Snyder Twiss
#                                                from the given initial conditions
#
# Returns a dict of numpy arrays (one entry per element boundary) with beta,
# alpha, phase advance mu, and dispersion along s, plus the fractional tunes
# for rings. The underlying lattice-level building blocks are available as
# ``sim.lattice.transfer_map(ref)`` (end-to-end map) and
# ``sim.lattice.map_trace(ref)`` (per-element cumulative trace).
impactx_pybind.ImpactX.twiss = _compute_twiss_lattice
