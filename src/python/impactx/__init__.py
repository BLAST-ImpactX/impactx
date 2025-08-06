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
from .extensions.KnownElementsList import (
    register_KnownElementsList_extension,
)
from .extensions.ImpactXParticleContainer import (
    register_ImpactXParticleContainer_extension,
)

# at this place we can enhance Python classes with additional methods written
# in pure Python or add some other Python logic

# MAD-X file reader for beamline lattice elements
register_KnownElementsList_extension(impactx_pybind.elements.KnownElementsList)

# MAD-X file reader for reference particle
RefPart.load_file = read_beam  # noqa

# Pure Python extensions to ImpactX types
register_ImpactXParticleContainer_extension(impactx_pybind.ImpactXParticleContainer)
