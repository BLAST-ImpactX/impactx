import pytest


def test_synergia_version():
    """Make sure that the version number generated at configuration time
    is the expected one.
    """
    from synergia import version

    assert version.__version__.split("-")[0] == "2024.04.17"


def test_python_version():
    """Make sure that the Python interpreter version as found by CMake at
    configuration time matches the version found by the runtime, when this
    test is executed.
    """
    from synergia import version
    from sys import version_info as vi

    assert version.python_interp_version == f"{vi.major}.{vi.minor}.{vi.micro}"


def test_mpi_version():
    """Make sure that the MPI library found at runtime is the same as that
    which was found by CMake at configuration time.
    """
    from synergia import version
    from mpi4py import MPI

    # For some reason, the module function MPI.Get_library_version returns
    # a string with a trailing null character, which we remove.
    assert (
        version.mpi_library_version.split()
        == MPI.Get_library_version().rstrip("\x00").split()
    )


def test_hdf5_version():
    from synergia import version
    from h5py import version as v

    from_syn_major, from_syn_minor, *_ = version.hdf5_library_version_tuple
    from_h5py_major, from_h5py_minor, *_ = v.hdf5_version_tuple
    assert from_syn_major == from_h5py_major
    assert from_syn_minor == from_h5py_minor
