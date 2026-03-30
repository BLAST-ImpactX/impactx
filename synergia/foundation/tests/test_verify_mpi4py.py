import pytest


def test_mpi4py_installed():
    """Make sure that mpi4py is installed."""
    try:
        import mpi4py
    except ImportError:
        print("failed to import mpi4py")
        assert False

    try:
        from mpi4py import MPI
    except ImportError:
        print("failed to import MPI from mpi4py")
        assert False
