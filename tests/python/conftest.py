# -*- coding: utf-8 -*-

import gc
import os

import pytest

import amrex.space3d as amr
import impactx

if impactx.Config.have_mpi:
    from mpi4py import MPI  # noqa

    print("loaded mpi4py")
else:
    print("NO mpi4py load")

# base path for input files
basepath = os.getcwd()


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "manages_amrex: test manages its own AMReX init/finalize lifecycle"
    )


# what the setup and the call phase of the current test raised
raised_by_test = []


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    """Remember what the setup and the call phase of a test raised.

    This is an error in the test or in a fixture, and the outcome of
    ``pytest.fail``, ``pytest.skip``, ``pytest.xfail`` and of a ``pytest.raises``
    or ``pytest.warns`` that did not see what it expects.
    """
    if call.when == "setup":
        raised_by_test.clear()  # a new test: drop what the last one raised
    if call.excinfo is not None and call.when != "teardown":
        raised_by_test.append(call.excinfo.value)
    yield


def release_failed_test_objects():
    """Destroy the objects that a failed test left behind.

    The traceback of what a test raised holds the local variables of the test,
    e.g., the simulation, its particle containers and MultiFabs. pytest keeps
    it beyond the test: an ``Exception`` of the call phase in ``sys.last_exc``,
    ``sys.last_value`` and ``sys.last_traceback`` until the next test is called,
    and whatever a test raised in a reference cycle until the garbage collector
    runs. These objects must be destroyed before AMReX is finalized: freed
    later, they release their memory into an AMReX that is gone or was
    initialized anew for the next test, which aborts the whole pytest run on
    GPU and finalizes AMReX under the next test on CPU.

    A test that raised nothing is left alone: what it keeps in a reference
    cycle, e.g., with the ``excinfo`` of a ``pytest.raises``, lives until the
    garbage collector runs. Not released either are the variables of
    module-level code, e.g., of a script that a test runs with ``exec``, and
    what the teardown of another fixture raised.
    """
    pending = list(raised_by_test)

    # the frames in the tracebacks of these exceptions, of the exceptions they
    # are chained to and of the exceptions they group: callers before callees
    frames = []
    seen = set()
    while pending:
        exc = pending.pop()
        if exc is None or id(exc) in seen:
            continue
        seen.add(id(exc))
        tb = exc.__traceback__
        while tb is not None:
            frames.append(tb.tb_frame)
            tb = tb.tb_next
        pending += [exc.__cause__, exc.__context__]
        # an ExceptionGroup holds the exceptions that it groups in a tuple
        grouped = getattr(exc, "exceptions", ())
        if isinstance(grouped, tuple):
            pending += [e for e in grouped if isinstance(e, BaseException)]
    if not frames:
        return  # the test raised nothing

    # Callees before their callers, as when Python unwinds the stack: a callee
    # can hold what depends on the simulation of its caller, and a simulation
    # that is destroyed finalizes AMReX.
    while frames:
        frame = frames.pop()
        try:
            # drops the local variables of a frame that finished executing
            frame.clear()
        except RuntimeError:
            pass  # still executing, e.g., a frame of pytest itself
        else:
            # Up to Python 3.12, a frame caches the dict that f_locals
            # returned, pytest reads it to report the failure, and clear()
            # leaves it as is. Reading f_locals again updates this dict.
            frame.f_locals

    # objects in reference cycles are only destroyed by the garbage collector
    gc.collect()


@pytest.fixture(autouse=True, scope="function")
def amrex_init(request, tmpdir):
    # Allow tests to manage their own AMReX lifecycle
    if request.node.get_closest_marker("manages_amrex"):
        with tmpdir.as_cwd():
            yield
            release_failed_test_objects()
            # a failed test may not have got to finalize AMReX itself
            if amr.initialized():
                amr.finalize()
        return

    with tmpdir.as_cwd():
        # warning: with an external AMReX initialize, our ImpactX overwrite_amrex_parser_defaults is never called!
        amr.initialize(
            [
                # print AMReX status messages
                # consider also 0 (silent) and 2 (FabArray and TileArray/FB/Copy/FillPatch/CrsFineCache usage)
                "amrex.verbose=1",
                # disable verbose profiler plots at the end of each test
                "tiny_profiler.enabled=0",
                # throw exceptions and create core dumps instead of
                # AMReX backtrace files: allows to attach to
                # debuggers
                "amrex.throw_exception=1",
                "amrex.signal_handling=0",
                # abort GPU runs if out-of-memory instead of swapping to host RAM
                "amrex.abort_on_out_of_gpu_memory=1",
                # allocate GPU memory on-demand instead of pre-allocating 3/4th
                # to enable parallel test runs on the same GPU
                # https://amrex-codes.github.io/amrex/docs_html/RuntimeParameters.html?highlight=arena#memory
                "amrex.the_arena_init_size=0",
                # We override the default tiling option for particles, which is always
                # "false" in AMReX, to "false" if running on GPU and "true" on CPU.
                f"particles.do_tiling={0 if impactx.Config.have_gpu else 1}",
            ]
        )
        yield
        release_failed_test_objects()
        if amr.initialized():
            amr.finalize()
