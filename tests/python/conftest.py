# -*- coding: utf-8 -*-

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


@pytest.fixture(autouse=True, scope="function")
def amrex_init(tmpdir):
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
        if amr.initialized():
            amr.finalize()
