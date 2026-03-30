#!/usr/bin/env python

import sys

sys.path.append("../../..")
import local_paths

from mpi4py import MPI
from synergia.utils import Commxx
import synergia
import numpy


def test_commxx():
    commxx = Commxx()
    assert commxx.rank() == MPI.COMM_WORLD.Get_rank()
    assert commxx.size() == MPI.COMM_WORLD.Get_size()


# def test_decompose_1d_raw():
#     procs = 4
#     num_per_proc = 5
#     offsets, counts = synergia.utils.decompose_1d_raw(procs, procs * num_per_proc)
#     expected_offset = 0
#     for i in range(0, procs):
#         assert_equal(offsets[i], expected_offset)
#         assert_equal(counts[i], num_per_proc)
#         expected_offset += counts[i]
#
#
# def test_decompose_1d():
#     length = 23
#     offsets, counts = synergia.utils.decompose_1d(Commxx(), length)
#     assert_equal(offsets[0], 0)
#     assert_equal(counts[0], length)
#
#
# def test_decompose_1d_local():
#     length = 17
#     local_length = synergia.utils.decompose_1d_local(Commxx(), length)
#     assert_equal(length, local_length)
