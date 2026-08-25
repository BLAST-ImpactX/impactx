#!/usr/bin/env bash
#
# Copyright 2022-2026 The ImpactX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl

set -eu -o pipefail

sudo apt-get -qqq update
sudo apt-get install -y \
    build-essential     \
    ca-certificates     \
    ccache              \
    cmake               \
    gnupg               \
    libboost-dev        \
    libfftw3-dev        \
    libhdf5-dev         \
    ninja-build         \
    pkg-config          \
    python3             \
    python3-pip         \
    wget

# vir-simd
# TODO: back to the release tarball once vir/simd_vecmath.h is in one. It is
#       what makes the SIMD transcendentals call a vector math library instead
#       of evaluating them one lane at a time.
git clone --depth 1 --branch topic-vecmath https://github.com/ax3l/vir-simd.git vir-simd-src
cmake -S vir-simd-src -B vir-simd-build
sudo cmake --build vir-simd-build --target install

python3 -m pip install -U pip
python3 -m pip install -U build packaging setuptools[core] wheel
python3 -m pip install -U cmake
python3 -m pip install -U -r requirements.txt
python3 -m pip install -U -r src/python/impactx/dashboard/requirements.txt
python3 -m pip install -U -r examples/requirements.txt
python3 -m pip install -U -r tests/python/requirements.txt
python3 -m pip install -U pytest-codspeed

# extra tests
python3 -m pip install -U -r examples/requirements_torch_cpu.txt
python3 -m pip install -U openPMD-validator
