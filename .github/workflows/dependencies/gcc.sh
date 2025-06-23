#!/usr/bin/env bash
#
# Copyright 2021-2023 The ImpactX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl

set -eu -o pipefail

sudo apt -qqq update
sudo apt install -y     \
    build-essential     \
    ca-certificates     \
    ccache              \
    cmake               \
    gnupg               \
    libfftw3-dev        \
    libhdf5-dev         \
    ninja-build         \
    pkg-config          \
    python3             \
    python3-pip         \
    wget

python3 -m pip install -U pip pipx uv

uv pip install --system -U build cmake packaging setuptools[core] wheel
uv pip install --system -U -r requirements.txt
uv pip install --system -U -r src/python/impactx/dashboard/requirements.txt
uv pip install --system -U -r examples/requirements.txt
uv pip install --system -U -r tests/python/requirements.txt

# extra tests
uv pip install --system -U -r examples/requirements_torch_cpu.txt
uv pip install --system -U openPMD-validator
