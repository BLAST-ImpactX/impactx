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
    clang               \
    cmake               \
    gnupg               \
    libc++-dev          \
    libc++abi-dev       \
    libfftw3-dev        \
    libfftw3-mpi-dev    \
    libhdf5-openmpi-dev \
    libomp-dev          \
    libopenmpi-dev      \
    ninja-build         \
    pkg-config          \
    python3             \
    python3-pip         \
    wget

curl -LsSf https://astral.sh/uv/install.sh | sh
export PATH=$HOME/.local/bin:$PATH

sudo uv pip install --system -U build cmake packaging setuptools[core] wheel
sudo uv pip install --system -U -r requirements_mpi.txt
sudo uv pip install --system -U -r src/python/impactx/dashboard/requirements.txt
sudo uv pip install --system -U -r examples/requirements.txt
sudo uv pip install --system -U -r tests/python/requirements.txt

# cmake-easyinstall
#
sudo curl -L -o /usr/local/bin/cmake-easyinstall https://raw.githubusercontent.com/ax3l/cmake-easyinstall/main/cmake-easyinstall
sudo chmod a+x /usr/local/bin/cmake-easyinstall
export CEI_SUDO="sudo"
export CEI_TMP="/tmp/cei"
