#!/usr/bin/env bash
#
# Copyright 2021-2025 The ImpactX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl

set -eu -o pipefail

sudo apt -qqq update
sudo apt install -y     \
    build-essential     \
    ca-certificates     \
    cmake               \
    gnupg               \
    libhdf5-openmpi-dev \
    libhiredis-dev      \
    libopenmpi-dev      \
    libzstd-dev         \
    ninja-build         \
    openmpi-bin         \
    pkg-config          \
    wget

sudo apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/3bf863cc.pub
echo "deb https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64 /" \
    | sudo tee /etc/apt/sources.list.d/cuda.list

sudo apt update
sudo apt install -y              \
    cuda-command-line-tools-11-7 \
    cuda-compiler-11-7           \
    cuda-cupti-dev-11-7          \
    cuda-minimal-build-11-7      \
    cuda-nvml-dev-11-7           \
    cuda-nvtx-11-7               \
    libcufft-dev-11-7            \
    libcurand-dev-11-7           \
    libcusparse-dev-11-7
sudo ln -s cuda-11.7 /usr/local/cuda

# uv
#
curl -LsSf https://astral.sh/uv/install.sh | sh

# cmake-easyinstall
#
sudo curl -L -o /usr/local/bin/cmake-easyinstall https://raw.githubusercontent.com/ax3l/cmake-easyinstall/main/cmake-easyinstall
sudo chmod a+x /usr/local/bin/cmake-easyinstall
export CEI_SUDO="sudo"
export CEI_TMP="/tmp/cei"

# ccache 4.2+
#
CXXFLAGS="" cmake-easyinstall --prefix=/usr/local \
    git+https://github.com/ccache/ccache.git@v4.10.2 \
    -DCMAKE_BUILD_TYPE=Release        \
    -DENABLE_DOCUMENTATION=OFF        \
    -DENABLE_TESTING=OFF              \
    -DWARNINGS_AS_ERRORS=OFF
