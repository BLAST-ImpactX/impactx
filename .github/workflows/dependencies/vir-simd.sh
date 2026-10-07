#!/usr/bin/env bash
#
# Copyright 2022-2026 The ImpactX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl

set -eu -o pipefail

# vir-simd, needed for ImpactX_SIMD=ON
#   built outside of the source tree, so it cannot end up in a commit
cd "$(mktemp -d)"
wget https://github.com/mattkretz/vir-simd/archive/refs/tags/v0.4.4.tar.gz
tar -xf v0.4.4.tar.gz
cmake -S vir-simd-0.4.4 -B vir-simd-build
sudo cmake --build vir-simd-build --target install
