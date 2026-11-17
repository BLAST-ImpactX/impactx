#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import os

import matplotlib.image as mpimg
import numpy as np

# one phase space frame after each of the 5 elements of the FODO cell
frames = sorted(os.listdir("diags/fodo_ascent_frames"))
print(f"frames: {frames}")
assert frames == [f"frame_{i:06d}.png" for i in range(5)]

# rendered images, which are not blank
for path in [os.path.join("diags/fodo_ascent_frames", f) for f in frames] + [
    "diags/fodo_ascent_beam3d.png"
]:
    pixels = mpimg.imread(path)
    print(f"{path}: {pixels.shape}")
    assert np.ptp(pixels[..., :3]) > 0.0

# the interactive animation embeds all frames
with open("diags/fodo_ascent.html") as f:
    html = f.read()
assert html.count("<img") >= 1
