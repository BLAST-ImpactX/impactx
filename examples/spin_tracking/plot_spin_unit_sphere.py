#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import argparse
import glob
import re
       
import matplotlib.pyplot as plt
import numpy as np
import openpmd_api as io
import pandas as pd
from matplotlib.ticker import MaxNLocator
# Required for 3D projection
from mpl_toolkits.mplot3d import Axes3D

# Collect data for the initial/final beam
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
last_step = list(series.iterations)[-1]
initial = series.iterations[1].particles["beam"].to_df()
final = series.iterations[last_step].particles["beam"].to_df()

# Create a new figure and add a 3D subplot
fig, axs = plt.subplots(1, 2, subplot_kw={'projection': '3d'})

# Define the spherical coordinates
# u and v are parameters for the surface
u = np.linspace(0, 2 * np.pi, 100) # Azimuth angle, ranges from 0 to 2*pi
v = np.linspace(0, np.pi, 100)     # Elevation angle, ranges from 0 to pi

# Create a meshgrid from u and v
u, v = np.meshgrid(u, v)

# Calculate the Cartesian coordinates for a unit sphere (radius r=1)
# x = cos(u) * sin(v)
# y = sin(u) * sin(v)
# z = cos(v)
x = np.cos(u) * np.sin(v)
y = np.sin(u) * np.sin(v)
z = np.cos(v)

ax1 = axs[0]
ax2 = axs[1]

# Plot the surface
# alpha controls opacity (0.3 makes it somewhat transparent)
ax1.plot_surface(x, y, z, color='b', alpha=0.3)
ax2.plot_surface(x, y, z, color='b', alpha=0.3)

# Set axis labels
# Rotate axes so that z is longitudinal and y is vertical:
ax1.set_xlabel('$S_z$', fontsize=14, weight='bold')
ax1.set_ylabel('$S_x$', fontsize=14, weight='bold')
ax1.set_zlabel('$S_y$', fontsize=14, weight='bold')
ax1.set_title('Initial Particle Spin')
ax2.set_xlabel('$S_z$', fontsize=14, weight='bold')
ax2.set_ylabel('$S_x$', fontsize=14, weight='bold')
ax2.set_zlabel('$S_y$', fontsize=14, weight='bold')
ax2.set_title('Final Particle Spin')

# Ensure the aspect ratio is equal so the sphere looks like a sphere and not an ellipsoid
# This involves manually setting limits based on the maximum extent of the data.
max_range = 1.0
ax1.set_xlim([-max_range, max_range])
ax1.set_ylim([-max_range, max_range])
ax1.set_zlim([-max_range, max_range])
ax2.set_xlim([-max_range, max_range])
ax2.set_ylim([-max_range, max_range])
ax2.set_zlim([-max_range, max_range])

# Set ticks
ax1.set_xticks([-1.0, -0.5, 0.0, 0.5, 1.0])
ax1.set_yticks([-1.0, -0.5, 0.0, 0.5, 1.0])
ax1.set_zticks([-1.0, -0.5, 0.0, 0.5, 1.0])
ax2.set_xticks([-1.0, -0.5, 0.0, 0.5, 1.0])
ax2.set_yticks([-1.0, -0.5, 0.0, 0.5, 1.0])
ax2.set_zticks([-1.0, -0.5, 0.0, 0.5, 1.0])

# scatter plot of initial particle spin on the sphere
ax1.scatter3D(
    initial.spin_z,
    initial.spin_x,  
    initial.spin_y,
    c='k',
    s=0.03,
)

# scatter plot of final particle spin on the sphere
ax2.scatter3D(
    final.spin_z,
    final.spin_x,
    final.spin_y,
    c='k',
    s=0.03,
)

# Adjust plot location and spacing
plt.subplots_adjust(left=0.05,wspace=0.3)

# Display the plot
plt.show()
