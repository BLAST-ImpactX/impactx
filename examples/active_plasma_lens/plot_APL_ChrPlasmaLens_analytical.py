#!/usr/bin/env python3
#
# Copyright 2022-2025 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell, Kyrre Sjobak
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

## Plots analytical Twiss parameters
#  as a function of APL gradient [T/m]

import matplotlib.pyplot as plt
import numpy as np
from run_APL import analytic_final_estimate

rigidity_Tm = -0.6688305274603505  # [T*m], 200MeV electrons
APL_length = 20e-3  # [m]

# From "zero" test, sigma_mid = 10 um
# beta_0 = 0.029408806267344052 #[m]
# alpha_0 = 2.5484916642663316  #[-]
# From focusing, sigma_mid = 100 um
beta_0 = 0.39264381163450024
alpha_0 = 0.025484916642663318

# negative g (focusing, rigidity is also negative)
beta = []
alpha = []
g_def = np.linspace(0, -20000)
for APL_g in g_def:
    (beta_end, alpha_end, gamma_end) = analytic_final_estimate(
        APL_g, rigidity_Tm, APL_length, beta_0, alpha_0
    )

    beta.append(beta_end)
    alpha.append(alpha_end)

plt.figure(1)
plt.plot(g_def, beta)
plt.xlabel("g_APL [T/m]")
plt.ylabel("sqrt(beta)")

plt.figure(2)
plt.plot(g_def, alpha)
plt.xlabel("g_APL [T/m]")
plt.ylabel("alpha")


# positive g (defocusing, rigidity is negative)
beta = []
alpha = []
g_def = np.linspace(0, 2000)
for APL_g in g_def:
    (beta_end, alpha_end, gamma_end) = analytic_final_estimate(
        APL_g, rigidity_Tm, APL_length, beta_0, alpha_0
    )

    beta.append(beta_end)
    alpha.append(alpha_end)

plt.figure(1)
plt.plot(g_def, beta)
plt.xlabel("g_APL [T/m]")
plt.ylabel("sqrt(beta)")
plt.grid()
plt.figure(2)
plt.plot(g_def, alpha)
plt.xlabel("g_APL [T/m]")
plt.ylabel("alpha")
plt.grid()

plt.show()
