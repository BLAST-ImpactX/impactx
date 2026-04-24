import numpy as np

from .bunch import *


def get_particles_numpy(self, group=ParticleGroup.regular):
    return np.array(self.get_host_particles(group), copy=False)


setattr(Bunch, "get_particles_numpy", get_particles_numpy)


def get_particle_masks_numpy(self, group=ParticleGroup.regular):
    return np.array(self.get_host_particle_masks(group), copy=False)


setattr(Bunch, "get_particle_masks_numpy", get_particle_masks_numpy)


def calculate_mean(b):
    return np.array(Core_diagnostics.calculate_mean_ka(b))


def calculate_abs_mean(b):
    return np.array(Core_diagnostics.calculate_abs_mean_ka(b))


def calculate_std(b, mean):
    # mean = Core_diagnostics.calculate_mean_ka(b)
    return np.array(Core_diagnostics.calculate_std_ka(b, mean))


def calculate_mom2(b, mean):
    # mom2 = Core_diagnostis.calculate_mom2_ka(b)
    return np.array(Core_diagnostics.calculate_mom2_ka(b, mean))


setattr(Core_diagnostics, "calculate_mean", staticmethod(calculate_mean))
setattr(Core_diagnostics, "calculate_abs_mean", staticmethod(calculate_abs_mean))
setattr(Core_diagnostics, "calculate_std", staticmethod(calculate_std))
setattr(Core_diagnostics, "calculate_mom2", staticmethod(calculate_mom2))
