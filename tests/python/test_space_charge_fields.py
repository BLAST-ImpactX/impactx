#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-
"""
Access space-charge fields from a Python callback hook.

This test runs a short constant-focusing channel with 3D space charge and
installs a ``sim.hook["before_slice"]`` callback that reads the scalar
potential ``sim.phi(lev=0)`` after each per-slice Poisson solve and saves a
PNG of the transverse mid-plane. The same pattern works for ``sim.rho`` and
``sim.space_charge_field``.

This script is also embedded into the ImpactX documentation via
``literalinclude``, see
``docs/source/usage/howto/python_field_data.rst``.
"""


def test_space_charge_phi_via_hook(save_png=True):
    """Plot the scalar potential phi via a before_slice callback hook."""
    # [doc-start] cfchannel-phi-hook
    import matplotlib.pyplot as plt

    from impactx import ImpactX, amr, distribution, elements

    sim = ImpactX()

    # numerical parameters
    sim.n_cell = [48, 48, 40]
    sim.tiny_profiler = False
    sim.particle_shape = 2
    sim.space_charge = "3D"
    sim.poisson_solver = "fft"
    sim.prob_relative = [1.1]
    sim.slice_step_diagnostics = True

    sim.init_grids()

    # 2 GeV proton beam, ~1 um normalized transverse rms emittance
    kin_energy_MeV = 2.0e3
    bunch_charge_C = 1.0e-8
    npart = 10000

    ref = sim.beam.ref
    ref.set_species("proton").set_kin_energy_MeV(kin_energy_MeV)

    distr = distribution.Waterbag(
        lambdaX=1.2154443728379865788e-3,
        lambdaY=1.2154443728379865788e-3,
        lambdaT=4.0956844276541331005e-4,
        lambdaPx=8.2274435782286157175e-4,
        lambdaPy=8.2274435782286157175e-4,
        lambdaPt=2.4415943602685364584e-3,
    )
    sim.add_particles(bunch_charge_C, distr, npart)

    # short CF cell so the test stays quick
    sim.lattice.extend(
        [
            elements.ConstF(name="cf1", ds=2.0, kx=1.0, ky=1.0, kt=1.0, nslice=5),
        ]
    )

    def plot_phi(sim):
        """Plot phi at the transverse mid-plane after an element finishes.

        ``sim.hook["after_element"]`` fires after the element's last per-slice
        Poisson solve, so ``sim.phi(lev=0)`` holds the freshly computed
        potential of the last slice when we read it.
        (``"before_slice"`` would fire *before* the slice's solve, when phi
        is from a previous slice's solve — or uninitialized on the first slice.)

        GPU/MPI-portable: we use global indexing and pyAMReX garthers to rank zero.
        """
        phi = sim.phi(lev=0)
        geom = sim.Geom(lev=0)
        half_z = sim.n_cell[2] // 2
        mm = 1e3

        # this triggers an MPI-collective gather
        plane = phi[:, :, half_z, 0]

        # other MPI ranks can return now
        if not amr.ParallelDescriptor.IOProcessor():
            return

        # now plot on the IOProcessor
        fig, ax = plt.subplots()
        im = ax.imshow(
            plane.T * 1.0,
            origin="lower",
            aspect="auto",
            extent=[
                geom.ProbLo(0) * mm,
                geom.ProbHi(0) * mm,
                geom.ProbLo(1) * mm,
                geom.ProbHi(0) * mm,
            ],
        )
        fig.colorbar(im, label=r"$\phi$  [V]")
        ax.set_xlabel(r"$x$  [mm]")
        ax.set_ylabel(r"$y$  [mm]")
        if save_png:
            fig.savefig(f"phi_step{sim.tracking_step:04d}.png")
        else:
            plt.show()
        plt.close(fig)

    sim.hook["after_element"] = plot_phi

    sim.track_particles()
    sim.finalize()
    # [doc-end] cfchannel-phi-hook


if __name__ == "__main__":
    test_space_charge_phi_via_hook(save_png=True)
