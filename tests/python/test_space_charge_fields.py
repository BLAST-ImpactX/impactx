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

    from impactx import ImpactX, distribution, elements

    sim = ImpactX()

    # numerical parameters
    sim.n_cell = [48, 48, 40]
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

    # capture phi statistics so we can sanity-check below
    phi_stats = []

    def plot_phi(sim):
        """Plot phi at the transverse mid-plane after an element finishes.

        ``sim.hook["after_element"]`` fires after the element's last per-slice
        Poisson solve, so ``sim.phi(lev=0)`` holds the freshly computed
        potential of the last slice when we read it.
        (``"before_slice"`` would fire *before* the slice's solve, when phi
        is from a previous slice's solve — or uninitialized on the first slice.)

        GPU/MPI-portable: we iterate over boxes and use ``to_numpy(copy=True)``
        to bring host-side data we can plot.
        """
        phi = sim.phi(lev=0)
        geom = sim.Geom(lev=0)
        dr = geom.data().CellSize()
        ng = phi.n_grow_vect
        half_z = sim.n_cell[2] // 2

        for mfi in phi:
            bx = mfi.validbox()
            arr = phi.array(mfi).to_numpy(copy=True)  # indices x, y, z, comp

            # shift box to zero-based local mfi index space
            half_z_local = half_z - bx.lo_vect[2]
            bx.shift(bx.lo_vect * -1)
            # skip tiles that don't contain the mid-plane
            if half_z_local < 0 or half_z_local >= arr.shape[2]:
                continue

            # strip ghost cells in the transverse plane
            plane = arr[ng[0] : -ng[0], ng[1] : -ng[1], half_z_local, 0]
            phi_stats.append(float(plane.max()))

            if save_png:
                fig, ax = plt.subplots()
                im = ax.imshow(
                    plane.T * 1.0,
                    origin="lower",
                    aspect="auto",
                    extent=[
                        bx.lo_vect[0] * dr[0],
                        (bx.lo_vect[0] + plane.shape[0]) * dr[0],
                        bx.lo_vect[1] * dr[1],
                        (bx.lo_vect[1] + plane.shape[1]) * dr[1],
                    ],
                )
                fig.colorbar(im, label=r"$\phi$  [V]")
                ax.set_xlabel(r"$x$  [m]")
                ax.set_ylabel(r"$y$  [m]")
                fig.savefig(f"phi_step{sim.tracking_step:04d}.png")
                plt.close(fig)

    sim.hook["after_element"] = plot_phi

    sim.track_particles()
    sim.finalize()
    # [doc-end] cfchannel-phi-hook

    # the hook should have fired at least once and seen a non-trivial potential
    assert len(phi_stats) > 0
    assert max(phi_stats) > 0.0


if __name__ == "__main__":
    import amrex.space3d as amr

    amr.initialize([])
    test_space_charge_phi_via_hook(save_png=True)
    if amr.initialized():
        amr.finalize()
