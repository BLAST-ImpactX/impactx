.. _usage-howto-python-field-data:

Accessing Field Data
====================

ImpactX exposes the **space-charge fields** computed on the AMR mesh:

  - the charge density :math:`\rho`
  - the scalar potential :math:`\phi`, and
  - the resulting force kick vector components,

directly to Python as `pyAMReX <https://pyamrex.readthedocs.io/en/latest/usage/api.html#amrex.space3d.MultiFab>`__ ``MultiFab`` objects.

Selecting a field
-----------------

Use the following accessors on the :py:class:`~impactx.ImpactX` instance, each of which returns
a ``MultiFab`` on the requested mesh-refinement level:

.. code-block:: python

   rho = sim.rho(lev=0)                              # charge density
   phi = sim.phi(lev=0)                              # scalar potential
   E_x = sim.space_charge_field(lev=0, comp="x")     # x-component of the space-charge force
   E_y = sim.space_charge_field(lev=0, comp="y")
   E_z = sim.space_charge_field(lev=0, comp="z")

Available mesh levels range from ``0`` to :py:attr:`sim.finest_level <impactx.ImpactX.finest_level>`.

Lifetime / Validity
-------------------

.. important::

   The space-charge fields are computed *during* the per-slice space-charge solve.
   Outside of that window the ``MultiFab`` data may be uninitialized, stale, or empty.
   You should therefore read them either:

   * From inside a :py:attr:`sim.hook <impactx.ImpactX.hook>` ``"after_element"`` callback,
     which fires after the element's last per-slice solve and sees the freshly computed
     fields of that last slice. Note: ``"before_slice"`` fires *before* the slice's solve,
     so on the very first slice ``phi``/``rho`` are uninitialized.
   * Or right after :py:meth:`~impactx.ImpactX.track_particles` returns, but before the simulation is finalized.

   See :ref:`usage-howto-python-extend-callbacks` for how to install such a hook.
   Space-charge must of course be enabled (see :py:attr:`sim.space_charge <impactx.ImpactX.space_charge>`)
   for the fields to be populated at all.

Mesh metadata
-------------

To interpret indices and array shapes, query the AMR geometry on the same level:

.. code-block:: python

   geom = sim.Geom(lev=0)
   prob_lo = geom.ProbLo()        # physical lower corner of the domain
   prob_hi = geom.ProbHi()        # physical upper corner
   dx      = geom.CellSize()      # cell sizes (dx, dy, dz)

   ba = sim.boxArray(lev=0)        # box decomposition
   dm = sim.DistributionMap(lev=0) # which MPI rank owns which box

Accessing field arrays
----------------------

There are several ways to read or modify a ``MultiFab``, with different trade-offs.
The underlying API is provided by pyAMReX, see the
`pyAMReX compute guide <https://pyamrex.readthedocs.io/en/latest/usage/compute.html#fields>`__
for the full reference.

.. tab-set::

   .. tab-item:: Pre-defined pyAMReX/AMReX methods

      Many reductions and bulk operations are exposed directly as ``MultiFab`` methods:

      .. code-block:: python

         rho_max = rho.max(comp=0)        # maximum value (per component)
         phi_min = phi.min(comp=0)        # minimum value
         rho.mult(2.0, 0, 1)              # scale in place

      These have low overhead and run on CPU or GPU as appropriate.

   .. tab-item:: NumPy-like global indexing

      Single elements and slices can be addressed with global ``(i, j, k)`` indices.
      This is convenient for inspection and debugging.

      .. code-block:: python

         # read a single line of cells (level 0, component 0)
         line = rho[:, ny // 2, nz // 2]

         # write a single value
         rho[i, j, k] = 0.0

      Negative-imaginary indices address ghost cells; see the pyAMReX compute guide linked above.

      In MPI-parallel simulations, this indexing interface will trigger MPI communication and can become costly.

   .. tab-item:: Explicit loop over boxes

      For high-performance, GPU-portable access: modifying a field over many cells, or coupling
      to ``cupy``/``numpy`` arrays in-place.
      Enables to iterate over the boxes owned by the current MPI rank:

      .. code-block:: python

         for mfi in rho:
             bx  = mfi.tilebox()
             arr = rho.array(mfi).to_xp()   # NumPy (CPU) or CuPy (GPU)
             arr[...] *= 0.5                # in-place modification

      :ref:`to_xp() <https://pyamrex.readthedocs.io/en/latest/usage/zerocopy.html>`__ (or the explicit ``to_numpy()``/``to_cupy()``) returns a zero-copy view
      into the box's data with the layout matching the local box, including ghost cells.

Full example: plotting the space-charge potential
-------------------------------------------------

The following script runs a short constant-focusing channel with 3D space charge
and installs an ``"after_element"`` hook that saves a PNG of the scalar potential
:math:`\phi` at the transverse mid-plane after the element's last per-slice
Poisson solve. The same pattern works for ``sim.rho`` and ``sim.space_charge_field``.

This and other examples are available under ``tests/python/``.

.. dropdown:: Full example: ``sim.hook["after_element"]`` reads ``sim.phi``

   .. literalinclude:: ../tests/python/test_space_charge_fields.py
      :language: python
      :start-after: # [doc-start] cfchannel-phi-hook
      :end-before: # [doc-end] cfchannel-phi-hook
      :dedent: 4

See also
--------

* :ref:`usage-howto-python-extend` — where and when to install callbacks for live field access.
* :ref:`usage-howto-python-particle-data` — how to access the beam particles.
* `pyAMReX compute guide <https://pyamrex.readthedocs.io/en/latest/usage/compute.html>`__ — full reference for ``MultiFab`` access.
