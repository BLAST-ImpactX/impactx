.. _usage-howto-python-extend:

Extend a Simulation with Python
===============================

Overview
--------

ImpactX's Python interface lets you weave custom Python code into a running simulation.
Through this interface, you can **access and modify simulation data**, e.g., particle phase-space coordinates, the reference particle, and space-charge fields, as the simulation is tracking the beam.
This facilitates a wide range of workflows, including:

   - **In-situ diagnostics** that compute quantities of interest at each turn or element without writing files to disk.
   - **Dynamic lattice manipulation**, e.g., changing element strengths between turns, ramping voltages, or feedback systems.
   - **Custom beam optics elements**, e.g., a non-linear kick or a tabulated map, plugged in via the :py:class:`impactx.elements.Programmable` element.
   - **AI/ML surrogate models** built with PyTorch or TensorFlow to emulate detailed transport in a single element.

Because ImpactX exposes particles and fields *without copies* through `pyAMReX <https://pyamrex.readthedocs.io/>`__,
per-step Python callbacks have very low overhead and stay GPU-resident when ImpactX runs on GPUs
(``cupy``, PyTorch, etc. can operate directly on the device data).

.. _usage-howto-python-extend-run:

How to run a simulation with Python extensions
----------------------------------------------

- **Install ImpactX with the Python interface**: when compiling from source, this means ``-DImpactX_PYTHON=ON``
  (see :ref:`building <install-developers>`). Pre-built :ref:`ImpactX packages <install-users>` also include this interface.

- **Write a Python script that drives the simulation**, registering callbacks or programmable elements before calling :py:meth:`~impactx.ImpactX.track_particles`.
  A minimal template looks like:

  .. code-block:: python

     from impactx import ImpactX, distribution, elements

     sim = ImpactX()
     sim.space_charge = False
     sim.init_grids()

     # ... set up reference particle, distribution, and lattice ...

     # register tracking hooks or programmable elements here:
     def after_each_period(sim):
         turn = sim.tracking_period
         # e.g. compute a custom diagnostic from sim.beam

     sim.hook["after_period"] = after_each_period

     # advance the simulation
     sim.track_particles()
     sim.finalize()

- **Run the script** like any Python script, optionally under MPI:

  .. code-block:: bash

     mpirun -np <n_ranks> python <script>.py

.. _usage-howto-python-extend-callbacks:

Callback Functions
------------------

ImpactX provides two complementary callback surfaces.

Tracking hooks
^^^^^^^^^^^^^^

The :py:attr:`~impactx.ImpactX.hook` mapping registers user-defined functions that are invoked at well-defined
points of the tracking loop. The supported locations are:

* ``"before_period"`` — before each lattice period (e.g. turn or channel period)
* ``"after_period"`` — after each lattice period
* ``"before_element"`` — before each lattice element is entered
* ``"after_element"`` — after each lattice element is exited
* ``"before_slice"`` — before each element slice (and before that slice's space-charge solve, when enabled; to read freshly-computed space-charge fields, prefer ``"after_element"``)

Each hook receives the :py:class:`~impactx.ImpactX` instance and may inspect or modify simulation state
through it (e.g., ``sim.beam``, ``sim.lattice``, ``sim.rho``, ``sim.phi``, ...).
Use :py:attr:`~impactx.ImpactX.tracking_step`, :py:attr:`~impactx.ImpactX.tracking_period`, and
:py:attr:`~impactx.ImpactX.tracking_element` to identify *where* in the loop the hook was called.
For instance, you can ``return`` if the current element does not match a type or name.

.. code-block:: python

   def after_each_turn(sim):
       turn = sim.tracking_period
       element = sim.tracking_element
       df = sim.beam.to_df(local=True)
       if df is not None:
           print(f"turn {turn} at element={element.name}: "
                 f"{len(df)} local particles, "
                 f"<x>={df['position_x'].mean():.3e} m")

   sim.hook["after_period"] = after_each_turn

Full example: :ref:`Acceleration by RF Cavities <examples-rfcavity-ref-part-hook>`.

Programmable lattice element
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For workflows that *replace* a beamline element entirely (rather than observing the beam between elements),
use the :py:class:`impactx.elements.Programmable` element. It exposes three property hooks:

* :py:attr:`~impactx.elements.Programmable.push`: push the whole particle container.
* :py:attr:`~impactx.elements.Programmable.ref_particle`: push only the reference particle.
* :py:attr:`~impactx.elements.Programmable.beam_particles`: called once per particle tile (high-performance);
  the reference particle is updated *before* this callback fires.

Typical usage (per-tile push) looks like this:

.. code-block:: python

   from impactx import elements

   def my_drift(pge, pti, refpart):
       """Push the beam particles as a drift over one slice."""
       soa = pti.soa().to_xp()   # NumPy (CPU) or CuPy (GPU)
       x = soa.real["position_x"]
       y = soa.real["position_y"]
       t = soa.real["position_t"]
       # ...

       slice_ds = pge.ds / pge.nslice
       betgam2 = refpart.pt**2 - 1.0

       x[:] += slice_ds * px[:]
       y[:] += slice_ds * py[:]
       t[:] += (slice_ds / betgam2) * pt[:]


   def my_ref_drift(pge, refpart):
       """Push the reference particle over one slice."""
       x = refpart.x
       px = refpart.px
       s = refpart.s

       slice_ds = pge.ds / pge.nslice

       # ...
       refpart.x = x + 1.23
       refpart.s = s + slice_ds


   dr = elements.Programmable(name="my_drift")
   dr.ds = 0.25
   dr.nslice = 5
   dr.beam_particles = lambda pti, refpart: my_drift(dr, pti, refpart)
   dr.ref_particle = lambda refpart: my_ref_drift(dr, refpart)

   sim.lattice.append(pge)

A complete script is provided in this :ref:`FODO Cell example <examples-fodo-programmable>`.
For a more complex example, see our :ref:`ML surrogate element example <examples-ml-surrogate>`.

.. _usage-howto-python-extend-data-access:

Accessing simulation data through Python
----------------------------------------

While the simulation is running, the Python code (i.e., the code inside callback functions) has read and
write access to the beam particles and the space-charge fields. The two sections below describe the
specific Python syntax for each.

.. toctree::
   :maxdepth: 1

   python_particle_data
   python_field_data

See also :ref:`dataanalysis` for file and streaming-based analysis of the openPMD output written by :py:class:`~impactx.elements.BeamMonitor`.
