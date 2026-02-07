.. _examples-quad-spin:

Spin Depolarization in a Quadrupole
===================================

This example illustrates the decay of the polarization vector (describing the mean of the three spin components) along the vertical y and longitudinal z directions for a beam undergoing
horizontal focusing in a quadrupole.

We use a 250 MeV proton beam with initial unnormalized rms emittance of 1 micron in the horizontal plane, and 2 micron in the vertical plane.

The beam propagates over one horizontal betatron period, to a location where the polarization vector is described by a simple expression.

In this test, the initial and final values of `mean_sx`, `mean_sy`, and `mean_sz` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_quad_spin.py`` or
* ImpactX **executable** using an input file: ``impactx input_quad_spin.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_quad_spin.py
          :language: python3
          :caption: You can copy this file from ``examples/spin_tracking/run_quad_spin.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_quad_spin.in
          :language: ini
          :caption: You can copy this file from ``examples/spin_tracking/input_quad_spin.in``.


Analyze
-------

The analysis can be run using **either** of the following scripts:


.. dropdown:: Script ``analysis_quad_spin.py``

   .. literalinclude:: analysis_quad_spin.py
      :language: python3
      :caption: You can copy this file from ``examples/spin_tracking/analysis_quad_spin.py``.

.. dropdown:: Script ``analysis_quad_spin_rbc.py``

   .. literalinclude:: analysis_quad_spin_rbc.py
      :language: python3
      :caption: You can copy this file from ``examples/spin_tracking/analysis_quad_spin_rbc.py``.


.. _examples-sbend-spin:

Spin Depolarization in a Dipole
===============================

This example illustrates the decay of the polarization vector (describing the mean of the three spin components) along the horizontal x and longitudinal z directions for a beam undergoing
bending in the x-z plane in a sector dipole.

We use a 2 GeV electron beam.  The beam parameters (in particular, the momentum and energy spread) are artificially large in order to enhance the effect.

The beam propagates over one period, as set by the design spin tune.  By increasing the number of slices, and turning on diagnostics, one can view precession of the polarization vector about
the vertical direction.  A clean precession over 1 period becomes visible when the beam size, momentum spread, and energy spread are set to small values.

In this test, the initial and final values of `mean_sx`, `mean_sy`, and `mean_sz` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_sbend_spin.py`` or
* ImpactX **executable** using an input file: ``impactx input_sbend_spin.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_sbend_spin.py
          :language: python3
          :caption: You can copy this file from ``examples/spin_tracking/run_sbend_spin.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_sbend_spin.in
          :language: ini
          :caption: You can copy this file from ``examples/spin_tracking/input_sbend_spin.in``.



Analyze
-------

The analysis can be run using **either** of the following script:


.. dropdown:: Script ``analysis_sbend_spin.py``

   .. literalinclude:: analysis_sbend_spin.py
      :language: python3
      :caption: You can copy this file from ``examples/spin_tracking/analysis_sbend_spin.py``.
