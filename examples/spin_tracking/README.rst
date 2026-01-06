.. _examples-quad-spin:

Spin Depolarization in a Quadrupole
===================================

This example illustrates the decay of the polarization vector (describing the mean of the three spin components) along the vertical y and longitudinal z directions for a beam undergoing
horizontal focusing in a quadrupole.

We use a 250 MeV proton beam with initial unnormalized rms emittance of 1 micron in the horizontal plane, and 2 micron in the vertical plane.

The beam propagates over one horizontal betatron period, to a location where the polarization vector is described by a simple expression.

In this test, the initial and final values of :math:`\polarization_x`, :math:`\polarization_y`, and :math:`\polarization_z` must agree with nominal values.


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

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_quad_spin.py``

   .. literalinclude:: analysis_quad_spin.py
      :language: python3
      :caption: You can copy this file from ``examples/spin_tracking/analysis_quad_spin.py``.
