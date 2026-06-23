.. _examples-solenoid-restart:

Solenoid channel (Restart)
==========================

This is the same example as the :ref:`proton beam undergoing 90 deg X-Y rotation in an ideal solenoid channel <examples-solenoid>`, but it restarts the resulting beam and rotates it another 3 channel periods to the initial X-Y conditions.

The solenoid magnetic field corresponds to B = 2 T.

The second moments of the transverse particle distribution after the solenoid channel are rotated by 90 (start) + 270 = 360 degrees:  the final transverse moments should coincide with the
initial transverse moments to within the level expected due to noise due to statistical sampling.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_solenoid.py`` or
* ImpactX **executable** using an input file: ``impactx input_solenoid.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_solenoid.py
          :language: python3
          :caption: You can copy this file from ``examples/solenoid_restart/run_solenoid.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_solenoid.in
          :language: ini
          :caption: You can copy this file from ``examples/solenoid_restart/input_solenoid.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_solenoid.py``

   .. literalinclude:: analysis_solenoid.py
      :language: python3
      :caption: You can copy this file from ``examples/solenoid_restart/analysis_solenoid.py``.
