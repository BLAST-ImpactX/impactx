.. _examples-phase-wrap:

Test of Periodic Phase-Wrapper
==============================

This example tests the application of a periodic (longitudinal) particle boundary condition.

For systems involving bunches that are long relative to the local size of the RF bucket, it is often necessary to capture the effect of particle
slippage between adjacent buckets. To handle this effectively without tracking multiple bunches, the longitudinal phase can be treated modulo 2pi.
This requires applying a wrap-around mapping to the temporal coordinate t, to handle particles that move outside the RF bucket.

In this test, the bucket length (in m) is specified as a property of the beam.  During tracking, two monitor elements are specified.
Since the boundary condition is applied once per slice, the boundary condition here is applied between the two monitors.

For the analysis comparison, the longitudinal coordinates found in the initial particle distribution are updated after application of the periodic
boundary condition in Python, and the result is compared against the final simulated particle data.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_phase_wrap.py`` or
* ImpactX **executable** using an input file: ``impactx input_phase_wrap.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_phase_wrap.py
          :language: python3
          :caption: You can copy this file from ``examples/phase_wrapper/run_phase_wrap.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_phase_wrap.in
          :language: ini
          :caption: You can copy this file from ``examples/phase_wrapper/input_phase_wrap.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_phase_wrap.py``

   .. literalinclude:: analysis_phase_wrap.py
      :language: python3
      :caption: You can copy this file from ``examples/phase_wrapper/analysis_phase_wrap.py``.
