.. _examples-rfcavity-wake-sc:

Acceleration by RF Cavities with Space Charge and Wakefields
============================================================

This test is based on examples-rfcavity, except that collective effects are present (both space charge and short-range wakefields).

A 300 pC bunch is accelerated through a sequence of 4 RF cavities.

We use a 230 MeV electron beam with initial normalized rms emittance of 1 um.

The lattice and beam parameters are based on Example 2 of the IMPACT-Z examples folder:

https://github.com/impact-lbl/IMPACT-Z/tree/master/examples/Example2

The final target beam energy and beam moments are based on simulation in
IMPACT-Z, without space charge.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run as:

* **Python** script: ``python3 run_rfcavity_wake_sc.py``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_rfcavity_wake_sc.py
          :language: python3
          :caption: You can copy this file from ``examples/rfcavity/run_rfcavity_wake_sc.py``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_rfcavity_wake_sc.py``

   .. literalinclude:: analysis_rfcavity_wake_sc.py
      :language: python3
      :caption: You can copy this file from ``examples/rfcavity/analysis_rfcavity_wake_sc.py``.
