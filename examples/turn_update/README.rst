.. _examples-turn-update:

Turn-Dependend Element Changes
==============================

10 periods of drift elements, where each drift is a stand-in for a "turn" map.
On each turn, the next turn map (drift) is modified based in particle properties.

A 2 GeV electron bunch with normalized transverse rms emittance of 50 pm.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run as:

* **Python** script: ``python3 run_turn_update.py`` or
For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_turn_update.py
          :language: python3
          :caption: You can copy this file from ``examples/turn_update/run_turn_update.py``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_turn_update.py``

   .. literalinclude:: analysis_turn_update.py
      :language: python3
      :caption: You can copy this file from ``examples/turn_update/analysis_turn_update.py``.
