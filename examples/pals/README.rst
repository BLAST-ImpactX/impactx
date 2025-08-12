.. _examples-fodo-pals:

FODO Cell with a PALS Lattice
=============================

This is a duplicate of the :ref:`FODO cell example <examples-fodo>`, except that the lattice file is specified in `PALS <https://pals-project.readthedocs.io>`__-compliant format.

The second moments of the particle distribution after the FODO cell should coincide with the second moments of the particle distribution before the FODO cell, to within the level expected due to noise due to statistical sampling.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run as:

* **Python** script: ``python3 run_fodo_pals.py``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_fodo_pals.py
          :language: python3
          :caption: You can copy this file from ``examples/pals/run_fodo_pals.py``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_fodo_pals.py``

   .. literalinclude:: analysis_fodo_pals.py
      :language: python3
      :caption: You can copy this file from ``examples/pals/analysis_fodo_pals.py``.


Visualize
---------
