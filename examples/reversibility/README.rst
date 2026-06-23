.. _examples-zero-bend-inverse:

Bending Dipole with Zero Field Strength
========================================

This example tests the effect of a bending dipole with zero field strength, as represented using the elements Sbend and ExactSbend.

In the limiting case of zero field strength, the effect of the dipole should be identical to that of a drift of the same length.

To test this, the map for a dipole is applied, and this is followed by applying the inverse map for a drift of the corresponding length.

We use a 2 GeV electron beam, identical to the distribution used in the test "examples-rotation".

The second moments of x, y, and t should be exactly unchanged.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_zero_bend_inverse.py`` or
* ImpactX **executable** using an input file: ``impactx input_zero_bend_inverse.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_zero_bend_inverse.py
          :language: python3
          :caption: You can copy this file from ``examples/reversibility/run_zero_bend_inverse.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_zero_bend_inverse.in
          :language: ini
          :caption: You can copy this file from ``examples/reversibility/input_zero_bend_inverse.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_zero_bend_inverse.py``

   .. literalinclude:: analysis_zero_bend_inverse.py
      :language: python3
      :caption: You can copy this file from ``examples/reversibility/analysis_zero_bend_inverse.py``.
