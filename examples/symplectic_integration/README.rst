.. _examples-exact-quad:

Symplectic Integration in an Exact Quadrupole
==============================================

This benchmark tests the use of the ExactQuad (quadrupole_exact) element for integrating through a quadrupole using the exact nonlinear Hamiltonian.

A 25 pC electron bunch with 100 MeV total energy, a small initial rms beam size of (3.9, 3.9, 1.0) microns, 2 mrad transverse divergence and 2.5% energy spread undergoes rapid expansion followed
by transverse focusing in a quadrupole doublet.  The parameters are chosen such that chromatic focusing effects are important.

In this test, the initial and final values of :math:`\lambda_x`, :math:`\lambda_y`, :math:`\lambda_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.

In addition, the Hamiltonian value is computed for each particle at the entrance and exit of the final quadrupole.  The change in the Hamiltonian value, taken relative to the standard deviation
:math:`\sigma_H` over particles, must be smaller than the allowed tolerance (here, taken to be 0.1%).

Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_exact_quad.py`` or
* ImpactX **executable** using an input file: ``impactx input_exact_quad.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_exact_quad.py
          :language: python3
          :caption: You can copy this file from ``examples/symplectic_integration/run_exact_quad.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_exact_quad.in
          :language: ini
          :caption: You can copy this file from ``examples/symplectic_integration/input_exact_quad.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_exact_quad.py``

   .. literalinclude:: analysis_exact_quad.py
      :language: python3
      :caption: You can copy this file from ``examples/symplectic_integration/analysis_exact_quad.py``.
