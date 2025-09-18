.. _examples-exact-quad:

Symplectic Integration in an Exact Quadrupole
==============================================

This benchmark tests the use of the ExactQuad (quadrupole_exact) element for integrating through a quadrupole using the exact nonlinear Hamiltonian.

A 25 pC electron bunch with 100 MeV total energy, a small initial rms beam size of (3.9, 3.9, 1.0) microns, 2 mrad transverse divergence and 2.5% energy spread undergoes rapid expansion followed
by transverse focusing in a quadrupole doublet.  The parameters are chosen such that chromatic focusing effects are important.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.

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


.. _examples-fodo-multipole:

FODO Channel with Quads Treated as Exact Multipoles
====================================================

This is identical to ``examples-fodo``, except that the quadrupoles have been replaced by ExactMultipole (multipole_exact) elements with only a nonzero quadrupole coefficient.
Its purpose is to test the limiting case of the ExactMultipole model corresponding to a simple quadrupole.

The analysis script is identical to the analysis script used for ``examples-fodo``.

Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_fodo_multipole.py`` or
* ImpactX **executable** using an input file: ``impactx input_fodo_multipole.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_fodo_multipole.py
          :language: python3
          :caption: You can copy this file from ``examples/symplectic_integration/run_fodo_multipole.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_fodo_multipole.in
          :language: ini
          :caption: You can copy this file from ``examples/symplectic_integration/input_fodo_multipole.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_fodo_multipole.py``

   .. literalinclude:: analysis_fodo_multipole.py
      :language: python3
      :caption: You can copy this file from ``examples/symplectic_integration/analysis_fodo_multipole.py``.


.. _examples-long-sextupole:

Symplectic Integration in a Long Sextupole
==============================================

This benchmark tests the use of the ExactMultipole (multipole_exact) element for integrating through a long sextupole using the exact nonlinear Hamiltonian.

An array of initial conditions corresponding to protons with 0.8 GeV total energy is tracked through a 0.5 m long sextupole.

In this test, each particle's final phase space vector is compared against numerical tracking results obtained in PTC (E. Stern).

All 6 phase space coordinates must agree within the specified tolerance.

Run
---

This example can be run as:

* **Python** script: ``python3 run_sextupole.py``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_sextupole.py
          :language: python3
          :caption: You can copy this file from ``examples/symplectic_integration/run_sextupole.py``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_sextupole.py``

   .. literalinclude:: analysis_sextupole.py
      :language: python3
      :caption: You can copy this file from ``examples/symplectic_integration/analysis_sextupole.py``.


.. _examples-exact-cfbend:

Symplectic Integration in an Exact Combined-Function Bend
=========================================================

This benchmark tests the use of the ``ExactCFbend`` (``cfbend_exact``) element for integrating through a combined-function bend using the exact nonlinear Hamiltonian.

This example tests the transport of a 2 GeV electron bunch through a combined function bending element with only the lowest-order (dipole) coefficient nonzero, representing the effect
of a pure dipole field.  This is compared with the effect of the ``ExactSbend`` element for equivalent field strength, by applying the inverse map (ds < 0 and phi < 0).

As a result, the second moments of x, y, and t and the associated emittances of the bunch (as well as individual particle coordinates) should all be exactly unchanged.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with
nominal values, to within the specified tolerance.

Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_exact_cfbend.py`` or
* ImpactX **executable** using an input file: ``impactx input_exact_cfbend.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_exact_cfbend.py
          :language: python3
          :caption: You can copy this file from ``examples/symplectic_integration/run_exact_cfbend.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_exact_cfbend.in
          :language: ini
          :caption: You can copy this file from ``examples/symplectic_integration/input_exact_cfbend.in``.
