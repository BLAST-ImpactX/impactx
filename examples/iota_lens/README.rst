.. _examples-iotalens:

A nonlinear focusing channel based on the IOTA nonlinear lens
=============================================================

A constant focusing channel with nonlinear focusing, using a string of thin
IOTA nonlinear lens elements alternating with constant focusing elements.

We use a 2.5 MeV proton beam, corresponding to the nominal IOTA proton energy.

The two functions H (Hamiltonian) and I (the second invariant) should remain unchanged for all particles.

In this test, the initial and final values of :math:`\mu_H`, :math:`\sigma_H`, :math:`\mu_I`, :math:`\sigma_I` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_iotalens.py`` or
* ImpactX **executable** using an input file: ``impactx input_iotalens.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_iotalens.py
          :language: python3
          :caption: You can copy this file from ``examples/iota_lens/run_iotalens.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_iotalens.in
          :language: ini
          :caption: You can copy this file from ``examples/iota_lens/input_iotalens.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_iotalens.py``

   .. literalinclude:: analysis_iotalens.py
      :language: python3
      :caption: You can copy this file from ``examples/iota_lens/analysis_iotalens.py``.


.. _examples-iotalens-sdep:

A nonlinear focusing channel based on the physical IOTA nonlinear magnet
========================================================================

A representation of the physical IOTA nonlinear magnetic element with realistic
s-dependence, obtained using a sequence of nonlinear lenses and drifts equivalent
to the use of a second-order symplectic integrator.

A thin linear lens is added at the exit of the nonlinear element, representing the
ideal map associated with the remainder of the lattice.

We use a 2.5 MeV proton beam, corresponding to the nominal IOTA proton energy.

The two functions H (Hamiltonian) and I (the second invariant) are evaluated at the
entrance to the nonlinear element, and then again after the thin lens (representing a
single period).  These values should be unchanged for all particles (to within acceptable tolerance).

In this test, the initial and final values of :math:`\mu_H`, :math:`\sigma_H`, :math:`\mu_I`, :math:`\sigma_I` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_iotalens_sdep.py`` or
* ImpactX **executable** using an input file: ``impactx input_iotalens_sdep.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_iotalens_sdep.py
          :language: python3
          :caption: You can copy this file from ``examples/iota_lens/run_iotalens_sdep.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_iotalens_sdep.in
          :language: ini
          :caption: You can copy this file from ``examples/iota_lens/input_iotalens_sdep.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_iotalens_sdep.py``

   .. literalinclude:: analysis_iotalens_sdep.py
      :language: python3
      :caption: You can copy this file from ``examples/iota_lens/analysis_iotalens_sdep.py``.


.. _examples-iotalens-sdep-aperture:

A nonlinear focusing channel based on the physical IOTA nonlinear magnet with aperture
======================================================================================

A representation of the physical IOTA nonlinear magnetic element with realistic
s-dependence, obtained using a sequence of nonlinear lenses and drifts equivalent
to the use of a second-order symplectic integrator.  Transverse aperture restrictions
are included.  The elliptical aperture dimensions are based on the physical magnet
design.

A thin linear lens is added at the exit of the nonlinear element, representing the
ideal map associated with the remainder of the lattice.

We use a 2.5 MeV proton beam, corresponding to the nominal IOTA proton energy.  The beam is transported over 3 periods.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.

The fraction of charge lost after 3 periods of transport must agree with the nominal value obtained from high-resolution simulation.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_iotalens_sdep_aperture.py`` or
* ImpactX **executable** using an input file: ``impactx input_iotalens_sdep_aperture.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_iotalens_sdep_aperture.py
          :language: python3
          :caption: You can copy this file from ``examples/iota_lens/run_iotalens_sdep_aperture.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_iotalens_sdep_aperture.in
          :language: ini
          :caption: You can copy this file from ``examples/iota_lens/input_iotalens_sdep_aperture.in``.
