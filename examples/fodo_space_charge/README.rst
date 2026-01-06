.. _examples-fodo-envelope-sc:

FODO Cell with 2D Space Charge using Envelope Tracking
======================================================

This example illustrates a 0.5 A proton beam with a kinetic energy of 6.7 MeV in a FODO cell,
with 2D space charge included.  The parameters are those described in:

R.D. Ryne et al, `"A Test Suite of Space-Charge Problems for Code Benchmarking," <https://accelconf.web.cern.ch/e04/PAPERS/WEPLT047.PDF>`__
in Proc. EPAC 2004, Lucerne, Switzerland:  KV Beam in a FODO Channel

The purpose of this example is to illustrate the use of envelope tracking mode with 2D space charge.

The second moments of the particle distribution after the FODO cell should coincide with the second moments of the particle distribution before the FODO cell, to within the level expected due to noise due to statistical sampling.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_fodo_envelope_sc.py`` or
* ImpactX **executable** using an input file: ``impactx input_fodo_envelope_sc.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_fodo_envelope_sc.py
          :language: python3
          :caption: You can copy this file from ``examples/fodo_space_charge/run_fodo_envelope_sc.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_fodo_envelope_sc.in
          :language: ini
          :caption: You can copy this file from ``examples/fodo_space_charge/input_fodo_envelope_sc.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_fodo_envelope_sc.py``

   .. literalinclude:: analysis_fodo_envelope_sc.py
      :language: python3
      :caption: You can copy this file from ``examples/fodo_space_charge/analysis_fodo_envelope_sc.py``.


.. _examples-fodo-Gaussian-sc:

FODO Cell with 3D Gaussian Space Charge Using Particle Tracking
===============================================================

This example illustrates a 1 nC electron beam with a kinetic energy of 100 MeV in a FODO cell,
with 3D Gaussian space charge included.  The parameters are those described in:

The purpose of this example is to illustrate the use of particle tracking mode with 3D space charge from a Gaussian density
distribution.

The second moments of the particle distribution after the FODO cell should coincide with the second moments of the particle distribution before the FODO cell with small noticeable growth.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_fodo_Gauss3D_sc.py`` or
* ImpactX **executable** using an input file: ``impactx input_fodo_Gauss3D_sc.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_fodo_Gauss3D_sc.py
          :language: python3
          :caption: You can copy this file from ``examples/fodo_space_charge/run_fodo_Gauss3D_sc.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_fodo_Gauss3D_sc.in
          :language: ini
          :caption: You can copy this file from ``examples/fodo_space_charge/input_fodo_Gauss3D_sc.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_fodo_Gauss3D_sc.py``

   .. literalinclude:: analysis_fodo_Gauss3D_sc.py
      :language: python3
      :caption: You can copy this file from ``examples/fodo_space_charge/analysis_fodo_Gauss3D_sc.py``.


.. _examples-fodo-2p5dGaussian-sc:

FODO Cell with 2.5D Gaussian Space Charge Using Particle Tracking
=================================================================

This example illustrates a 1 nC electron beam with a kinetic energy of 100 MeV in a FODO cell,
with 2.5D Gaussian space charge included.

The purpose of this example is to illustrate the use of particle tracking mode with 2.5D space charge from a Gaussian density distribution.

The second moments of the particle distribution after the FODO cell should coincide with the second moments of the particle distribution before the FODO cell with small noticeable growth.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.

Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_fodo_Gauss2p5D_sc.py`` or
* ImpactX **executable** using an input file: ``impactx input_fodo_Gauss2p5D_sc.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_fodo_Gauss2p5D_sc.py
          :language: python3
          :caption: You can copy this file from ``examples/fodo_space_charge/run_fodo_Gauss2p5D_sc.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_fodo_Gauss2p5D_sc.in
          :language: ini
          :caption: You can copy this file from ``examples/fodo_space_charge/input_fodo_Gauss2p5D_sc.in``.

Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_fodo_Gauss2p5D_sc.py``

   .. literalinclude:: analysis_fodo_Gauss2p5D_sc.py
      :language: python3
      :caption: You can copy this file from ``examples/fodo_space_charge/analysis_fodo_Gauss2p5D_sc.py``.


.. _examples-fodo-2d-sc:

FODO Cell with 2D Space Charge using FFT IGF Poisson Solver
===========================================================

This example illustrates a 0.5 A proton beam with a kinetic energy of 6.7 MeV in a FODO cell,
with 2D space charge included.  The parameters are those described in:

R.D. Ryne et al, `"A Test Suite of Space-Charge Problems for Code Benchmarking," <https://accelconf.web.cern.ch/e04/PAPERS/WEPLT047.PDF>`__
in Proc. EPAC 2004, Lucerne, Switzerland:  KV Beam in a FODO Channel

The purpose of this example is to illustrate the use of partilce tracking mode with 2D space charge.

The second moments of the particle distribution after the FODO cell should coincide with the second moments of the particle distribution before the FODO
cell, to within the specified tolerance.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and
:math:`\epsilon_t` must agree with nominal values.

Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_fodo_2d_sc.py`` or
* ImpactX **executable** using an input file: ``impactx input_fodo_2d_sc.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_fodo_2d_sc.py
          :language: python3
          :caption: You can copy this file from ``examples/fodo_space_charge/run_fodo_2d_sc.py``.

   .. tab-item:: Executable: Input File
       .. literalinclude:: input_fodo_2d_sc.in
          :language: ini
          :caption: You can copy this file from ``examples/fodo_space_charge/input_fodo_2d_sc.in``.

Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_fodo_2d_sc.py``

   .. literalinclude:: analysis_fodo_2d_sc.py
      :language: python3
      :caption: You can copy this file from ``examples/fodo_space_charge/analysis_fodo_2d_sc.py``.


.. _examples-fodo-2p5d-sc:

FODO Cell with 2.5D Space Charge using FFT IGF Poisson Solver
=============================================================

This example illustrates a 1 nC electron beam with a kinetic energy of 100 MeV in a FODO cell, with 2.5D space charge included.  The problem setup is identical
to ``examples-fodo-2p5dGaussian-sc`` above, except for the choice of space charge model.

The purpose of this example is to illustrate the use of particle tracking mode with 2.5D space charge, obtained using a particle-in-cell method based on a 2D Poisson solve with
transverse space charge weighted by the local value of beam current.

The second moments of the particle distribution after the FODO cell should coincide with the second moments of the particle distribution before the FODO cell with small noticeable growth.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.

Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_fodo_2p5D_sc.py`` or
* ImpactX **executable** using an input file: ``impactx input_fodo_2p5D_sc.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_fodo_2p5D_sc.py
          :language: python3
          :caption: You can copy this file from ``examples/fodo_space_charge/run_fodo_2p5D_sc.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_fodo_2p5D_sc.in
          :language: ini
          :caption: You can copy this file from ``examples/fodo_space_charge/input_fodo_2p5D_sc.in``.

Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_fodo_2p5D_sc.py``

   .. literalinclude:: analysis_fodo_2p5D_sc.py
      :language: python3
      :caption: You can copy this file from ``examples/fodo_space_charge/analysis_fodo_2p5D_sc.py``.


.. _examples-fodo-2p5d-kv-sc:

FODO Cell with 2.5D Space Charge using FFT IGF Poisson Solver (Matched)
=======================================================================

This example illustrates an 0.5 A proton beam with a kinetic energy of 6.7 MeV in a FODO cell, with 2.5D space charge included.  A KV beam distribution is used.  The problem setup is identical
to ``examples-fodo-2d-sc`` above, except for the following features:

1) The 2D FFT IGF Poisson solver is replaced by the 2.5D FFT IGF Poisson solver.
2) A long bunch with uniform current profile is used, and the total bunch charge is chosen to produce 0.5 A current in the bunch.
3) The energy spread of the bunch is set to zero, and the longitudinal space charge kick is turned off, in order to avoid debunching and retain 0.5 A beam current.

The purpose of this example is to test for consistency of the 2.5D solver with the 2D solver, in the case of a long, uniform bunch.

The second moments of the particle distribution after the FODO cell should coincide with the second moments of the particle distribution before the FODO cell.  (The beam is matched.)

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.

Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_fodo_2p5D_KV_sc.py`` or
* ImpactX **executable** using an input file: ``impactx input_fodo_2p5D_KV_sc.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_fodo_2p5D_KV_sc.py
          :language: python3
          :caption: You can copy this file from ``examples/fodo_space_charge/run_fodo_2p5D_KV_sc.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_fodo_2p5D_KV_sc.in
          :language: ini
          :caption: You can copy this file from ``examples/fodo_space_charge/input_fodo_2p5D_KV_sc.in``.

Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_fodo_2p5D_KV_sc.py``

   .. literalinclude:: analysis_fodo_2p5D_KV_sc.py
      :language: python3
      :caption: You can copy this file from ``examples/fodo_space_charge/analysis_fodo_2p5D_KV_sc.py``.
