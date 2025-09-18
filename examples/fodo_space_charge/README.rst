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


.. _examples-fodo-envelope-sc-gaussian:

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
