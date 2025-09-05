.. _examples-solenoid-softedge:

Soft-edge solenoid
===================

Proton beam propagating through a 6 m region containing a soft-edge
solenoid.

The solenoid model used is the default thin-shell model described in:
P. Granum et al, "Efficient calculations of magnetic fields of solenoids for simulations,"
NIMA 1034, 166706 (2022)
`DOI:10.1016/j.nima.2022.166706 <https://doi.org/10.1016/j.nima.2022.166706>`__

The solenoid is a cylindrical current sheet with a length of 1 m and a
radius of 0.1667 m, corresponding to an aspect ratio diameter/length = 1/3.
The peak magnetic field on-axis is 3 T.

We use a 250 MeV proton beam with initial unnormalized rms emittance of 1 micron
in the horizontal plane, and 2 micron in the vertical plane.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_solenoid_softedge.py`` or
* ImpactX **executable** using an input file: ``impactx input_solenoid_softedge.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_solenoid_softedge.py
          :language: python3
          :caption: You can copy this file from ``examples/solenoid_softedge/run_solenoid_softedge.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_solenoid_softedge.in
          :language: ini
          :caption: You can copy this file from ``examples/solenoid_softedge/input_solenoid_softedge.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_solenoid_softedge.py``

   .. literalinclude:: analysis_solenoid_softedge.py
      :language: python3
      :caption: You can copy this file from ``examples/solenoid_softedge/analysis_solenoid_softedge.py``.


.. _solenoid-softedge-solvable:

Exactly-solvable (non-uniform) soft-edge solenoid
==================================================

This benchmark checks the calculation of the linear map for a soft-edge solenoid with a non-uniform longitudinal
on-axis magnetic field profile, in the special case of a field profile for which the map is exactly-solvable.

The test involves 250 MeV protons propagating through a 2 m region with a solenoidal magnetic field.

The user can run this test for various values of the magnetic field strength by modifying the parameter "bscale" in both the input file and in the analysis script.  (The values in the two files must be consistent.)

In this test, all 36 elements of the 6x6 transport matrix must agree with predicted values to within numerical tolerance (currently 1e-7).


Run
---

This example can be run as:

* **Python** script: ``python3 run_solenoid_softedge.py``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_solenoid_softedge.py
          :language: python3
          :caption: You can copy this file from ``examples/solenoid_softedge/run_solenoid_softedge_solvable.py``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_solenoid_softedge_solvable.py``

   .. literalinclude:: analysis_solenoid_softedge_solvable.py
      :language: python3
      :caption: You can copy this file from ``examples/solenoid_softedge/analysis_solenoid_softedge_solvable.py``.
