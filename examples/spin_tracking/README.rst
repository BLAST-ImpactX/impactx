.. _examples-quad-spin:

Spin Depolarization in a Quadrupole
===================================

This example illustrates the decay of the polarization vector (describing the mean of the three spin components) along the vertical y and longitudinal z directions for a beam undergoing
horizontal focusing in a quadrupole.

We use a 250 MeV proton beam with initial unnormalized rms emittance of 1 micron in the horizontal plane, and 2 micron in the vertical plane.

The beam propagates over one horizontal betatron period, to a location where the polarization vector is described by a simple expression.

In this test, the initial and final values of :math:`\langle{s_x\rangle}`, :math:`\langle{s_y\rangle}`, :math:`\langle{s_z\rangle}` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_quad_spin.py`` or
* ImpactX **executable** using an input file: ``impactx input_quad_spin.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_quad_spin.py
          :language: python3
          :caption: You can copy this file from ``examples/spin_tracking/run_quad_spin.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_quad_spin.in
          :language: ini
          :caption: You can copy this file from ``examples/spin_tracking/input_quad_spin.in``.


Analyze
-------

The analysis can be run using **either** of the following scripts:


.. dropdown:: Script ``analysis_quad_spin.py``

   .. literalinclude:: analysis_quad_spin.py
      :language: python3
      :caption: You can copy this file from ``examples/spin_tracking/analysis_quad_spin.py``.

.. dropdown:: Script ``analysis_quad_spin_rbc.py``

   .. literalinclude:: analysis_quad_spin_rbc.py
      :language: python3
      :caption: You can copy this file from ``examples/spin_tracking/analysis_quad_spin_rbc.py``.


.. _examples-fodo-spin:

Spin Depolarization in a FODO Channel
=====================================

This example illustrates the decay of the polarization vector (describing the mean of the three spin components) for a matched beam in a FODO channel.

We use a 10 GeV electron beam, with an initially 6D Gaussian distribution in the phase space.  The FODO channel and Twiss parameters are otherwise identical to
those appearing in ``examples/fodo_channel``.

The beam spin undergoes rapid mixing and depolarization.

In this test, the initial and final values of :math:`\langle{s_x\rangle}`, :math:`\langle{s_y\rangle}`, :math:`\langle{s_z\rangle}` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_fodo_channel_spin.py`` or
* ImpactX **executable** using an input file: ``impactx input_fodo_channel_spin.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_fodo_channel_spin.py
          :language: python3
          :caption: You can copy this file from ``examples/spin_tracking/run_fodo_channel_spin.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_fodo_channel_spin.in
          :language: ini
          :caption: You can copy this file from ``examples/spin_tracking/input_fodo_channel_spin.in``.


.. _examples-sbend-spin:

Spin Depolarization in a Dipole
===============================

This example illustrates the decay of the polarization vector (describing the mean of the three spin components) along the horizontal x and longitudinal z directions for a beam undergoing
bending in the x-z plane in a sector dipole.

We use a 2 GeV electron beam.  The beam parameters (in particular, the momentum and energy spread) are artificially large in order to enhance the effect.

The beam propagates over one period, as set by the design spin tune.  By increasing the number of slices, and turning on diagnostics, one can view precession of the polarization vector about
the vertical direction.  A clean precession over 1 period becomes visible when the beam size, momentum spread, and energy spread are set to small values.

In this test, the initial and final values of :math:`\langle{s_x\rangle}`, :math:`\langle{s_y\rangle}`, :math:`\langle{s_z\rangle}` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_sbend_spin.py`` or
* ImpactX **executable** using an input file: ``impactx input_sbend_spin.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_sbend_spin.py
          :language: python3
          :caption: You can copy this file from ``examples/spin_tracking/run_sbend_spin.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_sbend_spin.in
          :language: ini
          :caption: You can copy this file from ``examples/spin_tracking/input_sbend_spin.in``.



Analyze
-------

The analysis can be run using:


.. dropdown:: Script ``analysis_sbend_spin.py``

   .. literalinclude:: analysis_sbend_spin.py
      :language: python3
      :caption: You can copy this file from ``examples/spin_tracking/analysis_sbend_spin.py``.


.. _examples-sol-spin:

Spin Depolarization in a Solenoid
=================================

This example illustrates the decay of the polarization vector (describing the mean of the three spin components) along the horizontal x and vertical y directions for a 2 GeV proton beam undergoing
transverse focusing and rotation in a solenoid.

The beam propagates over a distance:  :math:`\Delta s = 2\pi/(Gk_s)`,

where :math:`G` is the value of the gyromagnetic anomaly, :math:`k_s` denotes the solenoid focusing strength.
At this location, the polarization vector is described by a simple expression.

In this test, the initial and final values of :math:`\langle{s_x\rangle}`, :math:`\langle{s_y\rangle}`, :math:`\langle{s_z\rangle}` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_sol_spin.py`` or
* ImpactX **executable** using an input file: ``impactx input_sol_spin.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_sol_spin.py
          :language: python3
          :caption: You can copy this file from ``examples/spin_tracking/run_sol_spin.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_sol_spin.in
          :language: ini
          :caption: You can copy this file from ``examples/spin_tracking/input_sol_spin.in``.

Analyze
-------

The analysis can be run using:


.. dropdown:: Script ``analysis_sol_spin.py``

   .. literalinclude:: analysis_sol_spin.py
      :language: python3
      :caption: You can copy this file from ``examples/spin_tracking/analysis_sol_spin.py``.



.. _examples-reversibility-spin:

Element Reversibility with Spin
===============================

In the case of linear elements, including spin, the joint spin-orbit map has the following exact reversibility property.

The effect of setting ds -> -ds is equivalent to replacing the map by its inverse.

In this test, a beam is propagated forward through an element of length ds, following by the corresponding element with length -ds.  As a result, the composite map is the identity.  This provides a non-trivial test for consistency between the spin map and the orbit map.

In this test, the initial and final spin components :math:`s_x`, :math:`s_y`, and :math:`s_z` are compared.  The norm of the change in the spin vector must lie within a very small tolerance.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_reversibility_spin.py`` or
* ImpactX **executable** using an input file: ``impactx input_reversibility_spin.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_reversibility_spin.py
          :language: python3
          :caption: You can copy this file from ``examples/spin_tracking/run_reversibility_spin.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_reversibility_spin.in
          :language: ini
          :caption: You can copy this file from ``examples/spin_tracking/input_reversibility_spin.in``.

Analyze
-------

The analysis can be run using:


.. dropdown:: Script ``analysis_reversibility_spin.py``

   .. literalinclude:: analysis_reversibility_spin.py
      :language: python3
      :caption: You can copy this file from ``examples/spin_tracking/analysis_reversibility_spin.py``.


.. _examples-cfbend-spin:

Spin Limiting Cases of a Combined-Function Bend
===============================================

This example tests the spin dynamics of a combined-function bend in the quad-like and bend-like limiting cases when the curvature or the quadrupole gradient are small, respectively.
The beam parameters are identical to those used in examples-cfbend, with the addition of spin.

In this test, the beam is tracked through a combined-function dipole, and then tracked backward through an element with equivalent parameters (up to a small tolerance).
Due to the s-reversibility of the maps applied, the final phase space coordinates and spin variables should coincide with their initial values.

This test computes the norm of the change in the spin polarization vector, which must be zero to within a small tolerance.

Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_cfbend_spin.py`` or
* ImpactX **executable** using an input file: ``impactx input_cfbend_spin.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_cfbend_spin.py
          :language: python3
          :caption: You can copy this file from ``examples/spin_tracking/run_cfbend_spin.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_cfbend_spin.in
          :language: ini
          :caption: You can copy this file from ``examples/spin_tracking/input_cfbend_spin.in``.

Analyze
-------

The analysis can be run using:

.. dropdown:: Script ``analysis_cfbend_spin.py``

   .. literalinclude:: analysis_cfbend_spin.py
      :language: python3
      :caption: You can copy this file from ``examples/spin_tracking/analysis_cfbend_spin.py``.
