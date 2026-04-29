.. _examples-quadrupole-softedge:

Soft-Edge Quadrupole
====================

This is a modification of :ref:`the test for a matched electron beam propagating through a stable FODO cell <examples-fodo>`,
in which the quadrupoles have been replaced with soft-edge quadrupole
elements.  The on-axis field profile in this example has been chosen to correspond
to the hard-edge limit, so the two tests should coincide.

We use a 2 GeV electron beam with initial unnormalized rms emittance of 2 nm.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_quadrupole_softedge.py`` or
* ImpactX **executable** using an input file: ``impactx input_quadrupole_softedge.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_quadrupole_softedge.py
          :language: python3
          :caption: You can copy this file from ``examples/quadrupole_softedge/run_quadrupole_softedge.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_quadrupole_softedge.in
          :language: ini
          :caption: You can copy this file from ``examples/quadrupole_softedge/input_quadrupole_softedge.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_quadrupole_softedge.py``

   .. literalinclude:: analysis_quadrupole_softedge.py
      :language: python3
      :caption: You can copy this file from ``examples/quadrupole_softedge/analysis_quadrupole_softedge.py``.


.. _examples-fodo-softedge:

FODO Cell Using Soft-Edge Quadrupoles
=======================================

This is a modification of the example above, in which the quadrupoles have been replaced by soft-edge quadrupoles with an adjustable magnetic gap parameter,
which allows the user to tune the rate of fringe field decay near the magnetic edge.  In the limit :math:`g\rightarrow 0`, this model approaches the model of an ideal (hard-edge) quadrupole.

The on-axis gradient of the quadrupoles is given by an analytical function, namely:

:math:`\frac{k_0}{2}\left\{\tanh\left(\frac{z+z_0}{g}\right)-\tanh\left(\frac{z-z_0}{g}\right)\right\}`

The effect of the soft fringe field decay on the linear focusing is small but visible.  In this test, the beam moments are compared against those used in the hard-edge limit.

The initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run as:

* **Python** script: ``python3 run_fodo_softedge.py``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_fodo_softedge.py
          :language: python3
          :caption: You can copy this file from ``examples/quadrupole_softedge/run_fodo_softedge.py``.

Analyze

-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_fodo_softedge.py``

   .. literalinclude:: analysis_fodo_softedge.py
      :language: python3
      :caption: You can copy this file from ``examples/quadrupole_softedge/analysis_fodo_softedge.py``.
