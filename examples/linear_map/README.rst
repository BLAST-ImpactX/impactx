.. _examples-linear-map:

Iteration of a User-Defined Linear Map
======================================

This example illustrates the application of a user-defined linear map via a matrix.

Here, the linear map represents an abstract symplectic transformation of the beam in 6D phase space.
If desired, the user may interpret the matrix as the one-turn map of a storage ring or circular collider.

The (fractional) tunes (Qx, Qy, Qt) of the map are given by (0.139, 0.219, 0.0250).
We use a 45.6 GeV electron beam that is invariant under the action of the linear map (matched).
The horizontal and vertical unnormalized emittances are 0.27 nm and 1.0 pm, respectively.

These parameters are based on the `single-beam parameters of FCC-ee (Z-mode) <https://twiki.cern.ch/twiki/bin/view/FCC/FCCeeParameters_CDRBaseline-1_0>`__.
(`backup <https://web.archive.org/web/20250000000000*/https://twiki.cern.ch/twiki/bin/view/FCC/FCCeeParameters_CDRBaseline-1_0>`__).

The second moments of the phase space variables should be unchanged under application of the map.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.

In addition, the tunes associated with a single particle orbit are extracted, and must agree with the values given above.

Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_map.py`` or
* ImpactX **executable** using an input file: ``impactx input_map.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_map.py
          :language: python3
          :caption: You can copy this file from ``examples/linear_map/run_map.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_map.in
          :language: ini
          :caption: You can copy this file from ``examples/linear_map/input_map.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_map.py``

   .. literalinclude:: analysis_map.py
      :language: python3
      :caption: You can copy this file from ``examples/linear_map/analysis_map.py``.
