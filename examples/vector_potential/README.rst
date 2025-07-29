.. _examples-fodo-vector-potential:

Symplectic Integration in a FODO Cell Based on User-Provided Vector Potential
=============================================================================

This benchmark tests the implementation of the element MagnetostaticVectorPotential, which uses symplectic integration based on the exact nonlinear Hamiltonian for a user-provided vector 
potential.  The physical example is identical to examples-fodo, except that the vector potential in the quadrupoles is provided in formulae by the user.

The analysis script is identical to the analysis script used for ``examples-fodo``.

Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_fodo_vector_potential.py`` or
* ImpactX **executable** using an input file: ``impactx input_fodo_vector_potential.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_fodo_vector_potential.py
          :language: python3
          :caption: You can copy this file from ``examples/vector_potential/run_fodo_vector_potential.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_fodo_vector_potential.in
          :language: ini
          :caption: You can copy this file from ``examples/vector_potential/input_fodo_vector_potential.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_fodo_vector_potential.py``

   .. literalinclude:: analysis_fodo_vector_potential.py
      :language: python3
      :caption: You can copy this file from ``examples/vector_potential/analysis_vector_potential.py``.


.. _examples-exact-quad-vector-potential:

Symplectic Integration in Exact Quads Based on User-Provided Vector Potential
=============================================================================

This benchmark tests the implementation of the element MagnetostaticVectorPotential, which uses symplectic integration based on the exact nonlinear Hamiltonian for a user-provided vector 
potential.  The physical example is identical to examples-exact-quad, except that the vector potential in the quadrupoles is provided in formulae by the user.

The analysis script is identical to the analysis script used for ``examples-exact-quad``.

Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_exact_quad_vector_potential.py`` or
* ImpactX **executable** using an input file: ``impactx input_exact_quad_vector_potential.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_exact_quad_vector_potential.py
          :language: python3
          :caption: You can copy this file from ``examples/vector_potential/run_exact_quad_vector_potential.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_exact_quad_vector_potential.in
          :language: ini
          :caption: You can copy this file from ``examples/vector_potential/input_exact_quad_vector_potential.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_exact_quad_vector_potential.py``

   .. literalinclude:: analysis_exact_quad_vector_potential.py
      :language: python3
      :caption: You can copy this file from ``examples/vector_potential/analysis_exact_quad_vector_potential.py``.



.. _examples-solenoid-vector-potential:

Symplectic Integration Through a Soft-Edge Solenoid Based on User-Provided Vector Potential
============================================================================================

This benchmark tests the implementation of the element MagnetostaticVectorPotential, which uses symplectic integration based on the exact nonlinear Hamiltonian for a user-provided vector 
potential.  The physical example is identical to examples-solenoid-softedge, except that the vector potential of the solenoid is provided in formulae by the user.

The analysis script is identical to the analysis script used for ``examples-solenoid-softedge``.

Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_solenoid_vector_potential.py`` or
* ImpactX **executable** using an input file: ``impactx input_solenoid_vector_potential.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_solenoid_vector_potential.py
          :language: python3
          :caption: You can copy this file from ``examples/vector_potential/run_solenoid_vector_potential.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_solenoid_vector_potential.in
          :language: ini
          :caption: You can copy this file from ``examples/vector_potential/input_solenoid_vector_potential.in``.
