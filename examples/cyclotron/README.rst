.. _examples-cyclotron:

Cyclotron
=========

This demonstrates a simple cyclotron as published by Ernest O. Lawrence and M. Stanley Livingston, *The Production of High Speed Light Ions Without the Use of High Voltages*, Phys. Rev. **40**, 19 (1932).
`DOI: 10.1103/PhysRev.40.19 <https://doi.org/10.1103/PhysRev.40.19>`__


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_cyclotron.py`` or
* ImpactX **executable** using an input file: ``impactx input_cyclotron.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_cyclotron.py
          :language: python3
          :caption: You can copy this file from ``examples/cyclotron/run_cyclotron.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_cyclotron.in
          :language: ini
          :caption: You can copy this file from ``examples/cyclotron/input_cyclotron.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_cyclotron.py``

   .. literalinclude:: analysis_cyclotron.py
      :language: python3
      :caption: You can copy this file from ``examples/cyclotron/analysis_cyclotron.py``.


Visualize
---------

.. note::

   TODO :)


.. _examples-cyclotron-loss:

Cyclotron with Dynamical Losses
===============================

This is identical to the example above, except that the initial energy spread of the beam is increased.  The value is so large that some particles are not successfully transported.
These particles are considered "lost".

In this test, the initial values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.

In addition, the fraction of charge that is lost must agree with the nominal value of 7.7754%, to within a specified tolerance.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_cyclotron_loss.py`` or
* ImpactX **executable** using an input file: ``impactx input_cyclotron_loss.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_cyclotron_loss.py
          :language: python3
          :caption: You can copy this file from ``examples/cyclotron/run_cyclotron_loss.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_cyclotron_loss.in
          :language: ini
          :caption: You can copy this file from ``examples/cyclotron/input_cyclotron_loss.in``.

Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_cyclotron_loss.py``

   .. literalinclude:: analysis_cyclotron_loss.py
      :language: python3
      :caption: You can copy this file from ``examples/cyclotron/analysis_cyclotron_loss.py``.
