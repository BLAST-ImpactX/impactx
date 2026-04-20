.. _examples-simple-booster:

Fermilab Booster (Simple Lattice)
=================================

This example models the `Fermilab Booster <https://www.fnal.gov/pub/science/particle-accelerators/accelerator-complex.html>`__ synchrotron as an ImpactX lattice.
The lattice is expressed as a Python fragment that instantiates the equivalent ImpactX elements from a translated MAD-X description.

.. note::

   TODO for Eric: a short physics description of the lattice and what this example demonstrates.


Run
---

This example can be run as:

* **Python** script: ``python3 run_simple_booster.py``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

      .. dropdown:: Booster Lattice: Python List
         :color: light
         :icon: info
         :animate: fade-in-slide-down

         .. literalinclude:: booster_impactx_lattice.py
            :language: python3
            :caption: You can copy this file from ``examples/simple_booster/booster_impactx_lattice.py``.

      .. literalinclude:: run_simple_booster.py
         :language: python3
         :caption: You can copy this file from ``examples/simple_booster/run_simple_booster.py``.

   .. tab-item:: Reference: MAD-X Lattice

      The original MAD-X description is kept for reference.

      .. dropdown:: Booster Lattice: MAD-X File
         :color: light
         :icon: info
         :animate: fade-in-slide-down

         .. literalinclude:: sbbooster-cooked-rfon.madx
            :language: text
            :caption: You can copy this file from ``examples/simple_booster/sbbooster-cooked-rfon.madx``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_simple_booster.py``

   .. literalinclude:: analysis_simple_booster.py
      :language: python3
      :caption: You can copy this file from ``examples/simple_booster/analysis_simple_booster.py``.


Visualize
---------

We run the following script to plot the lattice survey:

.. dropdown:: Script ``plot_simple_booster.py``

   .. literalinclude:: plot_simple_booster.py
      :language: python3
      :caption: You can copy this file from ``examples/simple_booster/plot_simple_booster.py``.
