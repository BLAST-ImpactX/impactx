
.. _examples-kick-Gauss2p5D:

Momentum Kick in the 2.5D Gaussian Space Charge Model
======================================================

The test compares the space charge momentum kick experienced by particles within a long, initially cold, Gaussian bunch against the 2.5D analytical model in:

J. Qiang, Phys. Rev. Accel. Beams 28, 114602 (2025), eqs. (31-32).

We use a cold (zero emittance) 100 MeV electron bunch in a drift space, with a single space charge slice.

Here :math:`\sigma_x = \sigma_y`, so the predicted momentum kick can be evaluated in closed form.

In this test, the final momentum of every particle is compared against its predicted value.  Both the maximum error and its rms value over the bunch population are returned.  Both values must be zero within a small tolerance.

The test can be run using either the Gauss2p5D space charge model or the 2.5D PIC space charge model.


Run
---

This example can be run as:

* **Python** script: ``python3 run_sc_kick_Gauss2p5D.py`` or ``python3 run_sc_kick_PIC2p5D.py``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Gauss2p5D space charge solver

      .. literalinclude:: run_sc_kick_Gauss2p5D.py
         :language: python3
         :caption: You can copy this file from ``examples/space_charge_gaussian_kick/run_sc_kick_Gauss2p5D.py``.

   .. tab-item:: 2p5D (PIC) space charge solver

      .. literalinclude:: run_sc_kick_PIC2p5D.py
         :language: python3
         :caption: You can copy this file from ``examples/space_charge_gaussian_kick/run_sc_kick_PIC2p5D.py``.

Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_sc_kick_Gauss2p5D.py``

   .. literalinclude:: analysis_sc_kick_Gauss2p5D.py
      :language: python3
      :caption: You can copy this file from ``examples/space_charge_gaussian_kick/analysis_sc_kick_Gauss2p5D.py``.

