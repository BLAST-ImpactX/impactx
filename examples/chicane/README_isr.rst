.. _examples-chicane-isr:

Chicane with ISR
================

This is the :ref:`Berlin-Zeuthen magnetic bunch compression chicane <examples-chicane>` example, but this time with incoherent synchrotron radiation (CSR) modelled in the bending magnets.
Coherent synchrotron radiation (CSR) is turned off.

`All parameters can be found online <https://www.desy.de/csr/>`__.
A 5 GeV electron bunch with normalized transverse rms emittance of 1 um undergoes longitudinal compression by a factor of 10 in a standard 4-bend chicane.
The rms pulse length should decrease by the compression factor (10).

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_chicane_isr.py`` or
* ImpactX **executable** using an input file: ``impactx input_chicane_isr.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_chicane_isr.py
          :language: python3
          :caption: You can copy this file from ``examples/chicane/run_chicane_isr.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_chicane_isr.in
          :language: ini
          :caption: You can copy this file from ``examples/chicane/input_chicane_isr.in``.

Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_chicane_isr.py``

   .. literalinclude:: analysis_chicane_isr.py
      :language: python3
      :caption: You can copy this file from ``examples/chicane/analysis_chicane_isr.py``.
