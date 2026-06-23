.. _examples-linac-segment:

PIP-II Linac Segment
====================

This example illustrates a small segment of the PIP-II linac, which includes the last RF cavity in the MEBT section (a quarter-wave resonator cavity for longitudinal bunching) and the first
RF cavity in the succeeding cryomodule (a superconducting half-wave resonator for acceleration).

For simplicity, aperture restrictions have been removed, and there are no active kickers.  RF cavity and solenoid field maps are represented using the Fourier coefficients of the longitudinal fields on-axis.

The nominal beam current is 4.845 mA (29.83 pC at 162.5 MHz), which assumes ~3% beam loss relative to the 5 mA current entering the MEBT.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_linac_segment.py`` or
* ImpactX **executable** using an input file: ``impactx input_linac_segment.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_linac_segment.py
          :language: python3
          :caption: You can copy this file from ``examples/linac_segment/run_linac_segment.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_linac_segment.in
          :language: ini
          :caption: You can copy this file from ``examples/linac_segment/input_linac_segment.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_linac_segment.py``

   .. literalinclude:: analysis_linac_segment.py
      :language: python3
      :caption: You can copy this file from ``examples/linac_segment/analysis_linac_segment.py``.
