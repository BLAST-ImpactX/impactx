.. _examples-htu-beamline:

Dynamics in the HTU Beamline
============================

This example involves tracking a 25 pC electron bunch with 100 MeV total energy through the `BELLA Hundred-Terawatt Undulator (HTU) beamline at LBNL <https://doi.org/10.1117/12.3056776>`__.

The bunch is generated in practice from a laser-plasma accelerator stage, resulting in a small intial rms beam size (3.9, 3.9, 1.0) microns, 2 mrad transverse divergence and 2.5% energy spread.  Due
to the large energy spread, chromatic focusing effects are important.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run as:

* **Python** script: ``python3 run_impactx.py`` or

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_impactx.py
          :language: python3
          :caption: You can copy this file from ``examples/htu_beamline/run_impactx.py``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_htu_beamline.py``

   .. literalinclude:: analysis_htu_beamline.py
      :language: python3
      :caption: You can copy this file from ``examples/htu_beamline/analysis_htu_beamline.py``.
