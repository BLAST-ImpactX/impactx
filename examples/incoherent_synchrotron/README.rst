.. _examples-bend-isr:

A Single Bend with ISR
=======================

This test illustrates the effects of incoherent synchrotron radiation (ISR) on a high-energy electron bunch in a bending dipole.  The effects of CSR are turned off.

An initially cold (zero emittance), 100 GeV electron bunch with a 0.1 mm rms beam size propagates in an ideal sector bend, through a region of uniform magnetic field.

Due to the stochastic emission of synchrotron radiation, the bunch experiences a mean loss of energy as well as a growth in rms energy spread.  Due to the presence of nonzero dispersion,
this also results in a growth of the bunch rms emittance.

The resulting quantities are compared against the expressions given in:

N. Yampolsky and B. Carlsten, "Beam debunching due to ISR-induced energy diffusion," Nucl. Instrum. and Methods in Phys. Res. A,
870, 156-162 (2017), equations (17-20).

In this test, the values of :math:`\langle{p_x\rangle}`, :math:`\langle{p_y\rangle}`, :math:`\langle{p_t\rangle}`, :math:`\sigma_{p_x}`, :math:`\sigma_{p_y}`, and :math:`\sigma_{p_t}` must agree with nominal values.

In addition, the final values of :math:`\langle{p_t\rangle}` and :math:`\sigma_{p_t}` must agree with predicted values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_bend_isr.py`` or
* ImpactX **executable** using an input file: ``impactx input_bend_isr.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_bend_isr.py
          :language: python3
          :caption: You can copy this file from ``examples/incoherent_synchrotron/run_bend_isr.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_bend_isr.in
          :language: ini
          :caption: You can copy this file from ``examples/incoherent_synchrotron/input_bend_isr.in``.

Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_bend_isr.py``

   .. literalinclude:: analysis_bend_isr.py
      :language: python3
      :caption: You can copy this file from ``examples/incoherent_synchrotron/analysis_bend_isr.py``.
