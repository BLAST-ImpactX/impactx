.. _examples-dogleg:

Dogleg
======

This is a 2-bend dogleg lattice obtained by taking the first 1/2 of the Berlin-Zeuthen magnetic bunch compression chicane:
https://www.desy.de/csr/

The primary purpose is to benchmark the reduced beam diagnostics in lattice regions with nonzero dispersion.

All parameters can be found online.  A 5 GeV electron bunch with normalized transverse rms emittance of 1 um is used.

The final expected dispersion is 267 mm.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must
agree with nominal values.

In addition, the initial and final values of :math:`\alpha_x`, :math:`\alpha_y`, :math:`\beta_x`, :math:`\beta_y`, :math:`\dispersion_x`, and :math:`\dispersion_px` must
agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_dogleg.py`` or
* ImpactX **executable** using an input file: ``impactx input_dogleg.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_dogleg.py
          :language: python3
          :caption: You can copy this file from ``examples/dogleg/run_dogleg.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_dogleg.in
          :language: ini
          :caption: You can copy this file from ``examples/dogleg/input_dogleg.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_dogleg.py``

   .. literalinclude:: analysis_dogleg.py
      :language: python3
      :caption: You can copy this file from ``examples/dogleg/analysis_dogleg.py``.


.. _examples-dogleg-reverse:

Dogleg in Reverse
=================

This is the reverse of the 2-bend dogleg lattice, obtained by taking the second 1/2 of the Berlin-Zeuthen magnetic bunch compression chicane:
https://www.desy.de/csr/

The primary purpose is to demonstrate the initialization of a beam with a nonzero x-pt correlation in a lattice region with nonzero dispersion.

In this example, the x-pt correlation of the beam is removed at the dogleg exit, and the dispersion goes to zero.

The initial dispersion is taken to be -267 mm.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must
agree with nominal values.

In addition, the initial and final values of :math:`\alpha_x`, :math:`\alpha_y`, :math:`\beta_x`, :math:`\beta_y`, :math:`\dispersion_x`, and :math:`\dispersion_px` must
agree with nominal values.

Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_dogleg_reverse.py`` or
* ImpactX **executable** using an input file: ``impactx input_dogleg_reverse.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_dogleg_reverse.py
          :language: python3
          :caption: You can copy this file from ``examples/dogleg/run_dogleg_reverse.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_dogleg_reverse.in
          :language: ini
          :caption: You can copy this file from ``examples/dogleg/input_dogleg_reverse.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_dogleg_reverse.py``

   .. literalinclude:: analysis_dogleg_reverse.py
      :language: python3
      :caption: You can copy this file from ``examples/dogleg/analysis_dogleg_reverse.py``.


.. _examples-dogleg-jitter:

Dogleg with Energy Jitter
=========================

This is identical to the :ref:`dogleg example <examples-dogleg>`, except the initial beam distribution has been given a 2.5% offset in the value of mean energy.

The primary purpose is to demonstrate the use of a beam centroid offset to study the effects of, e.g. energy jitter.

The 2.5% energy offset couples through the lattice R16 (dispersion) to result in a mean horizontal x-offset at the end of the dogleg.


In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must
agree with nominal values.

In addition, the initial and final values of :math:`\alpha_x`, :math:`\alpha_y`, :math:`\beta_x`, :math:`\beta_y`, :math:`\dispersion_x`, and :math:`\dispersion_px` must
agree with nominal values.

Finally, the values of :math:`\mean_pt`, :math:`\mean_x`, and :math:`\mean_px` must agree with predicted values.

Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_dogleg_jitter.py`` or
* ImpactX **executable** using an input file: ``impactx input_dogleg_jitter.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_dogleg_jitter.py
          :language: python3
          :caption: You can copy this file from ``examples/dogleg/run_dogleg_jitter.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_dogleg_jitter.in
          :language: ini
          :caption: You can copy this file from ``examples/dogleg/input_dogleg_jitter.in``.

Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_dogleg_jitter.py``

   .. literalinclude:: analysis_dogleg_jitter.py
      :language: python3
      :caption: You can copy this file from ``examples/dogleg/analysis_dogleg_jitter.py``.
