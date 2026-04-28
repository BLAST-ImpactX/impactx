.. _examples-fodo:

Simple Booster
==============

Simplified model of the Fermilab Booster with 24 cells and 22 RF cavities.
The initial total voltage at injection is 200 KV and the Booster operated
a harmonic number of 84.

The Booster is a 474.202 m long rapid cycling synchrotron built in 1969-71
designed to rapidly accelerate
protons from a kinetic energy of 0.8 GeV to 8 GeV although this example
just runs at the injection energy.
The architecture of the accelerator is a modified FODO design and is usually
described as "FOFDOOD" which results from splitting the F and D magnets in two,
then elongating the space between the two D magnets.
The F and D magnets combine both bending and
focussing(defocussing) functions to economomize on space which is quite limited.

The lattice is read from the file 'booster_impactx_lattice.py' which is
based on the original MAD-X file sbbooster-cooked-rfon.madx.

The matched Twiss parameters determined by both Synergia and MAD-X at entry are:

* :math:`\beta_\mathrm{x} = 33.73645362843065243` m
* :math:`\alpha_\mathrm{x} = -0.01.298673960026007664`
* :math:`\beta_\mathrm{y} =  5.252517912567207681` m
* :math:`\alpha_\mathrm{y} = 0.006089861210659328755`
* :math:`\mathrm{D}_\mathrm{x} = 3.785167992` m
* :math:`\mathrm{Dp}_\mathrm{x} = 0.001377568703`

The initial beam parameters follow the PIP-II Booster beam at injection
as described in the Conceptual Design
Report. The beam consists of protons at kinetic energy 800 MeV.

+------------------------+-----------------------------------------------+
| :math:`\epsilon_{x}`   | :math:`16 \pi \mathrm{mm-mr}` normalized 95%  |
+------------------------+-----------------------------------------------+
| :math:`\epsilon_{y}`   | :math:`16 \pi \mathrm{mm-mr}` normalized 95%  |
+------------------------+-----------------------------------------------+
| :math:`\epsilon_{L}`   | :math:`0.1 \mathrm{eV-s}` 97%                 |
+------------------------+-----------------------------------------------+



Run
---

This example can only be run with a python script:

* **Python** script: ``python3 run_simple_booster.py``


For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_simple_booster.py
          :language: python3
          :caption: You can copy this file from ``examples/simple_booster/run_simple_booster.py``. The file `booster_impactx_lattice.py` from the same directory is also required.



Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_booster_simple.py``

   .. literalinclude:: analysis_booster_simple.py
      :language: python3
      :caption: You can copy this file from ``examples/fodo/analysis_booster_simple.py``.

The second moments of the particle distribution after the FODO cell
should coincide with the second moments of the particle distribution
before the FODO cell, to within the level expected due to noise due to statistical sampling.


Visualize
---------

You can run the following script to visualize the beam evolution over time:

.. dropdown:: Script ``plot_simple_booster.py``

   .. literalinclude:: plot_simple_boooster.py
      :language: python3
      :caption: You can copy this file from ``examples/fodo/plot_simple_booster.py``.

.. figure:: https://user-images.githubusercontent.com/1353258/180287840-8561f6fd-278f-4856-abd8-04fbdb78c8ff.png
   :alt: focusing, defocusing and preserved emittance in our FODO cell benchmark.

   FODO transversal beam width and emittance evolution

.. figure:: https://user-images.githubusercontent.com/1353258/180287845-eb0210a7-2500-4aa9-844c-67fb094329d3.png
   :alt: focusing, defocusing and phase space rotation in our FODO cell benchmark.

   FODO transversal beam width and phase space evolution


.. _examples-fodo-envelope:

FODO Cell Using Envelope Tracking
=================================

This identical to the FODO example, except that envelope tracking is used instead of particle tracking.

Stable FODO cell with a zero-current phase advance of 67.8 degrees.

The matched Twiss parameters at entry are:

* :math:`\beta_\mathrm{x} = 2.82161941` m
* :math:`\alpha_\mathrm{x} = -1.59050035`
* :math:`\beta_\mathrm{y} = 2.82161941` m
* :math:`\alpha_\mathrm{y} = 1.59050035`

We use a 2 GeV electron beam with initial unnormalized rms emittance of 2 nm.

The second moments of the particle distribution after the FODO cell should coincide with the second moments of the particle distribution before the FODO cell, to within the level expected
due to$

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_fodo_envelope.py`` or
* ImpactX **executable** using an input file: ``impactx input_fodo_envelope.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_fodo_envelope.py
          :language: python3
          :caption: You can copy this file from ``examples/fodo/run_fodo_envelope.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_fodo_envelope.in
          :language: ini
          :caption: You can copy this file from ``examples/fodo/input_fodo_envelope.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_fodo_envelope.py``

   .. literalinclude:: analysis_fodo_envelope.py
      :language: python3
      :caption: You can copy this file from ``examples/fodo/analysis_fodo_envelope.py``.


.. _examples-fodo-exact:

FODO Cell Using Nonlinear Tracking
===================================

This is identical to the example ``examples-fodo``, except that fully nonlinear tracking is used based on the exact relativistic Hamiltonian.

The kinematic nonlinear effects are essentially negligible, so this is primarily a test that the nonlinear elements correctly reproduce the results of linear tracking.

The second moments of the particle distribution after the FODO cell should coincide with the second moments of the particle distribution before the FODO cell, to within the level expected due to
noise due to the finite particle population.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_fodo_exact.py`` or
* ImpactX **executable** using an input file: ``impactx input_fodo_exact.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_fodo_exact.py
          :language: python3
          :caption: You can copy this file from ``examples/fodo/run_fodo_exact.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_fodo_exact.in
          :language: ini
          :caption: You can copy this file from ``examples/fodo/input_fodo_exact.in``.


Analyze

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_fodo_exact.py``

   .. literalinclude:: analysis_fodo_exact.py
      :language: python3
      :caption: You can copy this file from ``examples/fodo/analysis_fodo_exact.py``.
