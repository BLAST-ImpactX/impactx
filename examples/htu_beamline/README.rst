.. _examples-htu-beamline:

Dynamics in the HTU Beamline
============================

This example involves tracking a 25 pC electron bunch with 100 MeV total energy through the BELLA Hundred-Terawatt Undulator (HTU) beamline at LBNL `[1] <https://doi.org/10.1103/vh62-gz1p>`__, `[2] <https://doi.org/10.1117/12.3056776>`__.

.. figure:: https://gist.githubusercontent.com/ax3l/cf3dcb517c8d5db293a257a332c87084/raw/5511edec168b6a33d806c11b9dadd01e99c72d52/htu.png
   :alt: Survey plot of the BELLA Hundred-Terawatt Undulator (HTU) beamline at LBNL.

   Survey plot of the BELLA Hundred-Terawatt Undulator (HTU) beamline at LBNL.
   This plot can be generated with ``sim.lattice.plot_survey(ref=ref)`` (see :py:meth:`impactx.elements.KnownElementsList.plot_survey`).

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

      .. dropdown:: HTU Beamline: Lattice File
         :color: light
         :icon: info
         :animate: fade-in-slide-down

         .. literalinclude:: htu_lattice.py
            :language: python3
            :caption: You can copy this file from ``examples/htu_beamline/htu_lattice.py``.

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



.. _htu-beamline-off-energy:

Off-energy transport in the HTU Beamline
==========================================

This is a modification of the example htu-beamline, as follows.  The beam total energy is increased from 100 MeV to 150 MeV, and the chicane :math:`R_{56}` is decreased from 200 microns to 0.

The purpose of this example is to compare two distinct methods for transporting highly off-energy beams:

1) keep the 100 MeV reference energy, but displace the mean energy of the initial beam distribution by setting the parameter mean_pt
2) modify the reference energy to 150 MeV, and rescale the input parameters of the initial beam distribution to yield a distribution with the same second moments as 1)

The bunch is generated in practice from a laser-plasma accelerator stage, resulting in a small intial rms beam size (3.9, 3.9, 1.0) microns, 2 mrad transverse divergence and 2.5% energy spread.  Due
to the large energy spread, chromatic focusing effects are important.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.

The two methods described in 1) and 2) are physically equivalent and the test verifies that the two methods produce physically equivalent output.

Run
---

This example can be run as:

* **Python** script: ``python3 run_impactx_offenergy1.py`` or ``python3 run_impactx_offenergy2.py``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

      .. dropdown:: HTU Beamline: Lattice File
         :color: light
         :icon: info
         :animate: fade-in-slide-down

         .. literalinclude:: htu_lattice.py
            :language: python3
            :caption: You can copy this file from ``examples/htu_beamline/htu_lattice.py``.

      .. literalinclude:: run_impactx_offenergy1.py
         :language: python3
         :caption: You can copy this file from ``examples/htu_beamline/run_impactx_offenergy1.py``.

      .. literalinclude:: run_impactx_offenergy2.py
         :language: python3
         :caption: You can copy this file from ``examples/htu_beamline/run_impactx_offenergy2.py``.

The two scripts are physically equivalent, corresponding to methods 1) and 2) above, and should produce physically equivalent output.

Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_htu_beamline_offenergy.py``

   .. literalinclude:: analysis_htu_beamline_offenergy.py
      :language: python3
      :caption: You can copy this file from ``examples/htu_beamline/analysis_htu_beamline_offenergy.py``.
