.. _examples-scraping:

Expanding Beam Scraping Against a Vacuum Pipe
=============================================

This example describes a coasting bunch, expanding transversely and encountering the aperture defined by the vacuum pipe.
Space charge is neglected, making the problem analytically soluble.

We use a cold (zero emittance) 250 MeV proton bunch whose
initial distribution is a uniformly-populated cylinder of transverse radius :math:`r_b = 2 \mathrm{mm}` with zero momentum spread.

The beam propagates in a drift with a vacuum chamber radius of :math:`R = 3.5 \mathrm{mm}`.

To generate an expanding beam, a linear map is first applied.  This map applies a radial kick to each particle that is proportional to the particle's initial distance from the axis.
This induces a phase space correlation within the beam, such that :math:`p_x = k \cdot x` and :math:`p_y = k \cdot y`, similar to what would be induced by a space charge kick.

The beam remains cylindrical with zero emittance during its evolution in a 6 m drift.
In the absence of an aperture, the beam radius evolves as :math:`r_b(s) = r_b(1 + k\cdot s)`.
In the presence of an aperture, particles are lost during the transverse expansion.  The fraction of charge remaining after a distance s is given by:

.. math::

   \frac{Q_s}{Q_0} = \min\left[1,R^2/(r_b^2(1+s\cdot k)^2)\right].

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.

In addition, the initial and final values of :math:`\sigma_{p_x}`, :math:`\sigma_{p_y}`, and :math:`\sigma_{p_t}` must agree with nominal values.

Finally, the fraction of charge lost against the aperture at the exit of the drift must agree with nominal values.

The physical problem is defined by four relevant parameters, defined within ``run_scraping.py``, that can be modified by the user:

.. code-block:: python

   # problem parameters
   beam_radius = r_b
   aperture_radius = R
   correlation_k = k
   drift_distance = s

These parameters should also be modified inside ``analysis_scraping.py`` for testing.


Run
---

This example can be run as a Python script (``python3 run_scraping.py``) or with an app with an input file (``impactx input_scraping.in``).
Each can also be prefixed with an `MPI executor <https://www.mpi-forum.org>`__, such as ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python Script

       .. literalinclude:: run_scraping.py
          :language: python3
          :caption: You can copy this file from ``examples/scraping_beam/run_scraping.py``.

   .. tab-item:: App Input File

       .. literalinclude:: input_scraping.in
          :language: ini
          :caption: You can copy this file from ``examples/scraping_beam/input_scraping.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_scraping.py``

   .. literalinclude:: analysis_scraping.py
      :language: python3
      :caption: You can copy this file from ``examples/scraping_beam/analysis_scraping.py``.
