.. _examples-fourier-coeffs:

Generation of Fourier coefficients from on-axis data
=====================================================

This example illustrates the computation of Fourier coefficients that are used to represent on-axis field or gradient data for soft-edge elements.

Given data in the file ``onaxis_data.in``, execution of the Python script ``python3 fcoef.py`` results in the following output:

``fcoef.out`` - a file containing a list of cosine and sine Fourier coefficients that can be used in ImpactX to define a soft-edge element
``onaxis_data.out`` - a file containing the reconstructed on-axis signal, together with its first and second derivatives

The signal is represented in the form:

g(z) = c0/2 + sum_{j=1}^{nmax}cj*cos(2*pi*j*(z-zmid)/L) + sum_{j=1}^{nmax}sj*sin(2*pi*j*(z-zmid)/L),

where zmid = (zmin+zmax)/2 is the longitudinal location of the midpoint, and L = zmax - zmin is the total length of the z-domain.

The benchmark test uses these coefficients to define and track through a soft-edge quadrupole element.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run as:

* **Python** script: ``python3 run_quadrupole_fcoef.py`` or

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_quadrupole_fcoef.py
          :language: python3
          :caption: You can copy this file from ``examples/fourier_coefficients/run_quadrupole_fcoef.py``.

Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_quadrupole_fcoef.py``

   .. literalinclude:: analysis_quadrupole_fcoef.py
      :language: python3
      :caption: You can copy this file from ``examples/fourier_coefficients/analysis_quadrupole_fcoef.py``.


Visualize
---------

You can run the following script to visualize the reconstruction of the data on-axis:

.. dropdown:: Script ``plot_fcoef.py``

   .. literalinclude:: plot_fcoef.py
      :language: python3
      :caption: You can copy this file from ``examples/fourier_coefficients/plot_fcoef.py``.
