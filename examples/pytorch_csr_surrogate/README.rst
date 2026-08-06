.. _examples-ml-csr-surrogate:

Neural Network CSR Surrogate: Chicane
=====================================

This example runs the :ref:`Berlin-Zeuthen magnetic bunch compression chicane <examples-chicane>` with coherent synchrotron radiation (CSR), replacing the built-in analytic CSR wake model with a neural network surrogate that is coupled in through the user-facing Python interface ``sim.csr_kick_model``.
The approach follows:

- Edelen A L, Mayes C E, Emma C and Roussel R.
  **Neural Network Solver for Coherent Synchrotron Radiation Wakefield Calculations in Accelerator-Based Charged Particle Beams**.
  13th International Particle Accelerator Conference (IPAC'22), WEPOMS013, 2022.
  `DOI:10.18429/JACoW-IPAC2022-WEPOMS013 <https://doi.org/10.18429/JACoW-IPAC2022-WEPOMS013>`__
  (`arXiv:2203.07542 <https://arxiv.org/abs/2203.07542>`__)

Per tracking slice inside each bend, ImpactX deposits the binned longitudinal charge profile :math:`\lambda(t)` and hands it, together with the element context (bend radius, position in the element, reference particle), to the user-defined kick model.
The model returns the per-bin longitudinal CSR kick force, which ImpactX applies to the beam.

In this self-contained example, the network is trained on the fly (a few seconds) against the steady-state CSR wake model: it learns the map from the *normalized* charge profile shape to the dimensionless per-bin kick shape :math:`G`.
The physical kick then follows from the exact scaling of the steady-state model,

.. math::

   K = \frac{Q \, \kappa(R)}{q_e \, \Delta^{4/3}} \, G(\text{shape}),
   \qquad
   \kappa(R) = \frac{q_e^2}{2 \pi \varepsilon_0 \, 3^{1/3} R^{2/3}},

with bunch charge :math:`Q`, bin size :math:`\Delta` and bend radius :math:`R`, so a single trained network generalizes over bunch charge, compression stage (bunch length), and bend radius.
In production use, a pre-trained model (e.g., trained on data from a high-fidelity solver, including transient and 2D/3D CSR effects) would be loaded from disk instead.

This example does not require ImpactX to be compiled with FFT support: the built-in analytic CSR model is fully replaced by the surrogate.

In this test, the initial values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values, and the final values must agree with the built-in CSR model results within the surrogate approximation error, including the CSR-induced horizontal emittance growth.

Run
---

This example can **only** be run with **Python**:

* **Python** script: ``python3 run_chicane_csr_ml.py``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. literalinclude:: run_chicane_csr_ml.py
   :language: python
   :caption: You can copy this file from ``examples/pytorch_csr_surrogate/run_chicane_csr_ml.py``.

Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_chicane_csr_ml.py``

   .. literalinclude:: analysis_chicane_csr_ml.py
      :language: python
      :caption: You can copy this file from ``examples/pytorch_csr_surrogate/analysis_chicane_csr_ml.py``.
