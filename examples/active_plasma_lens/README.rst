.. _examples-active_plasma_lens:

Active Plasma Lens
==================

These examples demonstrate the effect of an Active Plasma Lens (APL) on the beam.
The lattice contains this element and nothing else.
The length of the element is 20 mm, and it can be run in no-field, focusing, and defocusing mode.

We use a 200 MeV electron beam with an initial normalized rms emittance of 10 um.
The beam is set to have :math:`\alpha = 0` in the middle of the lens in the case of no field.
The beam size in the middle of the lens is set to 10 µm for the no-field examples (in order to have a strongly parabolic :math:`\beta`-function within the lens), and 100 µm for the focusing and defocusing examples.
A :math:`\sigma_{pt} = 10^{-3}` is also assumed.
Before the simulation, this beam is back-propagated to the lens entry assuming zero field.

Run
---

This example can be run as
* ``python3 run_APL_ChrPlasmaLens_zero.py`` (no field, ``ChrPlasmaLens``, tracking)
* ``python3 run_APL_ChrPlasmaLens_focusing.py`` (focusing field, ``ChrPlasmaLens``, tracking)
* ``python3 run_APL_ChrPlasmaLens_defocusing.py`` (defocusing field, ``ChrPlasmaLens``, tracking)

These all use the library ``run_APL.py`` internally to create the simulations.

Analyze
-------

We run the following scripts to analyze correctness of the output:
* ``python3 analysis_APL_ChrPlasmaLens_zero.py`` (no field, ``ChrPlasmaLens``, tracking)
* ``python3 analysis_APL_ChrPlasmaLens_focusing.py`` (focusing field, ``ChrPlasmaLens``, tracking)
* ``python3 analysis_APL_ChrPlasmaLens_defocusing.py`` (defocusing field, ``ChrPlasmaLens``, tracking)

These all use the library ``analysis_APL_ChrPlasmaLens.py`` internally.

Visualize
---------
You can run the following scripts to visualize the beam evolution over time (e.g. :math:`s`):
* ``python3 s_APL_ChrPlasmaLens_zero.py`` (no field, ``ChrPlasmaLens``, tracking)plot
* ``python3 plot_APL_ChrPlasmaLens_focusing.py`` (focusing field, ``ChrPlasmaLens``, tracking)
* ``python3 plot_APL_ChrPlasmaLens_defocusing.py`` (defocusing field, ``ChrPlasmaLens``, tracking)

These all use the library ``plot_APL_ChrPlasmaLens.py`` internally.

Additionally, it is also possible to run ``python3 plot_APL_ChrPlasmaLens_analytical.py``, which plots the expected Twiss :math:`\alpha` and :math:`\beta` functions at the end of the lens as a function of the lens gradient. This uses the stand-alone Twiss propagation function ``analytic_final_estimate()`` from ``run_APL.py``.
