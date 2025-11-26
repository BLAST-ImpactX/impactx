.. _examples-active_plasma_lens:

Active Plasma Lens
==================

These examples demonstrate the effect of an Active Plasma Lens (APL) on the beam.
The lattice contains this element and nothing else.
The length of the element is 20 mm, and it can be run in no-field, focusing, and defocusing mode.
It is implemented with two different elements, ```ChrPlasmaLens``` and ```ConstF```;
for ```ConstF``` we test both particle tracking and envelope tracking.

We use a 200 MeV electron beam with an initial normalized rms emittance of 10 µm.
The initial Twiss parameters for the simulations are set by first assuming the beam
to have :math:`\alpha = 0` in the middle of the lens, and then the we backwards-propagate
this analytically to the start of the lens, under the assumption of no field.
The beam is then forward-propagated using ImpactX and the different cases for
element strength and simulation method.

The beam size in the middle of the lens is set to 10 µm for the no-field examples
in order to have a strongly parabolic :math:`\beta`-function within the lens,
and 100 µm for the focusing and defocusing examples.
The :math:`\beta` and :math:`\gamma`-function in the middle of the lens is calculated
from this beam size and the assumed emittance.
A :math:`\sigma_{pt} = 10^{-3}` is also assumed for the tracking examples.

Run
---

This example can be run as Python scripts, with the names
indicating the element used (``ChrPlasmaLens`` or ``ConstF``),
the field in the lens (``zero``, ``focusing``, or ``defocusing``),
and simulation type (``tracking`` or ``envelope``).

These both run the simulation, and produce analytical reference parameters
which are used for comparison by the analysis scripts.

* ``python3 run_APL_ChrPlasmaLens_tracking_zero.py``
* ``python3 run_APL_ChrPlasmaLens_tracking_focusing.py``
* ``python3 run_APL_ChrPlasmaLens_tracking_defocusing.py``
* ``python3 run_APL_ConstF_tracking_zero.py``
* ``python3 run_APL_ConstF_tracking_focusing.py``
* ``python3 run_APL_ConstF_tracking_defocusing.py``
* ``python3 run_APL_ConstF_envelope_zero.py``
* ``python3 run_APL_ConstF_envelope_focusing.py``
* ``python3 run_APL_ConstF_envelope_defocusing.py``

These all use the library ``run_APL.py`` internally to setup and and run the simulations.

.. dropdown:: ``run_APL.py``

   .. literalinclude:: run_APL.py
      :language: python3
      :caption: You can copy this file from ``examples/active_plasma_lens/run_APL.py``.

The scripts used to start the simulations:

.. dropdown:: ``run_APL_ChrPlasmaLens_tracking_zero.py``

   .. literalinclude:: run_APL_ChrPlasmaLens_tracking_zero.py
      :language: python3
      :caption: You can copy this file from ``examples/active_plasma_lens/run_APL_ChrPlasmaLens_tracking_zero.py``.

.. dropdown:: ``run_APL_ChrPlasmaLens_tracking_focusing.py``

   .. literalinclude:: run_APL_ChrPlasmaLens_tracking_focusing.py
      :language: python3
      :caption: You can copy this file from ``examples/active_plasma_lens/run_APL_ChrPlasmaLens_tracking_focusing.py``.

.. dropdown:: ``run_APL_ChrPlasmaLens_tracking_defocusing.py``

   .. literalinclude:: run_APL_ChrPlasmaLens_tracking_defocusing.py
      :language: python3
      :caption: You can copy this file from ``examples/active_plasma_lens/run_APL_ChrPlasmaLens_tracking_defocusing.py``.

.. dropdown:: ``run_APL_ConstF_tracking_zero.py``

   .. literalinclude:: run_APL_ConstF_tracking_zero.py
      :language: python3
      :caption: You can copy this file from ``examples/active_plasma_lens/run_APL_ConstF_tracking_zero.py``.

.. dropdown:: ``run_APL_ConstF_tracking_focusing.py``

   .. literalinclude:: run_APL_ConstF_tracking_focusing.py
      :language: python3
      :caption: You can copy this file from ``examples/active_plasma_lens/run_APL_ConstF_tracking_focusing.py``.

.. dropdown:: ``run_APL_ConstF_tracking_defocusing.py``

   .. literalinclude:: run_APL_ConstF_tracking_defocusing.py
      :language: python3
      :caption: You can copy this file from ``examples/active_plasma_lens/run_APL_ConstF_tracking_defocusing.py``.

.. dropdown:: ``run_APL_ConstF_envelope_zero.py``

   .. literalinclude:: run_APL_ConstF_envelope_zero.py
      :language: python3
      :caption: You can copy this file from ``examples/active_plasma_lens/run_APL_ConstF_envelope_zero.py``.

.. dropdown:: ``run_APL_ConstF_envelope_focusing.py``

   .. literalinclude:: run_APL_ConstF_envelope_focusing.py
      :language: python3
      :caption: You can copy this file from ``examples/active_plasma_lens/run_APL_ConstF_envelope_focusing.py``.

.. dropdown:: ``run_APL_ConstF_envelope_defocusing.py``

   .. literalinclude:: run_APL_ConstF_envelope_defocusing.py
      :language: python3
      :caption: You can copy this file from ``examples/active_plasma_lens/run_APL_ConstF_envelope_defocusing.py``.

Analyze
-------

The following scripts can be used to analyze correctness of the output,
by comparing it to a reference output that is produced and outputed to
the standard output (terminal) from the run scripts.

The output should be the same accross elements (``ConstF`` or ``ChrPlasmaLens``),
but depend on the field in the lens (``zero``, ``focusing``, or ``defocusing``),
and simulation type (``tracking`` or ``envelope``).
The analysis scripts are therefore the same for both element types.

All analysis scripts look at the output from most recent simulation run in
the current working directory, i.e. the ``diags`` folder.

* ``python3 analysis_APL_tracking_zero.py``
* ``python3 analysis_APL_tracking_focusing.py``
* ``python3 analysis_APL_tracking_defocusing.py``
* ``python3 analysis_APL_envelope_zero.py``
* ``python3 analysis_APL_envelope_focusing.py``
* ``python3 analysis_APL_envelope_defocusing.py``

These all use the library ``analysis_APL.py`` internally to load data etc:

.. dropdown:: ``analysis_APL.py``

   .. literalinclude:: analysis_APL.py
      :language: python3
      :caption: You can copy this file from ``examples/active_plasma_lens/analysis_APL.py``.

The analysis scripts are:

.. dropdown:: ``analysis_APL_tracking_zero.py``

   .. literalinclude:: analysis_APL_tracking_zero.py
      :language: python3
      :caption: You can copy this file from ``examples/active_plasma_lens/analysis_APL_tracking_zero.py``.

.. dropdown:: ``analysis_APL_tracking_focusing.py``

   .. literalinclude:: analysis_APL_tracking_focusing.py
      :language: python3
      :caption: You can copy this file from ``examples/active_plasma_lens/analysis_APL_tracking_focusing.py``.

.. dropdown:: ``analysis_APL_tracking_defocusing.py``

   .. literalinclude:: analysis_APL_tracking_defocusing.py
      :language: python3
      :caption: You can copy this file from ``examples/active_plasma_lens/analysis_APL_tracking_defocusing.py``.

.. dropdown:: ``analysis_APL_envelope_zero.py``

   .. literalinclude:: analysis_APL_envelope_zero.py
      :language: python3
      :caption: You can copy this file from ``examples/active_plasma_lens/analysis_APL_envelope_zero.py``.

.. dropdown:: ``analysis_APL_envelope_focusing.py``

   .. literalinclude:: analysis_APL_envelope_focusing.py
      :language: python3
      :caption: You can copy this file from ``examples/active_plasma_lens/analysis_APL_envelope_focusing.py``.

.. dropdown:: ``analysis_APL_envelope_defocusing.py``

   .. literalinclude:: analysis_APL_envelope_defocusing.py
      :language: python3
      :caption: You can copy this file from ``examples/active_plasma_lens/analysis_APL_envelope_defocusing.py``.


Visualize
---------
You can run the following scripts to visualize the beam evolution over time (e.g. :math:`s`),
and compare to analytical expectation.
Here, For this, the output format is identical accross the element- and simulation-types,
only depending on the selected lens field (``zero``, ``focusing``, or ``defocusing``).

* ``python3 plot_APL_zero.py``
* ``python3 plot_APL_focusing.py``
* ``python3 plot_APL_defocusing.py``

These all use the library ``plot_APL.py`` internally,
as well as ``run_APL.py`` and ``analysis_APL.py`` which are described above.

.. dropdown:: ``plot_APL.py``

   .. literalinclude:: plot_APL.py
      :language: python3
      :caption: You can copy this file from ``examples/active_plasma_lens/plot_APL.py``.

The plotting scripts are given below, together with their graphical output.
For the plots, output from the ``ChrPlasmaLens_tracking`` simulations were used,
which shows some error in the envelope due to statistical fluctuations in the initial particle distribution.

.. dropdown:: ``plot_APL_zero.py`` and output figures:

   .. literalinclude:: plot_APL_zero.py
      :language: python3
      :caption: You can copy this file from ``examples/active_plasma_lens/plot_APL_zero.py``.

   .. figure:: APL_zero-sigma_REF.png
      :alt: The calculated :math:`\sigma` values throughout the lens without field, for both planes and a separate analytical estimate.

.. dropdown:: ``plot_APL_focusing.py`` and output figures:

   .. literalinclude:: plot_APL_focusing.py
      :language: python3
      :caption: You can copy this file from ``examples/active_plasma_lens/plot_APL_focusing.py``.

   .. figure:: APL_focusing-sigma_REF.png
      :alt: The calculated :math:`\sigma` values throughout the lens with focusing field, for both planes and a separate analytical estimate.

.. dropdown:: ``plot_APL_defocusing.py`` and output figures:

   .. literalinclude:: plot_APL_defocusing.py
      :language: python3
      :caption: You can copy this file from ``examples/active_plasma_lens/plot_APL_defocusing.py``.

   .. figure:: APL_defocusing-sigma_REF.png
      :alt: The calculated :math:`\sigma` values throughout the lens with defocusing field, for both planes and a separate analytical estimate.

Additionally, it is also possible to run ``python3 plot_APL_analytical.py``,
which plots the expected Twiss :math:`\alpha` and :math:`\beta` functions at the end of the lens
as a function of the lens gradient. This uses the stand-alone Twiss propagation function
``analytic_final_estimate()`` from ``run_APL.py``.

.. dropdown:: ``plot_APL_analytical.py`` and output figures:

   .. literalinclude:: plot_APL_analytical.py
      :language: python3
      :caption: You can copy this file from ``examples/active_plasma_lens/plot_APL_analytical.py``.

   .. figure:: APL_analytical_sqrtBeta_REF.png
      :alt: The analytically computed :math:`\sqrt{\beta}` value at the end of the lens (proportional to beam size), as a function of gradient.

   .. figure:: APL_analytical_alpha_REF.png
      :alt: The analytically computed :math:`\alpha` value at the end of the lens, as a function of gradient.
