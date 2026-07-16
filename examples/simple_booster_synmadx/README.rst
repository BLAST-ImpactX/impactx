.. _examples-simple-booster-synmadx:

Simple Booster (synmadx MAD-X Parser)
=====================================

This is the same simplified model of the Fermilab Booster as the
:ref:`Simple Booster example <examples-simple-booster>`,
with 24 cells and 22 (inactive) RF cavities,
run at the injection energy.

Instead of a hand-converted ImpactX lattice file, this example parses the original
MAD-X file ``sbbooster-cooked.madx`` at runtime with *synmadx*, a standalone MAD-X
lattice parser that was ported from
`Synergia <https://github.com/fnalacceleratormodeling/synergia2>`__.
synmadx is an alternative MAD-X parser, additionally to our official MAD-X importer,
:py:meth:`impactx.elements.KnownElementsList.load_file`.
The parsed Synergia lattice is then converted element-by-element into ImpactX
elements by the ``impactx.synmadx.syn2_to_impactx`` function, which also sets up
the RF cavities programmatically.

.. note::

   This example requires ImpactX to be compiled with the CMake option
   ``-DImpactX_SYNMADX=ON``, which needs `Boost <https://www.boost.org>`__ installed.

The matched Twiss parameters at entry, calculated with Synergia from the
``sbbooster-cooked.madx`` file, are:

* :math:`\beta_\mathrm{x} = 33.73645362843065243` m
* :math:`\alpha_\mathrm{x} = -0.01298673960026007664`
* :math:`\beta_\mathrm{y} = 5.252517912567207681` m
* :math:`\alpha_\mathrm{y} = 0.006089861210659328755`
* :math:`\mathrm{D}_\mathrm{x} = 3.187407765856291153` m
* :math:`\mathrm{Dp}_\mathrm{x} = 0.001136005067625678322`

The initial beam parameters follow the PIP-II Booster beam at injection
as described in the Conceptual Design Report.
The beam consists of protons at a kinetic energy of 800 MeV with emittances
specified as:

+------------------------+--------------------------------------+
| :math:`\epsilon_{x}`   | :math:`16 \pi` mm-mr normalized 95%  |
+------------------------+--------------------------------------+
| :math:`\epsilon_{y}`   | :math:`16 \pi` mm-mr normalized 95%  |
+------------------------+--------------------------------------+
| :math:`\epsilon_{L}`   | 0.1 eV-s 97%                         |
+------------------------+--------------------------------------+


Run
---

This example can only be run with a Python script:

* **Python** script: ``python3 run_parse_booster_lattice.py``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_parse_booster_lattice.py
          :language: python3
          :caption: You can copy this file from ``examples/simple_booster_synmadx/run_parse_booster_lattice.py``. The file ``sbbooster-cooked.madx`` from the same directory is also required.

   .. tab-item:: Lattice Converter

      .. literalinclude:: ../../src/python/impactx/synmadx/syn2_to_impactx.py
         :language: python3
         :caption: Implementation of the conversion functions that ship with ImpactX, directly available to users as ``impactx.synmadx.syn2_to_impactx`` and ``impactx.synmadx.unroll_impactx_lattice``.

   .. tab-item:: MAD-X: Script

       .. literalinclude:: sbbooster-cooked.madx
          :language: text
          :caption: Original MAD-X lattice file that describes the simple Booster model, available at ``examples/simple_booster_synmadx/sbbooster-cooked.madx``.


Visualize
---------

You can run the following script to visualize the beam evolution:

.. dropdown:: Script ``analyze_diags.py``

   .. literalinclude:: analyze_diags.py
      :language: python3
      :caption: You can copy this file from ``examples/simple_booster_synmadx/analyze_diags.py``.
