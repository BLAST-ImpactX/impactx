.. _examples-rfcavity:

Acceleration by RF Cavities
===========================

Beam accelerated through a sequence of 4 RF cavities (without space charge).

We use a 230 MeV electron beam with initial normalized rms emittance of 1 um.

The lattice and beam parameters are based on Example 2 of the IMPACT-Z examples folder:

https://github.com/impact-lbl/IMPACT-Z/tree/master/examples/Example2

The final target beam energy and beam moments are based on simulation in
IMPACT-Z, without space charge.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_rfcavity.py`` or
* ImpactX **executable** using an input file: ``impactx input_rfcavity.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_rfcavity.py
          :language: python3
          :caption: You can copy this file from ``examples/rfcavity/run_rfcavity.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_rfcavity.in
          :language: ini
          :caption: You can copy this file from ``examples/rfcavity/input_rfcavity.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_rfcavity.py``

   .. literalinclude:: analysis_rfcavity.py
      :language: python3
      :caption: You can copy this file from ``examples/rfcavity/analysis_rfcavity.py``.


.. _examples-rfcavity-ref-part:

Acceleration by RF Cavities (Reference Particle Tracking)
=========================================================

This test is identical to the test examples-rfcavity above, except that the code is run in reference orbit tracking mode.

It is used to validate the numerical integration of the reference energy gain in a chain of 4 RF cavities, and to demonstrate this tracking mode.

In this test, the initial and final reference values of :math:`s` and :math:`\gamma` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_rfcavity_ref_part.py`` or
* ImpactX **executable** using an input file: ``impactx input_rfcavity_ref_part.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_rfcavity_ref_part.py
          :language: python3
          :caption: You can copy this file from ``examples/rfcavity/run_rfcavity_ref_part.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_rfcavity_ref_part.in
          :language: ini
          :caption: You can copy this file from ``examples/rfcavity/input_rfcavity_ref_part.in``.

Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_rfcavity_ref_part.py``

   .. literalinclude:: analysis_rfcavity_ref_part.py
      :language: python3
      :caption: You can copy this file from ``examples/rfcavity/analysis_rfcavity_ref_part.py``.



.. _examples-rfcavity-ref-part-hook:

Acceleration by RF Cavities (Using a Callback Hook)
===================================================

This test is similar to the :ref:`test above <examples-rfcavity-ref-part>`, except that it illustrates the use of the callback hook feature :py:attr:`~impactx.ImpactX.hook` to provide an alternative method to set the RF cavity field amplitude(s) and phase(s).

The two functions ``get_synch_phase_Veff`` and ``get_phase_emax`` allow the user to convert between the pair of inputs (``Veff, phase_synch``) and (``escale,phase``).  Here ``escale`` and ``phase`` are the documented ImpactX inputs in :py:class:`~impactx.elements.RFCavity`, while ``Veff`` (:math:`V_{\rm eff}`) and ``phase_synch`` (:math:`\phi_s`) denote the effective voltage and synchronous phase of the cavity, defined here such that:

.. math::

    \Delta\gamma = V_{\rm eff}\cos(\phi_s).

The conversion is implemented under the approximation that the relativistic beta does not vary within the cavity.  (E.g., this is an excellent approximation for large values of gamma.)

For the definition of effective voltage and synchronous phase, see:

* T. Wangler, Principles of RF Linear Accelerators, John Wiley and Sons, Inc., New York, 1998, Sections 2.2 - 2.5
* D. Uriot and N. Pichoff, TraceWin Manual, CEA/SACLAY - DRF/Irfu/DACM - Didiear URIOT, updated May 31, 2023, "Synchronous phase definitions"

The equivalent voltage here is multiplied by the particle charge :math:`q` and divided by :math:`mc^2`, so this is a dimensionless quantity.

For each of the four RF cavities, the equivalent voltage and synchronous phase are computed.  Those values are then converted back to the default phase and escale inputs, whose values are updated within the cavity.

When implemented correctly, the dynamics of the reference particle within these cavities should be identical to the dynamics seen in :ref:`single reference particle tracking <examples-rfcavity-ref-part>`.

In this test, the initial and final reference values of :math:`s` and :math:`\gamma` must agree with nominal values.


Run
---

This example can be run as:

* **Python** script: ``python3 run_rfcavity_ref_part_hook.py`

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_rfcavity_ref_part_hook.py
          :language: python3
          :caption: You can copy this file from ``examples/rfcavity/run_rfcavity_ref_part_hook.py``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_rfcavity_ref_part_hook.py``

   .. literalinclude:: analysis_rfcavity_ref_part_hook.py
      :language: python3
      :caption: You can copy this file from ``examples/rfcavity/analysis_rfcavity_ref_part_hook.py``.
