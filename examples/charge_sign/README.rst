.. _examples-charge-sign:

Testing Charge and Field Sign Consistency
==========================================

ImpactX uses a normalized field (focusing) strength internally for tracking. For example, in a Quad a focusing strength k > 0 indicates horizontal focusing, independent of the charge
of the particle (just as in MAD-X).

For elements that support an option for user-provided field strength in SI units (e.g. ChrQuad and ExactSbend), varying the sign of the charge should provide tracking consistent
with the sign of the specified field.

This example tests the transport of a 1 GeV electron bunch through multiple field sign combinations.  The transport occurs through pairs of elements.  For example, the element
``quad1`` specifies a positive normalized field strength, and should provide horizontal focusing.  The element ``quad2inv`` specifies a quadrupole field gradient in units of T/m, which
should provide an equivalent horizontal focusing.  (For electrons, the required field gradient is negative.)

To test this, the map for ``quad1`` is applied (ds > 0), and this is followed by applying the inverse map for ``quad2inv`` (ds < 0), which should result in the identity map.  Similar
considerations apply to the other element pairs.

As a result, the second moments of x, y, and t and the associated emittances of the bunch should all be exactly unchanged.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_charge_sign.py`` or
* ImpactX **executable** using an input file: ``impactx input_charge_sign.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_charge_sign.py
          :language: python3
          :caption: You can copy this file from ``examples/charge_sign/run_charge_sign.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_charge_sign.in
          :language: ini
          :caption: You can copy this file from ``examples/charge_sign/input_charge_sign.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_charge_sign.py``

   .. literalinclude:: analysis_charge_sign.py
      :language: python3
      :caption: You can copy this file from ``examples/charge_sign/analysis_charge_sign.py``.
