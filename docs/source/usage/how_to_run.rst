.. _usage_run:

Run ImpactX
===========


How to select a tracking mode
-----------------------------

ImpactX can be run using any of three distinct tracking modes.  ImpactX's most powerful tracking mode makes use of symplectic particle tracking with collective effects included (space charge, CSR, wakefields, etc.).
Additionally, ImpactX provides two simplified tracking modes to aid scientists through every step, from beamline inception to operation:
only tracking of the reference particle orbit, or tracking of the beam envelope (6x6 covariance matrix) through linearized transport maps.

================== =============== =============== ==================
Mode               Use Case        Generality      Collective Effects
================== =============== =============== ==================
Particle Tracking  Full Dynamics   Most general    Supported
Reference Tracking Early Design    Reference orbit No
Envelope Tracking  Rapid Scans     Linearized      Not yet
================== =============== =============== ==================


How to select a user interface
------------------------------

ImpactX can be run in three major interfaces:

#. As a Python module (script),
#. As an executable application with a key-value input file, or
#. From a :ref:`graphical user interface (dashboard) <usage-dashboard>`.

While each of these interfaces provide access to the same physics models, you might pick one or the other for specific needs, experience and workflows:

=========== ================= ========================= =========== =========== ============== ===========
Interface   Experience        Good For                  Interactive Extensible  Jupyter        HPC Support
=========== ================= ========================= =========== =========== ============== ===========
Python      Beginner-Advanced Automation, AI/ML         Yes         Yes         Yes            Yes
Application Advanced          Minimal requirements      No          No          Terminal       Yes
Dashboard   Beginner          Learning, Control-Systems Yes         No          Yes            Not yet
=========== ================= ========================= =========== =========== ============== ===========


How to run
----------

After installing ImpactX, run :ref:`one of our examples <usage-examples>`, e.g., :ref:`a FODO cell <examples-fodo>`, like this:

.. tab-set::

   .. tab-item:: Python: Script

      .. code-block:: bash

         python run_fodo.py

   .. tab-item:: Executable: Input File

      .. code-block:: bash

         impactx input_fodo.in

   .. tab-item:: Graphical User Interface

      .. code-block:: bash

         impactx-dashboard


How to select a computing system
--------------------------------

ImpactX supports running on your laptop, alongside edge computers in control systems, in the browser as a dashboard or Jupyter application, on cloud computers, or supercomputers (HPC).
This enables leveraging ImpactX through any stage of a particle accelerator's life cycle.

==================== =================================
Compute Location     Good for
==================== =================================
Local Computer       Initial designs
Edge/Control System  Operations, Digital Twins
Supercomputer        Large Design Studies, ML training
==================== =================================

Additionally, ImpactX runs on CPUs (from AMD, Intel, IBM, ARM, etc.) and modern hardware like GPUs (from Nvidia, AMD or Intel).
As a rough guidance, they are best used for:

======== ========================================================= ======= =================
Hardware Resolution needs                                          AI/ML   Energy Efficiency
======== ========================================================= ======= =================
CPUs     Small (<10M particles) runs and coarse collective effects Yes     Low
GPUs     Many (>10M) particles and/or detailed collective effects  Fastest High
======== ========================================================= ======= =================
