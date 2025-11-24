.. _examples-polygon-aperture:

Polygon Aperture
================

A 2D transverse aperture of a closed polygon defined by the :math:`x` and
:math:`y` coordinates of the vertices.
The vertices must be ordered in the counter-clockwise direction and must close,
i.e. the
first and last coordinates must be the same.

This example takes a 800 MeV proton beam generated as a waterbag
distribution with :math:`\sigma_x`, :math:`\sigma_y` both equal to
2 mm impinging directly the mask.

Several variations are given, with the mask either transmitting or
blocking the particles, also with option rotation and transverse offset.

Run
---

This example of a transmitting mask can be run **either** as:

* **Python** script: ``python3 run_polygon_aperture.py`` or
* ImpactX **executable** using an input file: ``impactx input_polygon_aperture.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

      .. literalinclude:: run_polygon_aperture.py
         :language: python3
         :caption: You can copy this file from ``examples/polygon_aperture/run_polygon_aperture.py``.

   .. tab-item:: Executable: Input File

      .. literalinclude:: input_polygon_aperture.in
         :language: ini
         :caption: You can copy this file from ``examples/polygon_aperture/input_polygon_aperture.in``.

Other examples are

.. tab-set::

   .. tab-item:: Aperture that absorbs

      .. literalinclude:: run_polygon_aperture_absorb.py
         :language: python3
         :caption: You can copy this file from ``examples/polygon_aperture/run_polygon_aperture_absorb.py``.

   .. tab-item:: Aperture with rotation that absorbs

      .. literalinclude:: run_polygon_aperture_absorb_rotate.py
         :language: python3
         :caption: You can copy this file from ``examples/polygon_aperture/run_polygon_aperture_absorb_rotate.py``.

   .. tab-item:: Aperture with offset that absorbs

      .. literalinclude:: run_polygon_aperture_absorb_offset.py
         :language: python3
         :caption: You can copy this file from run_polygon_aperture_absorb_offset.py

   .. tab-item:: Aperture with offset and rotation that absorbs

      .. literalinclude:: run_polygon_aperture_absorb_rotate_offset.py
         :language: python3
         :caption: You can copy this file from ``examples/polygon_aperture/run_polygon_aperture_absorb_rotate_offset.py``.

Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_polygon_aperture.py``

   .. literalinclude:: analysis_polygon_aperture.py
      :language: python3
      :caption: You can copy this file from ``examples/polygon_aperture/analysis_polygon_aperture.py``.

The number of surviving particles is printed and checked.

Visualize
---------

You can run the following script to visualize aperture effect:

.. dropdown:: Script ``plot_polygon_aperture.py``

   .. literalinclude:: plot_polygon_aperture.py
      :language: python3
      :caption: You can copy this file from ``examples/polygon_aperture/plot_polygon_aperture.py``.

.. figure:: polygon_aperture.png
   :alt: Initial and transmitted particles through the example polygon aperture.
