.. _examples-fodo:

Polygon Aperture
=========

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

       .. literalinclude:: input_fodo.in
          :language: ini
          :caption: You can copy this file from ``examples/polygon_aperture/input_fodo.in``.

Other examples are
..
   * Aperture that absorbs
      :Python:
         ``examples/polygon_aperture/run_polygon_aperture_absorb.py``
      :Input file:
         ``examples/polygon_aperture/input_polygon_aperture_absorb.in``
   * Aperture with rotation that absorbs
      :Python:
         ``examples/polygon_aperture/run_polygon_aperture_absorb_rotate.py``
      :Input file:
         ``examples/polygon_aperture/input_polygon_aperture_absorb_rotate.in``
   * Aperture with offset that absorbs
      :Python:
         ``examples/polygon_aperture/run_polygon_aperture_absorb_offset.py``
      :Input file:
         ``examples/polygon_aperture/input_polygon_aperture_absorb_offset.in``
   * Aperture with offset and rotation that absorbs
      :Python:
         ``examples/polygon_aperture/run_polygon_aperture_absorb_rotate_offset.py
      :Input file:
         ``examples/polygon_aperture/input_polygon_aperture_absorb_rotate_offset.in``


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

.. figure:: https://private-user-images.githubusercontent.com/16342668/513572482-f9029000-6327-45b2-bb10-5b5edbe60373.png?jwt=eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3NjI5ODc2NTcsIm5iZiI6MTc2Mjk4NzM1NywicGF0aCI6Ii8xNjM0MjY2OC81MTM1NzI0ODItZjkwMjkwMDAtNjMyNy00NWIyLWJiMTAtNWI1ZWRiZTYwMzczLnBuZz9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFWQ09EWUxTQTUzUFFLNFpBJTJGMjAyNTExMTIlMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjUxMTEyVDIyNDIzN1omWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPTViYzRkMDU2ZTZjMDkyNWM0M2Q2YmU3NGVjMWJhN2RiNWRlZDc2ZDdlNmI0ODg1NTY1MzdjMmNlYjRmOTg1MWImWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0In0.Bmzd6hFUj6PK172blKBA0zuBSn5ItMtjy3-pRZP3fhQ
   :alt: Initial and transmitted particles through the example polygon aperture.
   