.. _examples-solenoid-softedge:

Soft-edge solenoid
==================

Proton beam propagating through a 6 m region containing a soft-edge
solenoid.

The solenoid model used is the default thin-shell model described in:
P. Granum et al, "Efficient calculations of magnetic fields of solenoids for simulations,"
NIMA 1034, 166706 (2022)
`DOI:10.1016/j.nima.2022.166706 <https://doi.org/10.1016/j.nima.2022.166706>`__

The solenoid is a cylindrical current sheet with a length of 1 m and a
radius of 0.1667 m, corresponding to an aspect ratio diameter/length = 1/3.
The peak magnetic field on-axis is 3 T.

We use a 250 MeV proton beam with initial unnormalized rms emittance of 1 micron
in the horizontal plane, and 2 micron in the vertical plane.

In this test, the initial and final values of :math:`\sigma_x`, :math:`\sigma_y`, :math:`\sigma_t`, :math:`\epsilon_x`, :math:`\epsilon_y`, and :math:`\epsilon_t` must agree with nominal values.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 run_solenoid_softedge.py`` or
* ImpactX **executable** using an input file: ``impactx input_solenoid_softedge.in``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_solenoid_softedge.py
          :language: python3
          :caption: You can copy this file from ``examples/solenoid_softedge/run_solenoid_softedge.py``.

   .. tab-item:: Executable: Input File

       .. literalinclude:: input_solenoid_softedge.in
          :language: ini
          :caption: You can copy this file from ``examples/solenoid_softedge/input_solenoid_softedge.in``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_solenoid_softedge.py``

   .. literalinclude:: analysis_solenoid_softedge.py
      :language: python3
      :caption: You can copy this file from ``examples/solenoid_softedge/analysis_solenoid_softedge.py``.


.. _solenoid-softedge-solvable:

Exactly-solvable (non-uniform) soft-edge solenoid
=================================================

This benchmark checks the calculation of the linear map for a soft-edge solenoid with a non-uniform longitudinal
on-axis magnetic field profile, in the special case of a field profile for which the map is exactly-solvable.

The test involves 250 MeV protons propagating through a 2 m region with a solenoidal magnetic field.

The user can run this test for various values of the magnetic field strength by modifying the parameter "bscale" in both the input file and in the analysis script.  (The values in the two files must be consistent.)

In this test, all 36 elements of the 6x6 transport matrix must agree with predicted values to within numerical tolerance (currently 1e-7).

.. dropdown:: Details
   :color: light
   :icon: info
   :animate: fade-in-slide-down

    **Introduction**

    Consider a Hamiltonian :math:`H` quadratic in the phase space variables :math:`\zeta=(x,p_x,y,p_y,t,p_t)`.
    Then :math:`H` can be written in the form:

    .. math::

        H(\zeta,z)=\frac{1}{2}\zeta^TS(z)\zeta,

    where :math:`S` is a :math:`6\times 6` symmetric matrix.
    The linear symplectic map :math:`M` describing transport from :math:`z_0\mapsto z` is the unique solution of:

    .. math::

        \frac{dM(z)}{dz}=JS(z)M(z),\quad\quad M(z_0)=I,

    where :math:`J` is the :math:`6\times 6` matrix of the symplectic form.

    From :ref:`R. Ryne's USPAS notes <solenoid-softedge-solvable-ryne>`, the Hamiltonian for linear transport in a solenoid is given by:

    .. math::

        H=H_2^{foc}+H_2^{rot},

    where

    .. math::

        H_2^{foc}&=\frac{p_t^2}{2\beta_0^2\gamma_0^2}+\frac{1}{2}(p_x^2+p_y^2)+\frac{\alpha^2}{2}(x^2+y^2), \\
        H_2^{rot}&=-\alpha(xp_y-yp_x),

    and :math:`\alpha` is
    related to the on-axis (solenoidal) magnetic field :math:`B_z` by:

    .. math::

        \alpha(z)=\frac{1}{2}\frac{B_z(z)}{B\rho}.

    We consider an on-axis solenoid field profile of the following special form:

    .. math::

        B_z(z)=\frac{B_0}{1+(z/g)^2},

    where :math:`B_0`
    is just the peak magnetic field on-axis, and :math:`g` is a length parameter describing the rate of decay of the field away from its peak at :math:`z=0`.

    Since :math:`H_2^{foc}` and
    :math:`H_2^{rot}` Poisson-commute, the linear map :math:`M` can be written in the factorized form:

    .. math::

        M=M^0R,

    where :math:`M^0`
    is the linear map associated with
    :math:`H_2^{foc}`
    , and :math:`R` is the linear map associated with :math:`H_2^{rot}` :ref:`(see ref) <solenoid-softedge-solvable-ryne>`.

    To express the linear map between points :math:`z_0` and :math:`z`, we define the following dimensionless parameters:

    .. math::

        \Delta\phi=\tan^{-1}\left(\frac{z}{g}\right)-\tan^{-1}\left(\frac{z_0}{g}\right),\quad\quad b=\frac{1}{2}\frac{B_0g}{B\rho},\quad\quad \beta=\sqrt{1+b^2}.

    Then we obtain for the focusing matrix:

    .. math::

        M^0_{11}&=M^0_{33}=\sqrt{\frac{g^2+z^2}{g^2+z_0^2}}\left\{\cos(\beta\Delta\phi)-\frac{z_0\sin(\beta\Delta\phi)}{g\beta}\right\}, \\
        M^0_{12}&=M^0_{34}=\sqrt{(g^2+z^2)(g^2+z_0^2)}\frac{\sin(\beta\Delta\phi)}{g\beta}, \\
        M^0_{21}&=M^0_{43}=\frac{\beta(z-z_0)\cos(\beta\Delta\phi)-g(\beta^2+zz_0/g^2)\sin(\beta\Delta\phi)}{\beta\sqrt{(g^2+z^2)(g^2+z_0^2)}}, \\
        M^0_{22}&=M^0_{44}=\sqrt{\frac{g^2+z_0^2}{g^2+z^2}}\left\{\cos(\beta\Delta\phi)+\frac{z\sin(\beta\Delta\phi)}{g\beta}\right\}, \\
        M^0_{56}&=\frac{z-z_0}{\beta_0^2\gamma_0^2},\quad\quad M^0_{55}=M^0_{66}=1.

    Likewise, for the rotation matrix we obtain:

    .. math::

        R_{11}&=R_{22}=\cos(b\Delta\phi),\quad\quad R_{13}=R_{24}=\sin(b\Delta\phi), \\
        R_{31}&=R_{42}=-\sin(b\Delta\phi),\quad\quad R_{33}=R_{44}=\cos(b\Delta\phi), \\
        R_{55}&=R_{66}=1.

    The reader can verify that the map equations and the symplectic condition are satisfied:

    .. math::

        &\frac{dM^0(z)}{dz}=JS^{foc}(z)M^0(z),\quad\quad M^0(z_0)=I, \quad\quad (M^0)^TJM^0=J\\
        &\frac{dR(z)}{dz}=JS^{rot}(z)R(z),\quad\quad R(z_0)=I, \quad\quad R^TJR=J.

    Here the :math:`6\times 6`
    symmetric matrices :math:`S^{foc}`
    and :math:`S^{rot}` are defined such that:

    .. math::

        H^{foc}(\zeta,z)=\frac{1}{2}\zeta^TS^{foc}(z)\zeta,\quad\quad H^{rot}(\zeta,z)=\frac{1}{2}\zeta^TS^{rot}(z)\zeta.


    **End-to-end map**

    A case of special interest occurs when the map is taken over a region symmetric about the solenoid midpoint at :math:`z=0`.
    In this case, we define :math:`z_0=-\lambda g`
    and :math:`z=\lambda g`.  Then we obtain:

    .. math::

        \Delta\phi=2\tan^{-1}\lambda,\quad\quad \lambda=\frac{L}{2g},

    where :math:`L` is the length of the field region over which the map is computed.  The map takes the form:

    .. math::

        M^0_{11}&=M^0_{33}=\cos(\beta\Delta\phi)+\frac{\lambda\sin(\beta\Delta\phi)}{\beta}, \\
        M^0_{12}&=M^0_{34}=g(1+\lambda^2)\frac{\sin(\beta\Delta\phi)}{\beta}, \\
        M^0_{21}&=M^0_{43}=\frac{2\beta\lambda\cos(\beta\Delta\phi)+(\lambda^2-\beta^2)\sin(\beta\Delta\phi)}{\beta g(1+\lambda^2)}, \\
        M^0_{22}&=M^0_{44}=\cos(\beta\Delta\phi)+\frac{\lambda\sin(\beta\Delta\phi)}{\beta}, \\
        M^0_{56}&=\frac{2\lambda g}{\beta_0^2\gamma_0^2},\quad\quad M^0_{55}=M^0_{66}=1.

    The expression for the rotation matrix :math:`R` is unchanged.  Note from the form of :math:`R` that the Larmor rotation angle :math:`\Omega` is given by:

    .. math::

        \Omega=b\Delta\phi=2b\tan^{-1}\lambda.

    This has a well defined limit as :math:`\lambda\rightarrow\infty`, which is:

    .. math::

        \lim_{\lambda\rightarrow\infty}\Omega= b\pi =\frac{\pi}{2}\frac{B_0g}{B\rho}.

    This is consistent with the Larmor angle we expect, namely:

    .. math::

        \Omega_{\infty}=\int_{-\infty}^{\infty}\alpha(z)dz=\frac{\pi}{2}\frac{B_0g}{B\rho}.

    .. _solenoid-softedge-solvable-ryne:

    **References**

    * R.D. Ryne, "Computational Methods in Accelerator Physics," USPAS Lecture Notes (2009)


Run
---

This example can be run as:

* **Python** script: ``python3 run_solenoid_softedge.py``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

       .. literalinclude:: run_solenoid_softedge.py
          :language: python3
          :caption: You can copy this file from ``examples/solenoid_softedge/run_solenoid_softedge_solvable.py``.


Analyze
-------

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_solenoid_softedge_solvable.py``

   .. literalinclude:: analysis_solenoid_softedge_solvable.py
      :language: python3
      :caption: You can copy this file from ``examples/solenoid_softedge/analysis_solenoid_softedge_solvable.py``.
