.. _theory-softedge-elements:

Models of Soft-Edge Elements
============================

Similarly to other beam dynamics codes, most beamline elements are modeled by using an ideal, hard-edge approximation for the potentials and fields.  In particular, this model assumes that the fields are independent of the longitudinal path-length coordinate :math:`s`.
However, ImpactX also supports several additional soft-edge element models.  The models include:

* **RFCavity** - soft-edge model of an RF cavity using on-axis field data :math:`(z,E_z(z))`
* **SoftSolenoid** - soft-edge model of a solenoid using on-axis field data :math:`(z,B_z(z))`
* **SoftQuadrupole** - soft-edge model of a solenoid using on-axis field gradient data :math:`(z,\partial B_y(z)/\partial x)`

For these elements, the user may specify the on-axis field or field gradient in one of two forms:

* using tabulated on-axis data for the field or gradient (depending on the element type)
* using a set of precomputed Fourier coefficients obtained from the on-axis data

The vector potential off-axis is then determined from Maxwell's equations.  In either of the above cases, the field is represented internally using a set of Fourier coefficients.  If :math:`g(z)` denotes the on-axis field profile, the Fourier coefficients are defined such that:

.. math::

   \begin{align}
        g(z) &= \frac{c_0}{2} + \sum_{j=1}^{\rm ncoef}c_j\cos\left(\frac{2\pi j(z-z_{\rm mid})}{L}\right)+\sum_{j=1}^{\rm ncoef}s_j\sin\left(\frac{2\pi j(z-z_{\rm mid})}{L}\right) \\
        z_{\rm mid}&=(z_{\rm max}+z_{\rm min})/2,\quad\quad L=z_{\rm max}-z_{\rm min}.
   \end{align}

An example illustrating the construction of the Fourier coefficients from on-axis data appears here:  ADD LINK HERE

The internal representation of the field is scaled based on the user-specified inputs.  In particular, the longitudinal coordinate is scaled to coincide with the element length.  Likewise, the field strength is scaled based on input parameters (called ``escale``, ``bscale``, or ``gscale`` for the three element types above).

Examples demonstrating the use of soft-edge elements include:  ADD LINK HERE
