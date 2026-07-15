.. _glossary:

Glossary
========

In daily communication, we tend to abbreviate a lot of terms.
It is important to us to make it easy to interact with the ImpactX community and thus, this list shall help to clarify often used terms.

Abbreviations
-------------

* **2FA:** `Two-factor-authentication <https://en.wikipedia.org/wiki/Multi-factor_authentication>`__
* **ABLASTR:** Accelerated BLAST Recipes, a library shared between the BLAST codes WarpX and ImpactX, providing, e.g., the space charge and CSR solvers
* **ADIOS2:** `The Adaptable Input Output System version 2 <https://adios2.readthedocs.io>`__, an I/O backend for openPMD output (``bp`` files)
* **AI:** artificial intelligence. ImpactX provides user-friendly interfaces suitable for AI/ML workflows, see, e.g., :ref:`AI-assisted input design <ai_input_design>`
* **AMR:** adaptive mesh-refinement, used in our space charge solvers to resolve the particle beam where it needs it most
* **API:** `application programming interface <https://en.wikipedia.org/wiki/API>`__, e.g., the :ref:`Python interface <usage-python>` of ImpactX
* **BC:** boundary condition (of a simulation), e.g., the longitudinal :ref:`particle boundary condition <running-cpp-parameters-particle-bc>`
* **BLAST:** `Beam, Plasma & Accelerator Simulation Toolkit <https://blast.lbl.gov>`__, the open source code suite that ImpactX is part of
* **CI:** continuous integration, automated tests that we perform before a proposed code-change is accepted; see PR
* **CPU:** `central processing unit <https://en.wikipedia.org/wiki/Central_processing_unit>`__, we usually mean a socket or generally the host-side of a computer (compared to the computer's "accelerator", i.e. GPU)
* **CS:** Courant-Snyder; see *Twiss parameters* in the terms below
* **CSR:** :ref:`coherent synchrotron radiation <running-cpp-parameters-collective-csr>`, a collective effect from radiation generated when a charged particle beam is bent
* **DOE:** `The United States Department of Energy <https://en.wikipedia.org/wiki/United_States_Department_of_Energy>`__, the largest sponsor of national laboratory research in the United States of America
* **DP:** double precision (64-bit floating point), the default ``ImpactX_PRECISION``; see also SP
* **FFT:** `fast Fourier transform <https://en.wikipedia.org/wiki/Fast_Fourier_transform>`__, e.g., used by the IGF space charge solver and the CSR calculation (``-DImpactX_FFT=ON``)
* **FODO:** a periodic focusing lattice cell consisting of a focusing (F) quadrupole, a drift (O), a defocusing (D) quadrupole and another drift (O). The "hello world" of our :ref:`examples <usage-examples>`
* **GPU:** originally graphics processing unit, now used for fast `general purpose computing (GPGPU) <https://en.wikipedia.org/wiki/Graphics_processing_unit#Stream_processing_and_general_purpose_GPUs_(GPGPU)>`__. Also called (hardware) accelerator of a computer.
* **HDF5:** `Hierarchical Data Format 5 <https://www.hdfgroup.org/solutions/hdf5/>`__, an I/O backend for openPMD output (``h5`` files)
* **HPC:** high-performance computing, e.g., on large computing clusters and supercomputers
* **IGF:** integrated Green function, the method used by the default (``fft``) :ref:`space charge Poisson solver <running-cpp-parameters-collective-spacecharge>`
* **IO:** input/output, usually files and/or data
* **IOTA:** the `Integrable Optics Test Accelerator <https://fast.fnal.gov/>`__ at Fermilab. ImpactX implements its nonlinear magnetic insert as a lattice element
* **IPO:** `interprocedural optimization <https://en.wikipedia.org/wiki/Interprocedural_optimization>`__, a collection of compiler optimization techniques that analyze the whole code to avoid duplicate calculations and optimize performance.  Also called link-time optimization (LTO), see below.
* **ISR:** :ref:`incoherent synchrotron radiation <running-cpp-parameters-collective-isr>`, causing radiation reaction (mean energy loss) of the bunch in bending elements, often relevant at high energies
* **JSON:** `JavaScript Object Notation <https://en.wikipedia.org/wiki/JSON>`__, a simple text format.
* **K-V:** Kapchinskij-Vladimirskij distribution, an :ref:`initial beam distribution <running-cpp-parameters-particle>` that is uniform on the surface of a 4D transverse phase space ellipsoid
* **LBNL:** `Lawrence Berkeley National Laboratory <https://www.lbl.gov>`__, the U.S. DOE national laboratory leading the development of ImpactX
* **LLM:** large language model. See our documentation on :ref:`AI-assisted input design <ai_input_design>` and :ref:`developing ImpactX with LLMs <developers-llm>`
* **LTO:** `link-time optimization <https://en.wikipedia.org/wiki/Interprocedural_optimization#WPO_and_LTO>`__, program optimizations for file-by-file compilation that optimize object files before linking them together to an executable.  Also called interprocedural optimization (IPO), see above.
* **MAD-X:** `Methodical Accelerator Design <https://mad.web.cern.ch/mad/>`__ (version X), a CERN code for accelerator design. ImpactX often follows MAD-X conventions for element strengths and can import MAD-X lattice files
* **MCP:** `Model Context Protocol <https://modelcontextprotocol.io>`__, an open standard to provide context (e.g., the ImpactX documentation) to LLMs
* **ML:** machine learning, e.g., trained neural networks can be used as ML surrogate models for lattice elements, sections, or collective effects
* **MLMG:** multi-level multigrid, the iterative solver used by the ``multigrid`` :ref:`space charge Poisson solver <running-cpp-parameters-collective-spacecharge>` option
* **MPI:** `message passing interface <https://en.wikipedia.org/wiki/Message_Passing_Interface>`__, used for parallelization over multiple nodes and/or multiple GPUs
* **MR:** mesh-refinement, see AMR
* **NERSC:** `National Energy Research Scientific Computing Center <https://www.nersc.gov/>`__, a supercomputing center located in Berkeley, CA (USA)
* **OTP:** `One-Time-Password <https://en.wikipedia.org/wiki/One-time_password>`__, see 2FA
* **PALS:** `Particle Accelerator Lattice Standard <https://github.com/pals-project/pals>`__, a community effort toward a common lattice description file format, which ImpactX can read
* **PR:** github pull request, a proposed change to the ImpactX code base
* **RF:** radio-frequency, e.g., an RF cavity element for particle acceleration and bunching
* **RMS:** `root mean square <https://en.wikipedia.org/wiki/Root_mean_square>`__, e.g., in RMS beam sizes and RMS emittances calculated from the beam moments
* **SIMD:** `single instruction, multiple data <https://en.wikipedia.org/wiki/Single_instruction,_multiple_data>`__, CPU vectorization used to speed up particle tracking (``-DImpactX_SIMD=ON``)
* **SP:** single precision (32-bit floating point), see also DP

Terms
-----

* **accelerator:** depending on context, either a *particle accelerator* in physics or a *hardware accelerator* (i.e. GPU) in computing
* **AMReX:** `C++ library for block-structured adaptive mesh-refinement <https://amrex-codes.github.io>`__, a primary dependency of ImpactX
* **aperture:** a thin collimator element applying a transverse boundary to the beam; particles outside are lost (see *scraping*)
* **beam monitor:** a zero-length lattice element that writes the particle beam to :ref:`openPMD output <dataanalysis-monitor>` whenever the beam passes it
* **bunch:** a longitudinally confined group of particles, also referred to as beam; the particle ensemble that ImpactX tracks through the lattice
* **chromatic effects:** the dependence of the transverse focusing optics on the particle energy. Chromatic element models retain the exact nonlinear dependence on the energy variable :math:`p_t` (see :ref:`assumptions <theory-assumptions>`)
* **collective effects:** effects that depend on the particle beam as a whole, arising from the interaction of the beam particles with each other through their self-fields and radiation, in contrast to single-particle optics in external fields. ImpactX :ref:`applies <running-cpp-parameters-collective>` collective effects such as *space charge* and the CSR *wakefield* between the *slices* of lattice elements
* **covariance matrix:** the 6x6 matrix :math:`\Sigma_{ij}=\langle\zeta_i\zeta_j\rangle` of :ref:`second moments of the phase space coordinates <theory-collective-beam-distribution-input>`; used to characterize the initial beam distribution and as the tracked state in *envelope tracking*
* **dashboard:** the :ref:`browser-based graphical user interface <usage-dashboard>` of ImpactX, which can also run in Jupyter notebooks
* **dispersion:** linear dependence of the transverse particle orbit on the particle energy deviation (important within bending dipoles), reported in the :ref:`reduced beam characteristics <dataanalysis-beam-characteristics>`
* **distribution:** the initial phase space density of the particle beam; ImpactX :ref:`implements <running-cpp-parameters-particle>` Gaussian, K-V, waterbag, Kurth, semi-Gaussian, triangle and thermal distributions, among others
* **eigenemittances:** the three invariant emittances of the (possibly coupled) 6D beam distribution, also known as mode emittances; they coincide with the RMS emittances when inter-plane correlations vanish
* **emittance:** a measure of the phase space area occupied by the beam in a coordinate plane; ImpactX reports unnormalized (geometric, 1-RMS) emittances :math:`\epsilon` and normalized emittances :math:`\epsilon_\mathrm{n} = (\beta\gamma)_\mathrm{ref} \cdot \epsilon` (see :ref:`beam distribution input <theory-collective-beam-distribution-input>`)
* **envelope tracking:** a simplified :ref:`tracking mode <usage_run-tracking-mode>` that evolves the beam envelope (6x6 covariance matrix) through linearized transport maps, for rapid parameter scans
* **hard-edge element:** an ideal (thick) element model in which the :ref:`fields are independent of s <theory-softedge-elements>` inside the element and terminate abruptly at its ends; see *soft-edge element*
* **lattice:** the ordered sequence of :ref:`beamline elements <running-cpp-parameters-lattice>` (drifts, magnets, cavities, ...) that the beam is tracked through; it can be repeated periodically, e.g., for rings or channels
* **linac:** a linear particle accelerator
* **LUMI:** a `European pre-exascale supercomputer at CSC <https://www.lumi-supercomputer.eu>`__ (Finland)
* **macroparticle:** a weighted simulation particle representing many physical beam particles
* **magnetic rigidity:** :math:`B\rho = p/q`, the momentum per unit charge of the reference particle; element strengths (e.g., the quadrupole strength :math:`k`) are normalized by it
* **openPMD:** `Open Standard for Particle-Mesh Data Files <https://www.openPMD.org>`__, a community meta-data project for scientific data; the output format of the *beam monitor* and input for the particle *source* element
* **Perlmutter:** a Berkeley Lab nobel laureate and a `Pre-Exascale supercomputer at NERSC <https://www.nersc.gov/systems/perlmutter/>`__
* **pyAMReX:** the `Python binding of AMReX <https://pyamrex.readthedocs.io>`__, used to access particle and field data from Python without copying
* **Python:** a popular scripted `programming language <https://www.python.org>`__ and flexible interface language to ImpactX
* **reference particle:** the particle that follows the reference trajectory; it is described by 8 phase space variables in a :ref:`global lab coordinate system <theory-coordinates-and-units>` and the beam particle coordinates are deviations from it
* **reference trajectory:** the :ref:`nominal machine orbit <theory-concepts>`, i.e., the ideal single-particle orbit used as part of the optics design; it relates the local beam coordinate system to the global lab coordinate system
* **ring:** a circular accelerator, e.g., a storage ring or synchrotron. Usually modeled by repeating the lattice over many turns
* **s:** the path length along the reference trajectory, used as the :ref:`independent dynamical variable <theory-concepts>` of tracking (instead of time)
* **scraping:** the removal of particles that are lost from the beam, e.g., when hitting an *aperture* or leaving the simulation domain. Lost particles are collected and can be written to openPMD output (``particles_lost``)
* **slice:** thick lattice elements are sliced into ``nslice`` segments. Collective effects such as space charge kicks are applied between slices
* **soft-edge element:** an element modeled with an :math:`s`-dependent, user-provided :ref:`on-axis field profile <theory-softedge-elements>`, e.g., RF cavities, soft-edge solenoids and soft-edge quadrupoles
* **space charge:** the collective effect of the beam's own electrostatic self-field. :ref:`Calculated <running-cpp-parameters-collective-spacecharge>` by solving the Poisson equation in the bunch frame, between slices of the lattice elements
* **spin tracking:** :ref:`updating the particle spin 3-vector <running-cpp-parameters-spin>` in the electromagnetic fields of each element, using methods based on the Thomas-BMT equation; the mean of the particle spin vectors is the *polarization vector*
* **symplectic:** the defining property of Hamiltonian (transfer) maps, which preserve the phase space structure; :ref:`particle tracking through lattice optics <theory-assumptions>` in ImpactX is symplectic
* **tracking modes:** ImpactX can :ref:`track <usage_run-tracking-mode>` individual particles (full dynamics), the beam envelope (linearized optics), or only the reference orbit (for early design)
* **transfer map:** the map that advances the phase space coordinates through a lattice element; to lowest order a 6x6 matrix (linear transfer map), which users can also supply directly (``linear_map``)
* **Twiss parameters:** the Courant-Snyder parameters :math:`\alpha`, :math:`\beta` and :math:`\gamma` describing the :ref:`phase space ellipse <theory-collective-beam-distribution-input>` of the beam, related by :math:`\gamma\beta - \alpha^2 = 1`
* **wakefield:** a collective self-field (other than space charge) that acts back on the beam, e.g., the CSR wakefield generated in bending elements
* **WarpX:** the `electromagnetic and electrostatic particle-in-cell code <https://warpx.readthedocs.io>`__ of BLAST, e.g., used to model source, injector, or plasma-based accelerator elements. ImpactX shares the ABLASTR library with WarpX and can exchange particle beams with WarpX through openPMD or Python

For plasma, laser and fusion related terminology, please also see the `WarpX glossary <https://warpx.readthedocs.io/en/latest/glossary.html>`__.
