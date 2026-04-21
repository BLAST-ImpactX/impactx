.. _install-developers:
.. _install-build-cmake:
.. _building-cmake:
.. _building-cmake-intro:

Build from Source
=================

`CMake <https://cmake.org>`__ is our primary build system.
If you are new to CMake, we recommend starting with `this concise tutorial <https://hsf-training.github.io/hsf-training-cmake-webpage/>`__ from the HEP Software Foundation.
For those primarily interested in building the project, focus on these key sections: `1. Introduction <https://hsf-training.github.io/hsf-training-cmake-webpage/01-intro/index.html>`__, `2. Building with CMake <https://hsf-training.github.io/hsf-training-cmake-webpage/02-building/index.html>`__, and `9. Finding Packages <https://hsf-training.github.io/hsf-training-cmake-webpage/09-findingpackages/index.html>`__.


.. grid:: 2
   :gutter: 2

   .. grid-item-card:: :octicon:`alert` Build on HPC
      :link: install-hpc
      :link-type: ref

      Please refer to the :ref:`install-hpc` section.

   .. grid-item-card:: Install Dependencies
      :link: install-build-dependencies
      :link-type: ref

      Software dependencies of ImpactX.

   .. grid-item-card:: Build the Code
      :link: install-build-code
      :link-type: ref

      Configuration, compilation and install.

   .. grid-item-card:: Build Options
      :link: install-build-options
      :link-type: ref

      All build configuration options.


.. _install-build-dependencies:

Install Dependencies
--------------------

To begin, obtain a copy of the ImpactX source code:

.. code-block:: bash

   git clone https://github.com/BLAST-ImpactX/impactx.git $HOME/src/impactx
   cd $HOME/src/impactx

ImpactX relies on :ref:`widely-used third-party software <install-dependencies>`.
Below, you'll find instructions for installing these dependencies using various package managers.
To ensure compatibility, pick **one** package manager for your development workflows.

.. toctree::
   :hidden:

   dependencies


Install with conda-forge
^^^^^^^^^^^^^^^^^^^^^^^^

Conda provides a convenient way to install dependencies across Linux, macOS, and Windows platforms.

`conda-forge <https://conda-forge.org/download/>`__ is a community-led collection of recipes, build infrastructure and distributions for the conda package manager, offering cross-platform compatibility at the user level.

.. tip::

   We recommend to deactivate that conda self-activates its ``base`` environment.
   This `avoids interference with the system and other package managers <https://collegeville.github.io/CW20/WorkshopResources/WhitePapers/huebl-working-with-multiple-pkg-mgrs.pdf>`__.

   .. code-block:: bash

      conda config --set auto_activate_base false

   In order to make sure that the conda configuration uses ``conda-forge`` as the only channel, which will help avoid issues with blocked ``defaults`` or ``anaconda`` repositories, please set the following configurations:

   .. code-block:: bash

      conda config --add channels conda-forge
      conda config --set channel_priority strict

.. tab-set::

   .. tab-item:: With MPI (only Linux/macOS)

      .. code-block:: bash

         conda create -n impactx-cpu-mpich-dev -c conda-forge boost ccache cmake compilers git "openpmd-api=*=mpi_mpich*" packaging pytest pytest-benchmark python python-build numpy pandas quantiphy scipy setuptools yt "fftw=*=mpi_mpich*" pkg-config matplotlib mamba ninja mpich pip virtualenv vir-simd wheel
         conda activate impactx-cpu-mpich-dev

         # compile ImpactX with -DImpactX_MPI=ON
         # for pip, use: export IMPACTX_MPI=ON

   .. tab-item:: Without MPI

      .. code-block:: bash

         conda create -n impactx-cpu-dev -c conda-forge boost ccache cmake compilers git openpmd-api packaging pytest pytest-benchmark python python-build numpy pandas quantiphy scipy setuptools yt fftw pkg-config matplotlib mamba ninja pip virtualenv vir-simd wheel
         conda activate impactx-cpu-dev

         # compile ImpactX with -DImpactX_MPI=OFF
         # for pip, use: export IMPACTX_MPI=OFF

For OpenMP support, you will further need:

.. tab-set::

   .. tab-item:: Linux

      .. code-block:: bash

         conda install -c conda-forge libgomp

   .. tab-item:: macOS or Windows

      .. code-block:: bash

         conda install -c conda-forge llvm-openmp

For Nvidia CUDA GPU support, you will need to have `a recent CUDA driver installed <https://developer.nvidia.com/cuda-downloads>`__ or you can lower the CUDA version of `the Nvidia cuda package <https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html#conda-installation>`__ and `conda-forge to match your drivers <https://docs.cupy.dev/en/stable/install.html#install-cupy-from-conda-forge>`__ and then add these packages:

.. code-block:: bash

   conda install -c nvidia -c conda-forge cuda cuda-nvtx-dev cupy

For the ImpactX browser/Jupyter dashboard dependencies, install from the ImpactX source directory:

.. code-block:: bash

   python3 -m pip install -r src/python/impactx/dashboard/requirements.txt


Install with Spack
^^^^^^^^^^^^^^^^^^

Spack provides another option for installing dependencies on Linux and macOS systems.

`Spack <https://spack.readthedocs.io>`__ is a flexible, user-level package manager designed primarily for Linux, with growing support for macOS and planned future support for Windows.

To begin, download a `WarpX Spack desktop development environment <https://github.com/BLAST-WarpX/warpx/blob/development/Tools/machines/desktop>`__ configuration (which also works for ImpactX).
For most desktop development work, we recommend using the OpenMP environment for CPUs, unless you have a supported GPU device.

* **Debian/Ubuntu** Linux:

  * OpenMP: ``system=ubuntu; compute=openmp`` (CPUs)
  * CUDA: ``system=ubuntu; compute=cuda`` (Nvidia GPUs)
  * ROCm: ``system=ubuntu; compute=rocm`` (AMD GPUs)
  * SYCL: *todo* (Intel GPUs)
* **macOS**: first, prepare with ``brew install gpg2; brew install gcc``

  * OpenMP: ``system=macos; compute=openmp``

If you already `installed Spack <https://spack.io>`__, we recommend to activate its `binary caches <https://spack.io/spack-binary-packages/>`__ for faster builds:

.. code-block:: bash

   spack mirror add rolling https://binaries.spack.io/develop
   spack buildcache keys --install --trust

Now install the ImpactX dependencies in a new ImpactX development environment:

.. code-block:: bash

   # download environment file
   curl -sLO https://raw.githubusercontent.com/BLAST-WarpX/warpx/development/Tools/machines/desktop/spack-${system}-${compute}.yaml

   # create new development environment
   spack env create impactx-${compute}-dev spack-${system}-${compute}.yaml
   spack env activate impactx-${compute}-dev

   # installation
   spack install
   python3 -m pip install jupyter matplotlib numpy openpmd-api openpmd-viewer pandas quantiphy scipy virtualenv yt

In new terminal sessions, re-activate the environment with

.. code-block:: bash

   spack env activate impactx-openmp-dev

again.
Replace ``openmp`` with the equivalent you chose.

Compile ImpactX with ``-DImpactX_MPI=ON``.
For ``pip``, use ``export IMPACTX_MPI=ON``.


Install with Brew
^^^^^^^^^^^^^^^^^

Brew can be used to install dependencies on Linux and macOS.

`Homebrew (Brew) <https://brew.sh>`__ is a user-level package manager primarily for `Apple macOS <https://en.wikipedia.org/wiki/MacOS>`__, but also supports Linux.

.. code-block:: bash

   brew update
   brew tap openpmd/openpmd
   brew install adios2      # for openPMD
   brew install ccache
   brew install cmake
   brew install fftw        # for IGF space charge, CSR
   brew install git
   brew install hdf5-mpi    # for openPMD
   brew install libomp
   brew unlink gcc
   brew link --force libomp
   brew install pkg-config  # for fftw
   brew install open-mpi
   brew install openblas
   brew install openpmd-api # for openPMD

Compile ImpactX with ``-DImpactX_MPI=ON``.
For ``pip``, use ``export IMPACTX_MPI=ON``.


Install with APT
^^^^^^^^^^^^^^^^

The `Advanced Package Tool (APT) <https://en.wikipedia.org/wiki/APT_(software)>`__ is a system-level package manager on Debian-based Linux distributions, including Ubuntu.

.. tab-set::

   .. tab-item:: With MPI (only Linux/macOS)

      .. code-block:: bash

         sudo apt update
         sudo apt install build-essential ccache cmake g++ git libfftw3-mpi-dev libfftw3-dev libhdf5-openmpi-dev libopenmpi-dev pkg-config python3 python3-dev python3-matplotlib python3-mpi4py python3-numpy python3-pandas python3-pip python3-scipy python3-venv

         # optional:
         # for CUDA, either install
         #   https://developer.nvidia.com/cuda-downloads (preferred)
         # or, if your Debian/Ubuntu is new enough, use the packages
         #   sudo apt install nvidia-cuda-dev libcub-dev

         # compile ImpactX with -DImpactX_MPI=ON
         # for pip, use: export IMPACTX_MPI=ON

   .. tab-item:: Without MPI

      .. code-block:: bash

         sudo apt update
         sudo apt install build-essential ccache cmake g++ git libfftw3-dev libhdf5-dev pkg-config python3 python3-dev python3-matplotlib python3-numpy python3-pandas python3-pip python3-scipy python3-venv

         # optional:
         # for CUDA, either install
         #   https://developer.nvidia.com/cuda-downloads (preferred)
         # or, if your Debian/Ubuntu is new enough, use the packages
         #   sudo apt install nvidia-cuda-dev libcub-dev

         # compile ImpactX with -DImpactX_MPI=OFF
         # for pip, use: export IMPACTX_MPI=OFF


.. _install-build-code:

Build the Code
--------------

.. _build-the-executable-with-cmake:

Build the Executable with CMake
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To build ImpactX from the source directory, execute these commands:

.. code-block:: bash

   # Configure the build system
   # Additional options available, such as:
   #   -DImpactX_PYTHON=ON
   #   -DCMAKE_INSTALL_PREFIX=$HOME/sw/impactx
   cmake -S . -B build

   # Compile using four parallel threads
   cmake --build build -j 4

**That's it!**
The ImpactX binary is now available in ``build/bin/`` and is ready to :ref:`run <usage_run>` with any :ref:`example input file <usage-examples>`.
You can either run the binary directly from this location or copy it to another directory.

For a system-wide installation, use the following command:

.. code-block:: bash

   # for default install paths, you will need administrator rights, e.g. with sudo:
   cmake --build build --target install

You can inspect and modify build options after running ``cmake -S . -B build`` with either

.. code-block:: bash

   ccmake build

or by adding arguments with ``-D<OPTION>=<VALUE>`` to the first CMake call.
For example, this builds ImpactX with Python bindings and Nvidia GPU (CUDA) support:

.. code-block:: bash

   cmake -S . -B build -DImpactX_PYTHON=ON -DImpactX_COMPUTE=CUDA


An executable ImpactX application binary will be created in ``build/bin/``.
Additionally, a `symbolic link <https://en.wikipedia.org/wiki/Symbolic_link>`__ named ``impactx`` can be found in that directory, which points to the last built ImpactX executable.

More details on running simulations are in the section :ref:`Run ImpactX <usage_run>`.
Alternatively, read on and also build our Python interface.


.. _install-build-python-cmake:

Build the Python Interface with CMake
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note::

   First, ensure your Python development environment is up-to-date:

   .. code-block:: bash

      python3 -m pip install -U pip
      python3 -m pip install -U build packaging setuptools[core] wheel pytest pytest-benchmark
      python3 -m pip install -U cmake
      python3 -m pip install -U -r examples/requirements.txt

To build the Python bindings, configure ImpactX to generate a library and install it using our ``pip_install`` *CMake target*:

.. code-block:: bash

   # Configure with Python support enabled
   cmake -S . -B build_py -DImpactX_PYTHON=ON

   # Build and install the Python package
   cmake --build build_py --target pip_install -j 4

**That's it!**
You can now :ref:`run a first ImpactX Python script <usage-python>` from our :ref:`examples <usage-examples>`.

Developers could now change the ImpactX source code and then call the build line again to refresh the Python installation.

.. tip::

   If you do *not* develop with :ref:`a user-level package manager <install-dependencies>`, e.g., because you rely on a HPC system's environment modules, then consider to set up a virtual environment via `Python venv <https://docs.python.org/3/library/venv.html>`__.
   Otherwise, without a virtual environment, you likely need to add the CMake option ``-DPY_PIP_INSTALL_OPTIONS="--user"``.


.. _install-build-python-pip:

Build the Python Interface with pip
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section is relevant for Python package management, mainly for maintainers or people that rather like to interact only with ``pip``.

One can build and install ``impactx`` from the root of the ImpactX source tree:

.. code-block:: bash

   python3 -m pip wheel -v .
   python3 -m pip install impactx*whl

This will call the CMake logic above implicitly.


.. _install-build-options:

Build Options
-------------

Configure your Compiler
^^^^^^^^^^^^^^^^^^^^^^^

To use a specific compiler instead of the system default, set the appropriate environment variables.
For instance, to use Clang/LLVM:

.. code-block:: bash

   export CC=$(which clang)
   export CXX=$(which clang++)

For CUDA development, specify both the CUDA compiler and the host C++ compiler:

.. code-block:: bash

   export CUDACXX=$(which nvcc)
   export CUDAHOSTCXX=$(which clang++)

We also support adding `additional compiler flags via environment variables <https://cmake.org/cmake/help/latest/manual/cmake-language.7.html#cmake-language-environment-variables>`__ such as `CXXFLAGS <https://cmake.org/cmake/help/latest/envvar/CXXFLAGS.html>`__/`LDFLAGS <https://cmake.org/cmake/help/latest/envvar/LDFLAGS.html>`__:

.. code-block:: bash

   # example: treat all compiler warnings as errors
   export CXXFLAGS="-Werror"

.. note::

   Please clean your build directory with ``rm -rf build/`` after changing the compiler.
   Now call ``cmake -S . -B build`` (+ further options) again to re-initialize the build configuration.


CMake
^^^^^

=============================== ============================================ ===========================================================
CMake Option                    Default & Values                             Description
=============================== ============================================ ===========================================================
``BUILD_TESTING``               **ON**/OFF                                   `Build tests <https://cmake.org/cmake/help/latest/module/CTest.html>`__
``CMAKE_BUILD_TYPE``            RelWithDebInfo/**Release**/Debug             `Type of build, symbols & optimizations <https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html>`__
``CMAKE_INSTALL_PREFIX``        system-dependent path                        `Install path prefix <https://cmake.org/cmake/help/latest/variable/CMAKE_INSTALL_PREFIX.html>`__
``CMAKE_VERBOSE_MAKEFILE``      ON/**OFF**                                   `Print all compiler commands to the terminal during build <https://cmake.org/cmake/help/latest/variable/CMAKE_VERBOSE_MAKEFILE.html>`__
``ImpactX_APP``                 **ON**/OFF                                   Build the ImpactX executable application
``ImpactX_COMPUTE``             NOACC/**OMP**/CUDA/SYCL/HIP                  On-node, accelerated computing backend
``ImpactX_FASTMATH``            ON/**OFF**                                   Enable fast-math optimizations
``ImpactX_FFT``                 ON/**OFF**                                   FFT-based solvers (IGF space charge, CSR)
``ImpactX_IPO``                 ON/**OFF**                                   Compile ImpactX with interprocedural optimization (aka LTO)
``ImpactX_MPI``                 **ON**/OFF                                   Multi-node support (message-passing)
``ImpactX_MPI_THREAD_MULTIPLE`` **ON**/OFF                                   MPI thread-multiple support, i.e. for ``async_io``
``ImpactX_OPENPMD``             **ON**/OFF                                   openPMD I/O (HDF5, ADIOS)
``ImpactX_PRECISION``           SINGLE/**DOUBLE**                            Floating point precision (single/double)
``ImpactX_PYTHON``              ON/**OFF**                                   Python bindings
``ImpactX_SIMD``                ON/**OFF**                                   CPU SIMD acceleration (requires ``vir-simd``)
``Python_EXECUTABLE``           (newest found)                               Path to Python executable
``PY_PIP_OPTIONS``              ``-v``                                       Additional options for ``pip``, e.g., ``-vvv;-q``
``PY_PIP_INSTALL_OPTIONS``                                                   Additional options for ``pip install``, e.g., ``--user;-q``
=============================== ============================================ ===========================================================

ImpactX can be configured in further detail with options from AMReX, which are documented in the AMReX manual:

* `general AMReX build options <https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#customization-options>`__
* `GPU-specific options <https://amrex-codes.github.io/amrex/docs_html/GPU.html#building-gpu-support>`__.

**Developers** might be interested in additional options that control dependencies of ImpactX.
By default, the most important dependencies of ImpactX are automatically downloaded for convenience:

============================= ============================================== ===========================================================
CMake Option                  Default & Values                               Description
============================= ============================================== ===========================================================
``BUILD_SHARED_LIBS``         ON/**OFF**                                     `Build shared libraries for dependencies <https://cmake.org/cmake/help/latest/variable/BUILD_SHARED_LIBS.html>`__
``ImpactX_CCACHE``            **ON**/OFF                                     Search and use CCache to speed up rebuilds.
``ImpactX_UNITY_BUILD``       ON/**OFF**                                     ImpactX library as unity build (single TU)
``ImpactX_ablastr_src``       *None*                                         Path to ABLASTR source directory (preferred if set)
``ImpactX_ablastr_repo``      ``https://github.com/BLAST-WarpX/warpx.git``   Repository URI to pull and build ABLASTR from
``ImpactX_ablastr_branch``    *we set and maintain a compatible commit*      Repository branch for ``ImpactX_ablastr_repo``
``ImpactX_ablastr_internal``  **ON**/OFF                                     Needs a pre-installed ABLASTR/WarpX library if set to ``OFF``
``ImpactX_amrex_src``         *None*                                         Path to AMReX source directory (preferred if set)
``ImpactX_amrex_repo``        ``https://github.com/AMReX-Codes/amrex.git``   Repository URI to pull and build AMReX from
``ImpactX_amrex_branch``      *we set and maintain a compatible commit*      Repository branch for ``ImpactX_amrex_repo``
``ImpactX_amrex_internal``    **ON**/OFF                                     Needs a pre-installed AMReX library if set to ``OFF``
``ImpactX_openpmd_src``       *None*                                         Path to openPMD-api source directory (preferred if set)
``ImpactX_openpmd_repo``      ``https://github.com/openPMD/openPMD-api.git`` Repository URI to pull and build openPMD-api from
``ImpactX_openpmd_branch``    ``0.17.0``                                     Repository branch for ``ImpactX_openpmd_repo``
``ImpactX_openpmd_internal``  **ON**/OFF                                     Needs a pre-installed openPMD-api library if set to ``OFF``
``ImpactX_pyamrex_src``       *None*                                         Path to pyAMReX source directory (preferred if set)
``ImpactX_pyamrex_repo``      ``https://github.com/AMReX-Codes/pyamrex.git`` Repository URI to pull and build pyAMReX from
``ImpactX_pyamrex_branch``    *we set and maintain a compatible commit*      Repository branch for ``ImpactX_pyamrex_repo``
``ImpactX_pyamrex_internal``  **ON**/OFF                                     Needs a pre-installed pyAMReX module if set to ``OFF``
``ImpactX_PYTHON_IPO``        **ON**/OFF                                     Build Python w/ interprocedural/link optimization (IPO/LTO)
``ImpactX_pybind11_src``      *None*                                         Path to pybind11 source directory (preferred if set)
``ImpactX_pybind11_repo``     ``https://github.com/pybind/pybind11.git``     Repository URI to pull and build pybind11 from
``ImpactX_pybind11_branch``   *we set and maintain a compatible commit*      Repository branch for ``ImpactX_pybind11_repo``
``ImpactX_pybind11_internal`` **ON**/OFF                                     Needs a pre-installed pybind11 library if set to ``OFF``
``ImpactX_TEST_CLEANUP``      ON/**OFF**                                     Clean up automated test directories
============================= ============================================== ===========================================================

For example, one can also build against a local AMReX copy.
Assuming AMReX' source is located in ``$HOME/src/amrex``, add the ``cmake`` argument ``-DImpactX_amrex_src=$HOME/src/amrex``.
Relative paths are also supported, e.g. ``-DImpactX_amrex_src=../amrex``.

Or build against an AMReX feature branch of a colleague.
Assuming your colleague pushed AMReX to ``https://github.com/WeiqunZhang/amrex/`` in a branch ``new-feature`` then pass to ``cmake`` the arguments: ``-DImpactX_amrex_repo=https://github.com/WeiqunZhang/amrex.git -DImpactX_amrex_branch=new-feature``.

If you want to develop against local versions of ABLASTR (from WarpX) and AMReX at the same time, pass for instance ``-DImpactX_ablastr_src=$HOME/src/warpx -DImpactX_amrex_src=$HOME/src/amrex``.

You can speed up the install further if you pre-install these dependencies, e.g. with a package manager.
Set ``-DImpactX_<dependency-name>_internal=OFF`` and add installation prefix of the dependency to the environment variable `CMAKE_PREFIX_PATH <https://cmake.org/cmake/help/latest/envvar/CMAKE_PREFIX_PATH.html>`__.
Please see the :ref:`introduction to CMake <install-build-cmake>` if this sounds new to you.

If you re-compile often, consider installing the `Ninja <https://github.com/ninja-build/ninja/wiki/Pre-built-Ninja-packages>`__ build system.
Pass ``-G Ninja`` to the CMake configuration call to speed up parallel compiles.


pip
^^^

Environment variables can be used to control the build step:

============================= ============================================ ================================================================
Environment Variable          Default & Values                             Description
============================= ============================================ ================================================================
``IMPACTX_COMPUTE``           NOACC/**OMP**/CUDA/SYCL/HIP                  On-node, accelerated computing backend
``IMPACTX_MPI``               ON/**OFF**                                   Multi-node support (message-passing)
``IMPACTX_PRECISION``         SINGLE/**DOUBLE**                            Floating point precision (single/double)
``IMPACTX_FFT``               ON/**OFF**                                   FFT-based solvers (IGF space charge, CSR)
``IMPACTX_SIMD``              ON/**OFF**                                   CPU SIMD acceleration (requires ``vir-simd``)
``IMPACTX_BUILD_SHARED_LIBS`` ON/**OFF**                                   Build shared libraries for dependencies
``IMPACTX_AMREX_SRC``         *None*                                       Absolute path to AMReX source directory (preferred if set)
``IMPACTX_AMREX_REPO``        *None (uses cmake default)*                  Repository URI to pull and build AMReX from
``IMPACTX_AMREX_BRANCH``      *None (uses cmake default)*                  Repository branch for ``IMPACTX_AMREX_REPO``
``IMPACTX_AMREX_INTERNAL``    **ON**/OFF                                   Needs a pre-installed AMReX library if set to ``OFF``
``IMPACTX_PYAMREX_SRC``       *None*                                       Absolute path to pyAMReX source directory (preferred if set)
``IMPACTX_PYAMREX_REPO``      *None (uses cmake default)*                  Repository URI to pull and build pyAMReX from
``IMPACTX_PYAMREX_BRANCH``    *None (uses cmake default)*                  Repository branch for ``IMPACTX_PYAMREX_REPO``
``IMPACTX_PYAMREX_INTERNAL``  **ON**/OFF                                   Needs a pre-installed pyAMReX library if set to ``OFF``
``PYIMPACTX_LIBDIR``          *None*                                       If set, search for pre-built ImpactX C++ libraries (see below)
============================= ============================================ ================================================================

Note that we currently change the ``IMPACTX_MPI`` default intentionally to ``OFF``, to simplify a first install from source.

Additional CMake options can be passed through ``pip`` using variables of the form ``IMPACTX_CMAKE_<NAME>=<VALUE>``, which will be forwarded as ``-D<NAME>=<VALUE>`` to CMake.
Use this to control, e.g., ABLASTR, openPMD-api, or pybind11 source/repo/branch selection when building via ``pip``.

Some hints and workflows follow.
Developers that want to test a change of the source code but did not change the ``impactx`` version number can force a reinstall via:

.. code-block:: bash

   python3 -m pip install --force-reinstall --no-deps -v .

Some developers like to code directly against a local copy of AMReX, changing both code-bases at a time:

.. code-block:: bash

   IMPACTX_AMREX_SRC=$PWD/../amrex python3 -m pip install --force-reinstall --no-deps -v .

Additional environment control as common for CMake (:ref:`see above <install-build-cmake>`) can be set as well, e.g. ``CC``, ``CXX``, and ``CMAKE_PREFIX_PATH`` hints.
So another sophisticated example might be: use Clang as the compiler, build with local source copies of ABLASTR and AMReX, support the FFT-based solvers, MPI and openPMD, and hint a parallel HDF5 installation in ``$HOME/sw/hdf5-parallel-1.10.4``:

.. code-block:: bash

   CC=$(which clang) CXX=$(which clang++) IMPACTX_CMAKE_ImpactX_ablastr_src=$PWD/../warpx IMPACTX_AMREX_SRC=$PWD/../amrex IMPACTX_FFT=ON IMPACTX_MPI=ON CMAKE_PREFIX_PATH=$HOME/sw/hdf5-parallel-1.10.4:$CMAKE_PREFIX_PATH python3 -m pip install --force-reinstall --no-deps -v .

Here we wrote this all in one line, but one can also set all environment variables in a development environment and keep the pip call nice and short as in the beginning.
Note that you need to use absolute paths for external source trees, because pip builds in a temporary directory, e.g. ``export IMPACTX_AMREX_SRC=$HOME/src/amrex``.

All of this can also be run from CMake.
This is the workflow most developers will prefer as it allows rapid re-compiles:

.. code-block:: bash

   # build ImpactX executables and libraries
   cmake -S . -B build_py -DImpactX_PYTHON=ON

   # build & install Python only
   cmake --build build_py -j 4 --target pip_install

There is also a ``--target pip_install_nodeps`` option that skips pip-based dependency checks.

Last but not least, you can uninstall ``impactx`` as usual with:

.. code-block:: bash

   python3 -m pip uninstall impactx
