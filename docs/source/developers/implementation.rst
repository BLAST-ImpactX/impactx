.. _developers-implementation:

Implementation Details
======================

.. note::

   TODO :-)

Kernel Launches
---------------

Compute kernels are launched through AMReX's portable constructs, which compile to loops on CPU and kernel launches on GPU.
The constructs differ in the promises they make to the compiler, and picking the wrong one compiles and runs but can produce silently wrong results:

* ``amrex::ParallelFor`` promises the compiler that loop iterations are **independent**: on CPU, the loop is marked with a SIMD pragma (e.g., ``#pragma GCC ivdep``) that hints to the compiler that it may vectorize across iterations.
  Use it only when no two iterations can touch the same memory location, e.g., per-particle pushes or per-cell field updates.
* ``amrex::For`` is identical on GPU but carries no SIMD pragma on CPU.
  Use it whenever different iterations may write the same location (e.g. deposition, histogram bins, etc.).
  Note that ``amrex::Gpu::Atomic`` and ``amrex::HostDevice::Atomic`` operations are plain, non-atomic updates on serial CPU builds and therefore do not make a ``ParallelFor`` loop safe.
* ``amrex::ParallelForSIMD<WIDTH>``: explicitly vectorized CPU loops.
  Rather than relying on compiler auto-vectorization (which ``ParallelFor`` only hints to the compiler), the loop runs in chunks of the compile-time SIMD width and the kernel receives an ``amrex::SIMDindex<WIDTH>`` to perform explicit SIMD loads and stores (with scalar iterations for the remainder).
  ImpactX uses this heavily for the per-particle element pushes: the ``impactx::ParallelFor`` wrapper (``src/elements/mixin/beamoptic.H``) dispatches to ``amrex::ParallelForSIMD<T_Element::simd_width>`` for elements marked as vectorized and falls back to ``amrex::ParallelFor`` otherwise.
  When compiled for GPU, explicit SIMD is inactive (``ImpactX_SIMD`` is a CPU feature) and the dispatch falls back to a regular ``amrex::ParallelFor`` kernel launch, where each GPU thread processes one particle with a scalar index.
* Whole-loop reductions (sums, maxima) should use the ``amrex::ReduceSIMD`` (or ``amrex::Reduce``) functions.
  When compiled for GPU, the ``ReduceSIMD`` code path is inactive and the standard ``amrex::Reduce`` device reduction is used.

Getting ``amrex::ParallelFor`` and atomics wrong is compiler-dependent and silent; see the analysis in `WarpX issue #7097 <https://github.com/BLAST-WarpX/warpx/issues/7097>`__ and the corresponding `WarpX developer documentation on portability <https://warpx.readthedocs.io/en/latest/developers/portability.html>`__ for details.
