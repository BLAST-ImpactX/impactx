/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_PARTICLES_WAKEFIELDS_CSRKICKMODEL_H_
#define IMPACTX_PARTICLES_WAKEFIELDS_CSRKICKMODEL_H_

#include "particles/ReferenceParticle.H"

#include <AMReX_GpuContainers.H>
#include <AMReX_REAL.H>

#include <functional>
#include <string>


namespace impactx::particles::wakefields
{
    /** Binned longitudinal beam profile passed to a user-provided CSR kick model
     *
     * All arrays share the same binning: bin i spans
     * [bin_min + i * bin_size, bin_min + (i+1) * bin_size) in the longitudinal
     * coordinate t = ct [m]. The arrays have num_bins + 1 entries; the last
     * entry is a guard bin that is usually (close to) zero.
     *
     * The arrays live in device memory (on GPU builds) and are only valid
     * during the model call. Retained data must be copied.
     */
    struct CSRProfile
    {
        /** line charge density lambda(t) [C/m] from NGP deposition;
         * sum(charge) * bin_size approximates the bunch charge [C] */
        amrex::Gpu::DeviceVector<amrex::Real> charge;

        /** charge-weighted mean of x per bin [m]; zero for empty bins */
        amrex::Gpu::DeviceVector<amrex::Real> mean_x;

        /** charge-weighted mean of y per bin [m]; zero for empty bins */
        amrex::Gpu::DeviceVector<amrex::Real> mean_y;

        amrex::Real bin_min = 0.0;   //!< lower edge of bin 0 in t = ct [m]
        amrex::Real bin_size = 0.0;  //!< bin spacing [m]
        int num_bins = 0;            //!< number of kick bins (= algo.csr_bins)
    };

    /** Per-element, per-slice context passed to a user-provided CSR kick model
     */
    struct CSRElementContext
    {
        std::string element_name;  //!< user-given element name (empty if unnamed)
        std::string element_type;  //!< element type, e.g. "Sbend"

        amrex::Real rc = 0.0;          //!< radius-of-curvature magnitude |R| [m]
        amrex::Real signed_rc = 0.0;   //!< signed radius of curvature R [m]
        amrex::ParticleReal ds = 0.0;  //!< element (arc) length [m]

        int nslice = 1;  //!< number of slices in the element
        int slice = 0;   //!< current slice index, 0-based, in forward orientation

        /** distance from the element entrance at kick time [m]
         * (= slice * slice_ds, also during back-tracking) */
        amrex::ParticleReal s = 0.0;

        amrex::ParticleReal slice_ds = 0.0;  //!< slice length [m]

        RefPart ref;  //!< reference particle (snapshot copy)
    };

    /** Per-bin kick returned by a user-provided CSR kick model
     *
     * All entries are forces in Newtons on bins 0 .. num_bins-1.
     * pt is required (size num_bins); px and py are optional (empty = not
     * applied). \see WakePush for how the kicks are applied.
     */
    struct CSRKick
    {
        amrex::Gpu::DeviceVector<amrex::Real> pt;  //!< longitudinal force F_t [N]
        amrex::Gpu::DeviceVector<amrex::Real> px;  //!< horizontal force F_x [N] (optional)
        amrex::Gpu::DeviceVector<amrex::Real> py;  //!< vertical force F_y [N] (optional)
    };

    /** User-provided CSR kick model
     *
     * Maps the binned longitudinal beam profile and the element context to
     * per-bin kicks, replacing the built-in analytic CSR wake model
     * (e.g., with a pre-trained ML surrogate).
     */
    using CSRKickModel = std::function<CSRKick (CSRProfile const &, CSRElementContext const &)>;

} // namespace impactx::particles::wakefields

#endif // IMPACTX_PARTICLES_WAKEFIELDS_CSRKICKMODEL_H_
