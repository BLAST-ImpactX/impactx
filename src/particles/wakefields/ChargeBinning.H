/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Alex Bojanich, Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef CHARGE_BINNING_H
#define CHARGE_BINNING_H

#include "particles/ImpactXParticleContainer.H"


namespace impactx::particles::wakefields
{
    /** Function to calculate charge density profile
     *
     * @param[in] myspc the particle species to deposit along s
     * @param[out] charge_distribution the output array (1D)
     * @param[in] bin_min lower end of the beam in s
     * @param[in] bin_size size of the beam in s divided by num_bins
     * @param[in] is_unity_particle_weight ignore the particle weighting per macro particle,
     *                                     otherwise treat each particle as one physical particle
     */
    void DepositCharge1D (
        impactx::ImpactXParticleContainer& myspc,
        amrex::Gpu::DeviceVector<amrex::Real> & charge_distribution,
        amrex::Real bin_min,
        amrex::Real bin_size,
        bool is_unity_particle_weight = false
    );

    /** Function to calculate the slope of the number (or charge) density
     *
     * @param[in] charge_distribution deposited charge in 1D along s
     * @param[out] slopes derivative of charge_distribution along s with len(charge_distribution)-1
     * @param[in] bin_size size of the beam in s divided by num_bins
     * @param[in] GetNumberDensity number density if true, otherwise charge density
     */
    void DerivativeCharge1D (
        amrex::Gpu::DeviceVector<amrex::Real> const & charge_distribution,
        amrex::Gpu::DeviceVector<amrex::Real> & slopes,
        amrex::Real bin_size,
        bool GetNumberDensity = true
    );

    /** Function to calculate the mean transverse positions (x and y) at each bin along s
     *
     * @param[in] myspc the particle species to deposit along s
     * @param[out] mean_x the output array for mean x positions (1D)
     * @param[out] mean_y the output array for mean y positions (1D)
     * @param[in] bin_min lower end of the beam in s
     * @param[in] bin_size size of the beam in s divided by num_bins
     * @param[in] is_unity_particle_weight ignore the particle weighting per macro particle,
     *                                     otherwise treat each particle as one physical particle
     */
    void MeanTransversePosition (
        impactx::ImpactXParticleContainer& myspc,
        amrex::Gpu::DeviceVector<amrex::Real> & mean_x,
        amrex::Gpu::DeviceVector<amrex::Real> & mean_y,
        amrex::Real bin_min,
        amrex::Real bin_size,
        bool is_unity_particle_weight = false
    );

    /** All-reduce (sum) a 1D histogram across MPI ranks, in place.
     *
     * The histogram is stored in (non-managed) device memory. A device pointer
     * is only handed to MPI when AMReX is configured for GPU-aware MPI;
     * otherwise the data is staged on the host for the collective. No-op for a
     * single rank.
     *
     * @param[in,out] data device-resident 1D array, summed across all ranks
     */
    void AllReduceSum1D (
        amrex::Gpu::DeviceVector<amrex::Real> & data
    );

    /** Reduce (sum) a 1D histogram onto a single MPI rank, in place.
     *
     * Same GPU-aware-MPI staging behavior as @ref AllReduceSum1D, but the
     * reduced result is only valid on rank @p root.
     *
     * @param[in,out] data device-resident 1D array, summed onto @p root
     * @param[in] root the MPI rank that receives the reduced result
     */
    void ReduceSum1D (
        amrex::Gpu::DeviceVector<amrex::Real> & data,
        int root
    );

    /** Broadcast a 1D histogram from the IO processor to every MPI rank, in place.
     *
     * On entry @p data is only valid on the IO processor; on exit every rank holds
     * a correctly-sized copy (receivers are resized to match). Same GPU-aware-MPI
     * staging behavior as @ref AllReduceSum1D. No-op for a single rank.
     *
     * @param[in,out] data device-resident 1D array, broadcast from the IO processor
     */
    void BcastFromIOProc1D (
        amrex::Gpu::DeviceVector<amrex::Real> & data
    );

} // namespace impactx::particles::wakefields

#endif // CHARGE_BINNING_H
