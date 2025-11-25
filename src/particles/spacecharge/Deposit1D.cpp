/* Copyright 2022-2025 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Ji Qiang
 * License: BSD-3-Clause-LBNL
 */
#include "Deposit1D.H"

#include "diagnostics/ReducedBeamCharacteristics.H"
#include "particles/wakefields/ChargeBinning.H"

#include <AMReX_REAL.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_ParmParse.H>

#include <cmath>


namespace impactx::particles::spacecharge
{

    void Deposit1D (
        ImpactXParticleContainer & pc,
        amrex::Real * beam_profile,
        amrex::Real * beam_profile_slope,
        amrex::Real bin_min,
        amrex::Real bin_max,
        int num_bins
    )
    {
        BL_PROFILE("impactx::spacecharge::Deposit1D");

        using namespace amrex::literals;

        // Set parameters for charge deposition
        bool const is_unity_particle_weight = false;
        bool const GetNumberDensity = true;

        amrex::Real const bin_size = (bin_max - bin_min) / (num_bins - 1);  // number of evaluation points
        // Allocate memory for the charge profile
        amrex::Gpu::DeviceVector<amrex::Real> charge_distribution(num_bins + 1, 0.0);
        // Call charge deposition function
        impactx::particles::wakefields::DepositCharge1D(pc, charge_distribution, bin_min, bin_size, is_unity_particle_weight);

        // Sum up all partial charge histograms to each MPI process to calculate
        // the global charge slope.
        amrex::ParallelAllReduce::Sum(
            charge_distribution.data(),
            charge_distribution.size(),
            amrex::ParallelDescriptor::Communicator()
        );

        // Call charge density derivative function
        amrex::Gpu::DeviceVector<amrex::Real> slopes(charge_distribution.size() - 1, 0.0);
        impactx::particles::wakefields::DerivativeCharge1D(charge_distribution, slopes, bin_size,GetNumberDensity);

        beam_profile = charge_distribution.data();
        beam_profile_slope = slopes.data();

    }

}  // namespace impactx::particles::spacecharge
