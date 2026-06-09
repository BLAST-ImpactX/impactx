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

    template <class T_PC>
    amrex::Gpu::DeviceVector<amrex::Real>
    Deposit1D (
        T_PC & pc,
        amrex::Real bin_min,
        amrex::Real bin_max,
        int num_bins
    )
    {
        BL_PROFILE("impactx::spacecharge::Deposit1D");

        using namespace amrex::literals;

        // Set parameters for charge deposition
        bool const is_unity_particle_weight = false;

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

        return charge_distribution;
    }

    // explicit instantiations for the compiled beam precisions
#ifdef IMPACTX_COMPILE_DOUBLE
    template amrex::Gpu::DeviceVector<amrex::Real>
    Deposit1D (
        ImpactXParticleContainerT<double> & pc,
        amrex::Real bin_min,
        amrex::Real bin_max,
        int num_bins
    );
#endif
#ifdef IMPACTX_COMPILE_SINGLE
    template amrex::Gpu::DeviceVector<amrex::Real>
    Deposit1D (
        ImpactXParticleContainerT<float> & pc,
        amrex::Real bin_min,
        amrex::Real bin_max,
        int num_bins
    );
#endif

}  // namespace impactx::particles::spacecharge
