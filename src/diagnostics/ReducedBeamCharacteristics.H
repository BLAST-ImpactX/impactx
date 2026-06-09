/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Marco Garten, Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_REDUCED_BEAM_CHARACTERISTICS_H
#define IMPACTX_REDUCED_BEAM_CHARACTERISTICS_H

#include "particles/ImpactXParticleContainer.H"
#include "particles/CovarianceMatrix.H"

#include <AMReX_REAL.H>

#include <string>
#include <unordered_map>


namespace impactx::diagnostics
{
    /** Compute moments of the beam distribution from particle data
     *
     * This uses an MPI Allreduce and returns a result on all ranks.
     * Reductions are accumulated in amrex::ParticleReal (double) regardless of
     * the beam storage precision.
     *
     * @tparam T_PC the beam particle container type (any compiled precision)
     */
    template <class T_PC>
    std::unordered_map<std::string, amrex::ParticleReal>
    reduced_beam_characteristics (T_PC const & pc);

    /** Compute moments of the beam distribution from covariance matrix
     */
    std::unordered_map<std::string, amrex::ParticleReal>
    reduced_beam_characteristics (Map6x6 const & cm, RefPart const & ref_part);

} // namespace impactx::diagnostics

#endif // IMPACTX_REDUCED_BEAM_CHARACTERISTICS
