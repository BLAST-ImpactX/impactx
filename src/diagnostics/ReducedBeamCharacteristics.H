/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
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
     *
     * @param pc the particle beam
     * @param covariance optional: if not nullptr, receives the 6x6 covariance matrix of the
     *                   phase space coordinates x, px, y, py, t, pt
     * @throws std::runtime_error if covariance is not nullptr and the beam has no particle
     *         weight or the reference particle is not initialized
     */
    std::unordered_map<std::string, amrex::ParticleReal>
    reduced_beam_characteristics (
        ImpactXParticleContainer const & pc,
        CovarianceMatrix * covariance = nullptr
    );

    /** Compute the covariance matrix of the beam distribution from particle data
     *
     * The 6x6 covariance matrix of the phase space coordinates x, px, y, py, t, pt, as used
     * in envelope tracking. This uses an MPI Allreduce and returns a result on all ranks.
     *
     * @param pc the particle beam
     * @return the covariance matrix
     * @throws std::runtime_error if the beam has no particle weight or the reference particle
     *         is not initialized
     */
    CovarianceMatrix
    covariance_matrix (ImpactXParticleContainer const & pc);

    /** Compute moments of the beam distribution from covariance matrix
     */
    std::unordered_map<std::string, amrex::ParticleReal>
    reduced_beam_characteristics (Map6x6 const & cm, RefPart const & ref_part);

} // namespace impactx::diagnostics

#endif // IMPACTX_REDUCED_BEAM_CHARACTERISTICS
