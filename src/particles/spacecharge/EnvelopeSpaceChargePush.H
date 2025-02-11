/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ENVELOPESPACECHARGEPUSH_H
#define IMPACTX_ENVELOPESPACECHARGEPUSH_H

#include "particles/ImpactXParticleContainer.H"
#include "particles/CovarianceMatrix.H"

#include <AMReX_REAL.H>


namespace impactx::spacecharge
{
     /** This function pushes the 6x6 beam covariance matrix for a slice
      *  of length ds, using the linear space charge fields in an rms
      *  equivalent 2D ellipse, as determined from the beam covariance matrix.
      *  Note: This is a reduced model of 2D space charge.
      *
      * @param[in] refpart reference particle
      * @param[inout] cm covariance matrix
      * @param[in] current beam current [A]
      * @param[in] ds step size [m]
      */
     void
     envelope_space_charge2D_push (
         RefPart const & refpart,
         Map6x6 & cm,
         amrex::ParticleReal & current,
         amrex::ParticleReal & ds
    );

} // namespace impactx

#endif // IMPACTX_ENVELOPESPACECHARGEPUSH_H
