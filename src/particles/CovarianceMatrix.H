/* Copyright 2022-2024 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DISTRIBUTION_COVARIANCE_MATRIX_H
#define IMPACTX_DISTRIBUTION_COVARIANCE_MATRIX_H

#include <AMReX_REAL.H>
#include <AMReX_SmallMatrix.H>


namespace impactx
{
    /** this is a 6x6 matrix */
    using Map6x6 = amrex::SmallMatrix<amrex::ParticleReal, 6, 6, amrex::Order::F, 1>;

    /** this is the 6x6 covariance matrix */
    using CovarianceMatrix = Map6x6;


    /** This struct stores the beam envelope attribues, including the 6x6 beam
     *  covariance matrix.  Used during envelope tracking mode.
     */
    struct Envelope
    {
        amrex::ParticleReal bunch_charge = 0.0; ///< total charge in C (for 3D space charge)
        amrex::ParticleReal beam_current = 0.0; ///< current in A (for 2D space charge)
        CovarianceMatrix cm; ///< the 6x6 beam covariance matrix

        /** Set envelope beam current for 2D space charge
         *
         * @param beam_current_A beam current (A)
         */
        Envelope &
        set_beam_current_A (amrex::ParticleReal const beam_current_A)
        {
            using namespace amrex::literals;

            beam_current = beam_current_A;

            return *this;
         }

        /** Get envelope beam current for 2D space charge
         *
         * @returns beam current in A
         */
        amrex::ParticleReal
        beam_current_A ()
        {
            using namespace amrex::literals;

            return beam_current;
        }

        /** Set 6x6 covariance matrix for envelope tracking
         *
         * @param covariance_matrix beam 6x6 covariance matrix
         */
        Envelope &
        set_covariance_matrix (const CovarianceMatrix& covariance_matrix)
        {
            cm = covariance_matrix;

            return *this;
         }

    };


} // namespace impactx::distribution

#endif // IMPACTX_DISTRIBUTION_COVARIANCE_MATRIX_H
