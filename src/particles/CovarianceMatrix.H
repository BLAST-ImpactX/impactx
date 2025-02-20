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
        CovarianceMatrix m_env; ///< the 6x6 beam covariance matrix
        amrex::ParticleReal m_beam_intensity = 0.0; ///< optional: charge in A (for 3D space charge) or current in A (for 2D space charge)

        /** Set envelope beam charge/current for 3D/2D space charge
         *
         * @param intensity beam charge (C) in 3D or beam current (A) in 2D
         */
        Envelope &
        set_beam_intensity (amrex::ParticleReal const intensity)
        {
            m_beam_intensity = intensity;

            return *this;
         }

        /** Get envelope beam charge/current for 3D/2D space charge
         *
         * @returns 3D: beam charge in C; 2D: beam current in A
         */
        amrex::ParticleReal
        beam_intensity () const
        {
            using namespace amrex::literals;

            return m_beam_intensity;
        }

        /** Set 6x6 covariance matrix for envelope tracking
         *
         * @param covariance_matrix beam 6x6 covariance matrix
         */
        Envelope &
        set_covariance_matrix (CovarianceMatrix const & covariance_matrix)
        {
            m_env = covariance_matrix;

            return *this;
         }

    };


} // namespace impactx::distribution

#endif // IMPACTX_DISTRIBUTION_COVARIANCE_MATRIX_H
