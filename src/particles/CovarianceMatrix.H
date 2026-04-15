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

    /** Return the 6x6 skew-symmetric symplectic form
     *
     * diag([[0,1],[-1,0]], [[0,1],[-1,0]], [[0,1],[-1,0]])
     * for the phase space ordering (x, px, y, py, t, pt).
     *
     * Note: commonly referred to as J in accelerator physics, but
     * properly called Omega to avoid confusion with complex structures.
     *
     * @see https://en.wikipedia.org/wiki/Symplectic_matrix
     * Refs.:
     * - K. R. Meyer and D. C. Offin, Introduction to Hamiltonian Dynamical Systems and the N-Body Problem, Springer-Verlag NY, 2009.
     * - A. J. Dragt, Lie Methods for Nonlinear Dynamics with Applications to Accelerator Physics, University of Maryland, 2020, http://www.physics.umd.edu/dsat/
     */
    inline Map6x6 symplectic_form ()
    {
        Map6x6 Omega = Map6x6::Zero();
        for (int i = 0; i < 3; ++i) {
            Omega(2*i+1, 2*i+2) =  1;
            Omega(2*i+2, 2*i+1) = -1;
        }
        return Omega;
    }


    /** This struct stores the beam envelope attributes, including the 6x6 beam
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

        /** Get the 6x6 covariance matrix for envelope tracking
         */
        CovarianceMatrix
        covariance_matrix () const
        {
            return m_env;
        }

    };


} // namespace impactx::distribution

#endif // IMPACTX_DISTRIBUTION_COVARIANCE_MATRIX_H
