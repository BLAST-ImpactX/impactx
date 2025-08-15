/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Andrew Myers, Marco Garten
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_INITIALIZATION_INITDISTRIBUTION_H
#define IMPACTX_INITIALIZATION_INITDISTRIBUTION_H

#include "particles/ImpactXParticleContainer.H"
#include "particles/CovarianceMatrix.H"
#include "particles/distribution/All.H"

#include <AMReX_Extension.H>  // for AMREX_RESTRICT
#include <AMReX_REAL.H>

#include <optional>  // for std::optional
#include <utility>  // for std::move


namespace impactx::initialization
{
    /** Read the reference particle configuration from inputs
     *
     * @param pp_dist the distribution group in the ParmParser
     * @return a configured reference particle
     */
    RefPart
    read_reference_particle (amrex::ParmParse const & pp_dist);

    /** Read a distribution configuration from inputs
     *
     * @param pp_dist the distribution group in the ParmParser
     * @return a configured distribution
     */
    distribution::KnownDistributions
    read_distribution (amrex::ParmParse const & pp_dist);

    /** Ignore the shape of a distribution and use the 2nd moments to create a covariance matrix
     *
     * @param distr the distribution to draw from
     * @param intensity the beam charge (in C, 3D space charge) or current (in A, 2D space charge)
     */
    Envelope
    create_envelope (
        distribution::KnownDistributions const & distr,
        std::optional<amrex::ParticleReal> intensity = std::nullopt
    );

    /** Initialize a single particle's data using the given distribution
     *
     * Note: we usually would just write a C++ lambda below in ParallelFor. But, due to restrictions
     * in NVCC as of 11.5, we cannot write a lambda in a lambda as we also std::visit the element
     * types of our lattice elements list.
     *    error #3206-D: An extended __device__ lambda cannot be defined inside a generic lambda expression("operator()").
     * Thus, we fall back to writing a C++ functor here, instead of nesting two lambdas.
     *
     * Nvidia bug report: 3458976
     * Minimal demonstrator: https://cuda.godbolt.org/z/39e4q53Ye
     *
     * @tparam T_Distribution This can be a \see Gaussian, \see Waterbag, \see Kurth6D, \see Thermal etc.
     */
    template <typename T_Distribution>
    struct InitSingleParticleData
    {
        /** Constructor taking in pointers to particle data
         *
         * @param distribution the type of distribution function to call
         * @param part_x the array to the particle position (x)
         * @param part_y the array to the particle position (y)
         * @param part_t the array to the particle position (t)
         * @param part_px the array to the particle momentum (x)
         * @param part_py the array to the particle momentum (y)
         * @param part_pt the array to the particle momentum (t)
         */
        InitSingleParticleData (
            T_Distribution distribution,
            amrex::ParticleReal* AMREX_RESTRICT part_x,
            amrex::ParticleReal* AMREX_RESTRICT part_y,
            amrex::ParticleReal* AMREX_RESTRICT part_t,
            amrex::ParticleReal* AMREX_RESTRICT part_px,
            amrex::ParticleReal* AMREX_RESTRICT part_py,
            amrex::ParticleReal* AMREX_RESTRICT part_pt
        )
        : m_distribution(std::move(distribution)),
          m_part_x(part_x), m_part_y(part_y), m_part_t(part_t),
          m_part_px(part_px), m_part_py(part_py), m_part_pt(part_pt)
        {
        }

        InitSingleParticleData () = delete;
        InitSingleParticleData (InitSingleParticleData const &) = default;
        InitSingleParticleData (InitSingleParticleData &&) = default;
        ~InitSingleParticleData () = default;

        /** Initialize the data for a single particle
         *
         * @param i particle index in the current box
         * @param engine a random number engine (with associated state)
         */
        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        void
        operator() (
            amrex::Long i,
            amrex::RandomEngine const & engine
        ) const
        {
            m_distribution(
                m_part_x[i],
                m_part_y[i],
                m_part_t[i],
                m_part_px[i],
                m_part_py[i],
                m_part_pt[i],
                engine
            );
        }

    private:
        T_Distribution const m_distribution;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_x;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_y;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_t;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_px;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_py;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_pt;
    };


    /** Initialize the input parameters for all distributions that read phase space ellipse parameters from the input
     *
     * This function sets the distribution parameters based on the provided phase space inputs.
     * The parameters include the phase space ellipse intersections for position (`lambda_x`, `lambda_y`), time (`lambda_t`),
     * momentum (`lambda_px`, `lambda_py`) and energy (`lambda_pt`), as well as the correlation terms (`mu_x_px`, `mu_y_py`, `mu_t_pt`).
     *
     * @param pp_dist The parameter parser object.
     * @param lambda_x Reference to the variable storing the axis intersection of the phase space ellipse for x position.
     * @param lambda_y Reference to the variable storing the axis intersection of the phase space ellipse for y position.
     * @param lambda_t Reference to the variable storing the axis intersection of the phase space ellipse for time.
     * @param lambda_px Reference to the variable storing the axis intersection of the phase space ellipse for momentum in x direction.
     * @param lambda_py Reference to the variable storing the axis intersection of the phase space ellipse for momentum in y direction.
     * @param lambda_pt Reference to the variable storing the axis intersection of the phase space ellipse for energy.
     * @param mu_x_px Reference to the variable storing the correlation coefficient between the x position and the x momentum.
     * @param mu_y_py Reference to the variable storing the correlation coefficient between the y position and the y momentum.
     * @param mu_t_pt Reference to the variable storing the correlation coefficient between the t position and the t momentum.
     * @param mean_x Reference to the variable storing the mean for position in the x direction.
     * @param mean_y Reference to the variable storing the mean for position in the y direction.
     * @param mean_t Reference to the variable storing the mean for position in the t direction.
     * @param mean_px Reference to the variable storing the mean for momentum in the x direction.
     * @param mean_py Reference to the variable storing the mean for momentum in the y direction.
     * @param mean_pt Reference to the variable storing the mean for momentum in the t direction.
     * @param dispersion_x Reference to the variable storing the beam-based dispersion in the x direction.
     * @param dispersion_px Reference to the variable storing the derivative of the beam-based dispersion in the x direction.
     * @param dispersion_y Reference to the variable storing the beam-based dispersion in the y direction.
     * @param dispersion_py Reference to the variable storing the derivative of the beam-based dispersion in the y direction.
     */
    void
    set_distribution_parameters_from_phase_space_inputs (
        amrex::ParmParse const & pp_dist,
        amrex::ParticleReal& lambda_x, amrex::ParticleReal& lambda_y, amrex::ParticleReal& lambda_t,
        amrex::ParticleReal& lambda_px,amrex::ParticleReal& lambda_py, amrex::ParticleReal& lambda_pt,
        amrex::ParticleReal& mu_x_px, amrex::ParticleReal& mu_y_py, amrex::ParticleReal& mu_t_pt,
        amrex::ParticleReal& mean_x, amrex::ParticleReal& mean_y, amrex::ParticleReal& mean_t,
        amrex::ParticleReal& mean_px, amrex::ParticleReal& mean_py, amrex::ParticleReal& mean_pt,
        amrex::ParticleReal& dispersion_x, amrex::ParticleReal& dispersion_px, amrex::ParticleReal& dispersion_y, amrex::ParticleReal& dispersion_py
    );


    /** Set the distribution parameters from Twiss inputs
     *
     * This function reads Courant-Snyder / Twiss parameters from the provided `ParmParse` object and calculates the distribution parameters for a particle beam.
     * It sets the values for the axis intercepts of the phase space ellipse (`lambda_x`, `lambda_y`, `lambda_t`, and `lambda_px`, `lambda_py`, `lambda_pt`) and correlation terms (`mu_x_px`, `mu_y_py`, `mu_t_pt`)  normally accepted as input.
     *
     * @param pp_dist The parameter parser object.
     * @param lambda_x Reference to the variable storing the axis intersection of the phase space ellipse for x position.
     * @param lambda_y Reference to the variable storing the axis intersection of the phase space ellipse for y position.
     * @param lambda_t Reference to the variable storing the axis intersection of the phase space ellipse for time.
     * @param lambda_px Reference to the variable storing the axis intersection of the phase space ellipse for momentum in x direction.
     * @param lambda_py Reference to the variable storing the axis intersection of the phase space ellipse for momentum in y direction.
     * @param lambda_pt Reference to the variable storing the axis intersection of the phase space ellipse for energy.
     * @param mu_x_px Reference to the variable storing the correlation coefficient between the x position and the x momentum.
     * @param mu_y_py Reference to the variable storing the correlation coefficient between the y position and the y momentum.
     * @param mu_t_pt Reference to the variable storing the correlation coefficient between the t position and the t momentum.
     * @param mean_x Reference to the variable storing the mean for position in the x direction.
     * @param mean_y Reference to the variable storing the mean for position in the y direction.
     * @param mean_t Reference to the variable storing the mean for position in the t direction.
     * @param mean_px Reference to the variable storing the mean for momentum in the x direction.
     * @param mean_py Reference to the variable storing the mean for momentum in the y direction.
     * @param mean_pt Reference to the variable storing the mean for momentum in the t direction.
     * @param dispersion_x Reference to the variable storing the beam-based dispersion in the x direction.
     * @param dispersion_px Reference to the variable storing the derivative of the beam-based dispersion in the x direction.
     * @param dispersion_y Reference to the variable storing the beam-based dispersion in the y direction.
     * @param dispersion_py Reference to the variable storing the derivative of the beam-based dispersion in the y direction.
     */
    void
    set_distribution_parameters_from_twiss_inputs (
        amrex::ParmParse const & pp_dist,
        amrex::ParticleReal& lambda_x, amrex::ParticleReal& lambda_y, amrex::ParticleReal& lambda_t,
        amrex::ParticleReal& lambda_px,amrex::ParticleReal& lambda_py, amrex::ParticleReal& lambda_pt,
        amrex::ParticleReal& mu_x_px, amrex::ParticleReal& mu_y_py, amrex::ParticleReal& mu_t_pt,
        amrex::ParticleReal& mean_x, amrex::ParticleReal& mean_y, amrex::ParticleReal& mean_t,
        amrex::ParticleReal& mean_px, amrex::ParticleReal& mean_py, amrex::ParticleReal& mean_pt,
        amrex::ParticleReal& dispersion_x, amrex::ParticleReal& dispersion_px, amrex::ParticleReal& dispersion_y, amrex::ParticleReal& dispersion_py
    );

}  // namespace impactx::initialization

#endif  // IMPACTX_INITIALIZATION_INITDISTRIBUTION_H
