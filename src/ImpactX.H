/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell, Ji Qiang
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACT_X_H
#define IMPACT_X_H

#include "elements/All.H"
#include "initialization/AmrCoreData.H"
#include "particles/distribution/All.H"
#include "tracking/TrackingState.H"

#include <AMReX_REAL.H>

#include <functional>
#include <list>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>


namespace impactx
{
    /** An ImpactX simulation
     *
     * This is the central ImpactX simulation class
     */
    class ImpactX
    {
      public:
        /** Construct an ImpactX simulation object
         */
        ImpactX ();

        // removed constructors/assignments
        ImpactX (ImpactX const&) = delete;
        ImpactX (ImpactX &&) = delete;
        void operator= (ImpactX const&) = delete;
        void operator= (ImpactX &&) = delete;

        /** Destruct the ImpactX simulation object
         */
        ~ImpactX ();

        /** Deallocate all contexts and data
         */
        void finalize ();

        /** Finalize elements
         */
        void finalize_elements ();

        /** Initialize AMReX blocks/grids for domain decomposition & space charge mesh.
         *
         * This must come first, before particle beams and lattice elements are
         * initialized.
         */
        void init_grids ();

        /** Initialize the particle beam distribution
         *
         * This clears and initialized the particle beam from input file and
         * command-line options, as parsed by amrex::ParmParse.
         */
        void initBeamDistributionFromInputs ();

        /** Initialize the list of lattice elements
         *
         * This clears and initialized the accelerator lattice elements from
         * input file and command-line options, as parsed by amrex::ParmParse.
         */
        void initLatticeElementsFromInputs ();

        /** Generate and add n particles to the particle container
         *
         * Will also resize the geometry based on the updated particle
         * distribution's extent and then redistribute particles in according
         * AMReX grid boxes.
         *
         * @param bunch_charge bunch charge (C)
         * @param distr distribution function to draw from (object)
         * @param npart number of particles to draw
         * @param spin_distr optional spin distribution
         */
        void
        add_particles (
            amrex::ParticleReal bunch_charge,
            distribution::KnownDistributions distr,
            amrex::Long npart,
            std::optional<distribution::SpinvMF> spin_distr = std::nullopt
        );

        /** Validate the simulation is ready to run the particle tracking loop via @see track_particles
         */
        void validate ();

        /** Check inputs for unused parameters (e.g. typos) after step 1
         *
         * @return true if successful
         */
        bool early_param_check ();

        /** Run the main simulation loop
         */
        void evolve ();

        /** Run the particle tracking simulation loop
         */
        void track_particles ();

        /** Run the envelope tracking simulation loop
         *
         * Modeled as linear transport of the covariance matrix.
         */
        void track_envelope ();

        /** Run the reference orbit tracking simulation loop
         */
        void track_reference (RefPart & ref);

        /** Query input for warning logger variables and set up warning logger accordingly
         *
         * Input variables are: ``always_warn_immediately`` and ``abort_on_warning_threshold``.
         */
        void init_warning_logger ();

        /** The AMReX core object that contains AMR hierarchy, \see amrex::AmrCore
         *
         * Inside here, we store our particle containers and space charge fields.
         * This is held through the precision-agnostic base; the concrete,
         * runtime-selected beam precision lives in initialization::AmrCoreData<T_PT>.
         * Use amr_data_visit() to reach the concrete, precision-templated data.
         */
        std::unique_ptr<initialization::AmrCoreDataBase> amr_data;

        /** Dispatch on the runtime beam precision and invoke a callable.
         *
         * The callable @p f is invoked with a reference to the concrete
         * initialization::AmrCoreData<T_PT> matching the current precision. This
         * is the single boundary where runtime precision crosses into templated
         * code; use a generic lambda, e.g. amr_data_visit([&](auto& data){ ... }).
         */
        template <class F>
        decltype(auto) amr_data_visit (F && f)
        {
            return initialization::visit_amr_data(*amr_data, std::forward<F>(f));
        }

        /** The runtime-selected beam precision.
         *
         * Defaults to the base precision the binary was compiled with; change it
         * with set_precision() before init_grids(). Also overridable via the
         * ``impactx.precision`` input (``single``/``double``).
         */
        initialization::Precision m_precision =
#ifdef IMPACTX_COMPILE_DOUBLE
            initialization::Precision::double_;
#else
            initialization::Precision::single;
#endif

        /** Select the beam floating-point precision at runtime.
         *
         * Must be called before init_grids(); the requested precision must be one
         * compiled into the binary (see the ImpactX_PRECISION build option).
         *
         * @param precision the beam precision to use
         */
        void set_precision (initialization::Precision precision);

        /** Resize the mesh, based on the extent of the bunch of particle
         *
         * This only changes the physical extent of the mesh, but not the
         * number of grid cells.
         */
        void ResizeMesh ();

        /** these are elements defining the accelerator lattice */
        std::list<elements::KnownElements> m_lattice;

        /** Was init_grids already called?
         *
         * Some operations, like resizing a simulation in terms of cells and changing blocking
         * factors are not possible after they were initialized in AMReX. This keeps track
         * of this.
         */
        bool initialized ()
        {
            return m_grids_initialized;
        }

        /** Call user-defined functions while the simulation is running.
         *
         * Check if an optional, user-defined function is registered by
         * name and if yes, call it.
         *
         * @see m_hook for function signature:
         *   void user_defined_hook (ImpactX *);
         */
        void call_hook (std::string const & name);

        std::unordered_map<
            std::string,
            std::function<void(ImpactX *)>
        > m_hook;  //! hooks that users can call

        /** A few states during tracking loops, useful to know in hooks (Python callback functions).
         */
        TrackingState m_tracking_state;

      private:
        /** Keeps track if init_grids was called.
         *
         * Some operations, like resizing a simulation in terms of cells and changing blocking
         * factors are not possible after they were initialized in AMReX. This keeps track
         * of this.
         */
        bool m_grids_initialized = false;
    };

} // namespace impactx

#endif // IMPACT_X_H
