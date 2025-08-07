/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell, Ji Qiang, Remi Lehe
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactX.H"
#include "initialization/InitAmrCore.H"
#include "particles/ImpactXParticleContainer.H"
#include "particles/Push.H"
#include "elements/All.H"
#include "envelope/spacecharge/EnvelopeSpaceChargePush.H"
#include "diagnostics/DiagnosticOutput.H"
#include "diagnostics/ReducedBeamCharacteristics.H"
#include "particles/wakefields/HandleWakefield.H"

#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Math.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#include <array>
#include <iostream>
#include <memory>


namespace impactx {
    ImpactX::ImpactX() {
        // todo: if amr.n_cells is provided, overwrite/redefine AmrCore object

        // todo: if charge deposition and/or space charge are requested, require
        //       amr.n_cells from user inputs
    }

    ImpactX::~ImpactX()
    {
        this->finalize();
    }

    void ImpactX::finalize ()
    {
        // loop over all beamline elements & finalize them
        finalize_elements();

        if (m_grids_initialized)
        {
            m_lattice.clear();

            // this one last
            amr_data.reset();

            if (amrex::Initialized())
                amrex::Finalize();

            // only finalize once
            m_grids_initialized = false;
        }
    }

    void ImpactX::finalize_elements ()
    {
        // loop over all beamline elements & finalize them
        for (auto & element_variant : m_lattice)
        {
            std::visit([](auto&& element){
                element.finalize();
            }, element_variant);
        }
    }

    void ImpactX::init_grids ()
    {
        BL_PROFILE("ImpactX::init_grids");

        amr_data = std::make_unique<initialization::AmrCoreData>(initialization::init_amr_core());
        amr_data->track_particles.m_particle_container = std::make_unique<ImpactXParticleContainer>(amr_data.get());
        amr_data->track_particles.m_particles_lost = std::make_unique<ImpactXParticleContainer>(amr_data.get());

        // query input for warning logger variables and set up warning logger accordingly
        init_warning_logger();

        // move old diagnostics out of the way
        bool diag_enable = true;
        amrex::ParmParse("diag").queryAdd("enable", diag_enable);
        if (diag_enable) {
            amrex::UtilCreateCleanDirectory("diags", true);
        }

        // the particle container has been set to track the same Geometry as ImpactX

        // this is the earliest point that we need to know the particle shape,
        // so that we can initialize the guard size of our MultiFabs
        amr_data->track_particles.m_particle_container->SetParticleShape();

        // init blocks / grids & MultiFabs
        amr_data->InitFromScratch(0.0);

        // prepare particle containers
        //   have to do this here, not in the constructor because grids have not
        //   been built when constructor was called.
        amr_data->track_particles.m_particle_container->prepare();
        amr_data->track_particles.m_particles_lost->prepare();

        // register shortcut
        amr_data->track_particles.m_particle_container->SetLostParticleContainer(amr_data->track_particles.m_particles_lost.get());

        // print AMReX grid summary
        if (amrex::ParallelDescriptor::IOProcessor()) {
            // verbosity
            amrex::ParmParse pp_impactx("impactx");
            int verbose = 1;
            pp_impactx.queryAddWithParser("verbose", verbose);

            if (verbose > 0) {
                std::cout << "\nGrids Summary:\n";
                amr_data->printGridSummary(std::cout, 0, amr_data->finestLevel());
            }
        }

        // keep track that init is done
        m_grids_initialized = true;
    }

    void ImpactX::evolve ()
    {
        BL_PROFILE("ImpactX::evolve");

        amrex::ParmParse pp_algo("algo");
        std::string track = "particles";
        pp_algo.queryAdd("track", track);

        if (track == "particles") {
            track_particles();
        }
        else if (track == "envelope") {
            track_envelope();
        }
        else if (track == "reference_orbit") {
            if (!amr_data->track_reference.m_ref.has_value())
            {
                throw std::runtime_error("evolve: Reference particle not set.");
            }
            track_reference(amr_data->track_reference.m_ref.value());
        }
        else {
            throw std::runtime_error("Unknown tracking algorithm: algo.track=" + track);
        }
    }

    void
    evaluate_lattice (ImpactX * sim, amrex::ParticleReal * error, amrex::ParticleReal const * q1_k, amrex::ParticleReal const * q2_k) {
        // Example from:
        // https://impactx.readthedocs.io/en/latest/usage/examples/fodo_space_charge/README.html

        // reference
        // RefPart & ref = sim->amr_data->track_envelope.m_ref.value();
        // auto & env = sim->amr_data->track_envelope.m_env.value();
        // copy
        RefPart ref = sim->amr_data->track_envelope.m_ref.value();
        Envelope env = sim->amr_data->track_envelope.m_env.value();

        int nslice = 50;

        auto dr1 = elements::Drift{7.44e-2, 0, 0, 0, 0, 0, nslice};
        auto q1 = elements::Quad{6.10e-2, *q1_k, 0, 0, 0, 0, 0, nslice};
        auto dr2 = elements::Drift{14.88e-2, 0, 0, 0, 0, 0, nslice};
        auto q2 = elements::Quad{6.10e-2, *q2_k, 0, 0, 0, 0, 0, nslice};

        /*
        // quadrupole triplet
        // https://impactx.readthedocs.io/en/latest/usage/examples/optimize_triplet/README.html
        sim->m_lattice = {
            elements::Drift{2.7},
            elements::Quad{0.1, q1_k},
            elements::Drift{1.4},
            elements::Quad{0.2, q2_k},
            elements::Drift{1.4},
            elements::Quad{0.1, q1_k},
            elements::Drift{2.7}
        };


        // int ns = 25;  // number of slices per ds in the element
        // for (auto & e : sim.m_lattice) { e.m_nslice = ns; }
        */

        //sim->evolve();
        //auto & pc = *sim->amr_data->m_particle_container;
        //dr1(pc, 0, 0);

        // envelope mode
        //std::cout << "before ref\n";
        //std::cout << "before cm\n";
        auto cm = env.m_env;

        auto & intensity = env.m_beam_intensity;
        intensity = 0.5;  // Ampere

        // initial values alpha and beta
        // goal: periodic solution where the resulting alpha and beta are the same as the initial ones.
        //       the resulting -q1_k should be equal to q2_k (but we intentionally do not enforce that
        //       during the simulation loop here, and leave that to the optimizer to find out).
        auto rbc_goal = diagnostics::reduced_beam_characteristics(cm, ref);
        amrex::ParticleReal const alpha_x_goal = rbc_goal["alpha_x"];
        amrex::ParticleReal const alpha_y_goal = rbc_goal["alpha_y"];
        amrex::ParticleReal const beta_x_goal = rbc_goal["beta_x"];
        amrex::ParticleReal const beta_y_goal = rbc_goal["beta_y"];
        std::array<amrex::ParticleReal, 4> alpha_beta_goal = {alpha_x_goal, alpha_y_goal, beta_x_goal, beta_y_goal};

        // sub-steps for space charge within the element
        for (int slice_step = 0; slice_step < nslice; ++slice_step)
        {
            envelope::spacecharge::space_charge2D_push(ref, cm, intensity, dr1.ds() / nslice);
            dr1(ref);
            dr1(cm, ref);
        }

        for (int slice_step = 0; slice_step < nslice; ++slice_step)
        {
            envelope::spacecharge::space_charge2D_push(ref, cm, intensity, q1.ds() / nslice);
            q1(ref);
            q1(cm, ref);
        }

        for (int slice_step = 0; slice_step < nslice; ++slice_step)
        {
            envelope::spacecharge::space_charge2D_push(ref, cm, intensity, dr2.ds() / nslice);
            dr2(ref);
            dr2(cm, ref);
        }

        for (int slice_step = 0; slice_step < nslice; ++slice_step)
        {
            envelope::spacecharge::space_charge2D_push(ref, cm, intensity, q2.ds() / nslice);
            q2(ref);
            q2(cm, ref);
        }

        for (int slice_step = 0; slice_step < nslice; ++slice_step)
        {
            envelope::spacecharge::space_charge2D_push(ref, cm, intensity, dr1.ds() / nslice);
            dr1(ref);
            dr1(cm, ref);
        }

        /*
        // particles
        std::unordered_map<std::string, amrex::ParticleReal> const rbc =
                diagnostics::reduced_beam_characteristics(*active_sim->amr_data->m_particle_container);

        return rbc.at("alpha_x");  // TOOD: alpha_x, alpha_y, beta_x, beta_y
        */

        // envelope
        /*
        std::unordered_map<std::string, amrex::ParticleReal> const rbc =
            diagnostics::reduced_beam_characteristics(cm, ref);
        return rbc.at("alpha_x");
        */

        auto rbc_is = diagnostics::reduced_beam_characteristics(cm, ref);
        amrex::ParticleReal const alpha_x = rbc_is["alpha_x"];
        amrex::ParticleReal const alpha_y = rbc_is["alpha_y"];
        amrex::ParticleReal const beta_x = rbc_is["beta_x"];
        amrex::ParticleReal const beta_y = rbc_is["beta_y"];
        std::array<amrex::ParticleReal, 4> alpha_beta_is = {alpha_x, alpha_y, beta_x, beta_y};

        *error =
            amrex::Math::powi<2>(alpha_beta_is[0] - alpha_beta_goal[0]) +
            amrex::Math::powi<2>(alpha_beta_is[1] - alpha_beta_goal[1]) +
            amrex::Math::powi<2>(alpha_beta_is[2] - alpha_beta_goal[2]) +
            amrex::Math::powi<2>(alpha_beta_is[3] - alpha_beta_goal[3]);
    }

    void my_run ()
    {
        ImpactX sim;

        amrex::ParmParse pp_algo("algo");
        std::string track = "envelope";
        pp_algo.add("track", track);

        sim.init_grids();

        // TODO: replace with beam params from https://impactx.readthedocs.io/en/latest/usage/examples/fodo_space_charge/README.html
        sim.initBeamDistributionFromInputs();

        // design the accelerator lattice
        // sim.initLatticeElementsFromInputs();

        // initial quad strengths
        //amrex::ParticleReal q1_k = -103.12574100336;
        //amrex::ParticleReal q2_k = -q1_k;
        amrex::ParticleReal q1_k = -85;
        amrex::ParticleReal q2_k = 80;
        amrex::Print() << "q1_k = " << q1_k << "\n"
                       << "q2_k = " << q2_k
                       << std::endl;

        amrex::ParticleReal derror_dq1_k = std::numeric_limits<amrex::ParticleReal>::quiet_NaN();
        amrex::ParticleReal derror_dq2_k = std::numeric_limits<amrex::ParticleReal>::quiet_NaN();

#define MYMODE 2

#if MYMODE == 0
        // non-differentiable run:
        amrex::Print() << "Gradient-free run\n";
        amrex::ParticleReal dq1_k = std::numeric_limits<amrex::ParticleReal>::quiet_NaN();
        amrex::ParticleReal dq2_k = std::numeric_limits<amrex::ParticleReal>::quiet_NaN();
        amrex::ParticleReal derror = std::numeric_limits<amrex::ParticleReal>::quiet_NaN();
        amrex::ParticleReal error = 0.0;
        evaluate_lattice(&sim, &error, &q1_k, &q2_k);

#elif MYMODE == 1
        // forward differentiable run
        amrex::Print() << "Forward differentiable run\n";
        amrex::ParticleReal error = 0.0;
        amrex::ParticleReal derror = 0.0;

        // in forward-mode, usually use one variable at a time, but enzyme_width allows batching
        // two first order derivatives.
        //   note: seeded direction, this is AD and NOT a finite difference
        amrex::ParticleReal dq1_k = 1.0;  // derror += 1.0 * df/dq1_k
        amrex::ParticleReal dq2_k = 1.0;  // derror += 1.0 * df/dq1_k
        amrex::ParticleReal zero = 0.0;  // derror += 0.0 * df/dq*_k
        __enzyme_fwddiff(
            &evaluate_lattice,
            enzyme_const, &sim,
            enzyme_width, 2,  // 2x batched: first df/dq1_k then df/dq2_k
            // variable, batch 1, batch 2:
            enzyme_dup, &error, &derror_dq1_k, &derror_dq2_k,
            enzyme_dup, &q1_k, &dq1_k, &zero,
            enzyme_dup, &q2_k, &zero, &dq2_k
        );

#elif MYMODE == 2
        // reverse differentiable run:
        amrex::Print() << "Reverse differentiable run\n";
        amrex::ParticleReal dq1_k = 0.0;  // accumulator, zero-init!
        amrex::ParticleReal dq2_k = 0.0;  // accumulator, zero-init!
        amrex::ParticleReal error = 0.0;
        amrex::ParticleReal derror = 1.0;  // variable to back-prop: derror * df/dx
        __enzyme_autodiff(
            &evaluate_lattice,
            enzyme_const, &sim,
            enzyme_dup, &error, &derror,
            enzyme_dup, &q1_k, &dq1_k,
            enzyme_dup, &q2_k, &dq2_k
        );
        derror_dq1_k = dq1_k;
        derror_dq2_k = dq2_k;
#endif

        amrex::Print() << "error = " << error << "\n"
                       << "derror_dq1_k = " << derror_dq1_k << "\n"
                       << "derror_dq2_k = " << derror_dq2_k << "\n"
                       << std::endl;

        sim.finalize();
    }
} // namespace impactx
