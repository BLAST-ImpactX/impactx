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
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>

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

    amrex::ParticleReal
    evaluate_lattice (ImpactX * sim, amrex::ParticleReal const & q1_k)  //, amrex::ParticleReal q2_k)
    {
        amrex::ParticleReal q2_k = 3.5;

        // ns = 10  // TODO: number of slices per ds in the element

        auto dr1 = elements::Drift{2.7};
        auto q1 = elements::Quad{0.1, q1_k};
        auto dr2 = elements::Drift{1.4};
        auto q2 = elements::Quad{0.2, q2_k};
        auto dr3 = elements::Drift{1.4};
        auto q3 = elements::Quad{0.1, q1_k};
        auto dr4 = elements::Drift{2.7};

        /*
        // quadrupole triplet
        // https://impactx.readthedocs.io/en/latest/usage/examples/optimize_triplet/README.html
        active_sim->m_lattice = {
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

        //active_sim->evolve();
        //auto & pc = *active_sim->amr_data->m_particle_container;
        //dr1(pc, 0, 0);

        // envelope mode
        //std::cout << "before ref\n";
        auto ref = sim->amr_data->track_envelope.m_ref.value();
        //std::cout << "before env\n";
        auto env = sim->amr_data->track_envelope.m_env.value();
        //std::cout << "before cm\n";
        auto cm = env.m_env;

        auto & intensity = env.m_beam_intensity;
        intensity = 1.0e-6;  // a big value

        // push reference particle in global coordinates
        dr1(ref);
        // push Covariance Matrix in external fields
        dr1(cm, ref);

        q1(ref);
        q1(cm, ref);

        dr2(ref);
        dr2(cm, ref);

        q2(ref);
        q2(cm, ref);

        envelope::spacecharge::space_charge2D_push(ref, cm, intensity, q2.ds());

        dr3(ref);
        dr3(cm, ref);

        envelope::spacecharge::space_charge3D_push(ref, cm, intensity, dr3.ds());

        q3(ref);
        q3(cm, ref);

        dr4(ref);
        dr4(cm, ref);

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
        return diagnostics::reduced_beam_characteristics(cm, ref);


        // not alpha_x but instead a simple x_ms * y_ms * t_ms
        //return cm(1,1) * cm(3,3) * cm(5,5);
    }

    void my_run ()
    {
        ImpactX sim;

        amrex::ParmParse pp_algo("algo");
        std::string track = "envelope";
        pp_algo.add("track", track);

        sim.init_grids();

        // TODO: replace with beam params from https://impactx.readthedocs.io/en/latest/usage/examples/optimize_triplet/README.html
        sim.initBeamDistributionFromInputs();

        // design the accelerator lattice
        // sim.initLatticeElementsFromInputs();

        // initial quad strengths
        amrex::ParticleReal q1_k = -3.0;
        //amrex::ParticleReal q2_k = 3.0;
        amrex::Print() << "q1_k = " << q1_k << std::endl;

#define MYMODE 2

#if MYMODE == 0
        // non-differentiable run:
        amrex::ParticleReal ddx = std::numeric_limits<amrex::ParticleReal>::quiet_NaN();
        amrex::ParticleReal const alpha_x = evaluate_lattice(&sim, q1_k);
        amrex::Print() << "final alpha_x = " << alpha_x << std::endl;

#elif MYMODE == 1
        // forward differentiable run
        //   note: seeded direction, this is AD and NOT a finite difference
        amrex::ParticleReal dq1_k = 1.0;
        amrex::ParticleReal ddx = __enzyme_fwddiff(
            &evaluate_lattice,
            enzyme_const, &sim,
            enzyme_dup, &q1_k, &dq1_k
        );

#elif MYMODE == 2
        // reverse differentiable run:
        amrex::ParticleReal ddx = 0.0;  // accumulator, zero-init!
        __enzyme_autodiff(
            &evaluate_lattice,
            enzyme_const, &sim,
            enzyme_dup, &q1_k, &ddx
        );
#endif

        amrex::Print() << "ddx = " << ddx << std::endl;

        sim.finalize();
    }
} // namespace impactx
