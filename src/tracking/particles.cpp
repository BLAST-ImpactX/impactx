/* Copyright 2022-2025 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactX.H"
#include "initialization/Algorithms.H"
#include "initialization/InitAmrCore.H"
#include "particles/CollectLost.H"
#include "particles/ParticleBoundary.H"
#include "particles/ImpactXParticleContainer.H"
#include "particles/Push.H"
#include "diagnostics/DiagnosticOutput.H"
#include "particles/spacecharge/HandleSpacecharge.H"
#include "particles/wakefields/HandleWakefield.H"
#include "particles/wakefields/HandleISR.H"

#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <cstdint>
#include <memory>


namespace impactx
{
    void ImpactX::track_particles ()
    {
        BL_PROFILE("ImpactX::track_particles");

        validate();

        // verbosity
        amrex::ParmParse pp_impactx("impactx");
        int verbose = 1;
        pp_impactx.queryAddWithParser("verbose", verbose);

        // a global step for diagnostics including space charge slice steps in elements
        //   before we start the tracking loop, we are in "step 0" (initial state)
        int64_t & step = m_tracking_state.m_step;
        step = 0;

        // period in the lattice (e.g., turns)
        int & period = m_tracking_state.m_period;
        period = 0;

        // check typos in inputs after step 1
        bool early_params_checked = false;

        // shortcuts
        auto & pc = amr_data->track_particles.m_particle_container;

        // diags
        amrex::ParmParse pp_diag("diag");
        bool diag_enable = true;
        pp_diag.queryAdd("enable", diag_enable);
        if (verbose > 0) {
            amrex::Print() << " Diagnostics: " << diag_enable << "\n";
        }

        pc->reset_beam_moments_history();

        if (diag_enable)
        {
            int file_min_digits = 6;
            pp_diag.queryAddWithParser("file_min_digits", file_min_digits);

            // print initial reference particle to file
            diagnostics::DiagnosticOutput(amr_data->track_particles.m_particle_container->GetRefParticle(),
                                          "diags/ref_particle",
                                          step);

            // print the initial values of reduced beam characteristics
            diagnostics::DiagnosticOutput(*amr_data->track_particles.m_particle_container,
                                          "diags/reduced_beam_characteristics");

        }

        auto space_charge = get_space_charge_algo();
        if (verbose > 0) {
            amrex::Print() << " Space Charge effects: " << to_string(space_charge) << "\n";
        }

        amrex::ParmParse const pp_algo("algo");
        bool csr = false;
        pp_algo.query("csr", csr);
        bool isr = false;
        pp_algo.query("isr", isr);
        bool spin = false;
        pp_algo.query("spin", spin);

        if (verbose > 0) {
            amrex::Print() << " CSR effects: " << csr << "\n";
            amrex::Print() << " ISR effects: " << isr << "\n";
            amrex::Print() << " Spin tracking: " << spin << "\n";
        }

        // periods through the lattice
        int num_periods = 1;
        amrex::ParmParse("lattice").queryAddWithParser("periods", num_periods);

        for (period=0; period < num_periods; ++period) {
            // optional, user-defined function call
            m_tracking_state.m_element = &m_lattice.front();
            call_hook("before_period");

            // loop over all beamline elements
            for (auto &element_variant: m_lattice) {
                // update element edge of the reference particle
                amr_data->track_particles.m_particle_container->SetRefParticleEdge();

                // optional, user-defined function call
                m_tracking_state.m_element = &element_variant;
                call_hook("before_element");

                // number of slices used for the application of space charge
                int nslice = 1;
                amrex::ParticleReal slice_ds; // in meters
                std::visit([&nslice, &slice_ds](auto &&element) {
                    nslice = element.nslice();
                    slice_ds = element.ds() / nslice;
                }, element_variant);

                // sub-steps for space charge within the element
                for (int slice_step = 0; slice_step < nslice; ++slice_step) {
                    BL_PROFILE("ImpactX::evolve::slice_step");
                    step++;
                    if (verbose > 0) {
                        amrex::Print() << "\n++++ Starting step=" << step
                                       << " slice_step=" << slice_step;
                    }

                    // optional, user-defined function call
                    call_hook("before_slice");

                    // Wakefield calculation: call wakefield function to apply wake effects
                    particles::wakefields::HandleWakefield(*amr_data->track_particles.m_particle_container, element_variant, slice_ds);

                    // ISR calculation: call ISR function to apply incoherent synchrotron radiation effects
                    particles::wakefields::HandleISR(*amr_data->track_particles.m_particle_container, element_variant, slice_ds);

                    // Space-charge calculation
                    particles::spacecharge::HandleSpacecharge(amr_data, [this](){ this->ResizeMesh(); }, slice_ds);

                    // push all particles with external maps
                    push(*amr_data->track_particles.m_particle_container, element_variant, step, period);

                    // Apply optional particle boundary conditions
                    particles::ParticleBoundary(*amr_data->track_particles.m_particle_container);

                    // move "lost" particles to another particle container
                    collect_lost_particles(*amr_data->track_particles.m_particle_container);

                    // just prints an empty newline at the end of the slice_step
                    if (verbose > 0) {
                        amrex::Print() << "\n";
                    }

                    // slice-step diagnostics
                    bool slice_step_diagnostics = false;
                    pp_diag.queryAdd("slice_step_diagnostics", slice_step_diagnostics);

                    if (amr_data->track_particles.m_particle_container->store_beam_moments) {
                        amr_data->track_particles.m_particle_container->record_beam_moments();
                    }

                    if (diag_enable && slice_step_diagnostics) {
                        // print slice step reference particle to file
                        diagnostics::DiagnosticOutput(amr_data->track_particles.m_particle_container->GetRefParticle(),
                                                      "diags/ref_particle",
                                                      step,
                                                      true);

                        // print slice step reduced beam characteristics to file
                        diagnostics::DiagnosticOutput(*amr_data->track_particles.m_particle_container,
                                                      "diags/reduced_beam_characteristics",
                                                      step,
                                                      true);

                    }

                    // inputs: unused parameters (e.g. typos) check after step 1 has finished
                    if (!early_params_checked) { early_params_checked = early_param_check(); }

                } // end in-element space-charge slice-step loop

                // optional, user-defined function call
                call_hook("after_element");

            } // end beamline element loop

            // optional, user-defined function call
            call_hook("after_period");
        } // end periods though the lattice loop

        // avoid dangling references if users manipulate the lattice
        m_tracking_state.set_no_element();

        if (diag_enable)
        {
            // print final reference particle to file
            diagnostics::DiagnosticOutput(amr_data->track_particles.m_particle_container->GetRefParticle(),
                                          "diags/ref_particle_final",
                                          step);

            // print the final values of the reduced beam characteristics
            diagnostics::DiagnosticOutput(*amr_data->track_particles.m_particle_container,
                                          "diags/reduced_beam_characteristics_final",
                                          step);

            // output particles lost in apertures
            if (amr_data->track_particles.m_particles_lost->TotalNumberOfParticles() > 0)
            {
                std::string openpmd_backend = "default";
                pp_diag.queryAdd("backend", openpmd_backend);

                elements::diagnostics::BeamMonitor output_lost("particles_lost", openpmd_backend, "g");
                output_lost(*amr_data->track_particles.m_particles_lost, 0, 0);
                output_lost.finalize();
            }
        }
    }
} // namespace impactx
