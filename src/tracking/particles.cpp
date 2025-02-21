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
#include "particles/ImpactXParticleContainer.H"
#include "particles/Push.H"
#include "diagnostics/DiagnosticOutput.H"
#include "particles/spacecharge/HandleSpacecharge.H"
#include "particles/wakefields/HandleWakefield.H"

#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

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
        //   before we start the evolve loop, we are in "step 0" (initial state)
        int step = 0;

        // check typos in inputs after step 1
        bool early_params_checked = false;

        amrex::ParmParse pp_diag("diag");
        bool diag_enable = true;
        pp_diag.queryAdd("enable", diag_enable);
        if (verbose > 0) {
            amrex::Print() << " Diagnostics: " << diag_enable << "\n";
        }

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
            amrex::Print() << " Space Charge effects: " << amrex::getEnumNameString(space_charge) << "\n";
        }
        if (space_charge == SpaceChargeAlgo::True_2D)
        {
            throw std::runtime_error("2D space charge effects are not yet implemented for particle tracking.");
        }

        amrex::ParmParse const pp_algo("algo");
        bool csr = false;
        pp_algo.query("csr", csr);
        if (verbose > 0) {
            amrex::Print() << " CSR effects: " << csr << "\n";
        }

        // periods through the lattice
        int num_periods = 1;
        amrex::ParmParse("lattice").queryAddWithParser("periods", num_periods);

        for (int period=0; period < num_periods; ++period) {
            // loop over all beamline elements
            for (auto &element_variant: m_lattice) {
                // update element edge of the reference particle
                amr_data->track_particles.m_particle_container->SetRefParticleEdge();

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
                        amrex::Print() << " ++++ Starting step=" << step
                                       << " slice_step=" << slice_step << "\n";
                    }

                    // Wakefield calculation: call wakefield function to apply wake effects
                    particles::wakefields::HandleWakefield(*amr_data->track_particles.m_particle_container, element_variant, slice_ds);

                    // Space-charge calculation
                    particles::spacecharge::HandleSpacecharge(amr_data, [this](){ this->ResizeMesh(); }, slice_ds);

                    // push all particles with external maps
                    Push(*amr_data->track_particles.m_particle_container, element_variant, step, period);

                    // move "lost" particles to another particle container
                    collect_lost_particles(*amr_data->track_particles.m_particle_container);

                    // just prints an empty newline at the end of the slice_step
                    if (verbose > 0) {
                        amrex::Print() << "\n";
                    }

                    // slice-step diagnostics
                    bool slice_step_diagnostics = false;
                    pp_diag.queryAdd("slice_step_diagnostics", slice_step_diagnostics);

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

            } // end beamline element loop
        } // end periods though the lattice loop

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

        // loop over all beamline elements & finalize them
        finalize_elements();
    }
} // namespace impactx
