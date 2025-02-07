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
#include "initialization/InitAmrCore.H"
#include "particles/CollectLost.H"
#include "particles/ImpactXParticleContainer.H"
#include "particles/Push.H"
#include "diagnostics/DiagnosticOutput.H"
#include "particles/spacecharge/ForceFromSelfFields.H"
#include "particles/spacecharge/GatherAndPush.H"
#include "particles/spacecharge/PoissonSolve.H"
#include "particles/transformation/CoordinateTransformation.H"
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

        amrex::ParmParse const pp_algo("algo");
        bool space_charge = false;
        pp_algo.query("space_charge", space_charge);
        if (verbose > 0) {
            amrex::Print() << " Space Charge effects: " << space_charge << "\n";
        }

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

                    // Space-charge calculation: turn off if there is only 1 particle
                    if (space_charge &&
                        amr_data->track_particles.m_particle_container->TotalNumberOfParticles(true, false)) {

                        // transform from x',y',t to x,y,z
                        transformation::CoordinateTransformation(
                                *amr_data->track_particles.m_particle_container,
                                CoordSystem::t);

                        // Note: The following operation assume that
                        // the particles are in x, y, z coordinates.

                        // Resize the mesh, based on `m_particle_container` extent
                        ResizeMesh();

                        // Redistribute particles in the new mesh in x, y, z
                        amr_data->track_particles.m_particle_container->Redistribute();

                        // charge deposition
                        amr_data->track_particles.m_particle_container->DepositCharge(
                            amr_data->track_particles.m_rho,
                            amr_data->refRatio()
                        );

                        // poisson solve in x,y,z
                        spacecharge::PoissonSolve(
                            *amr_data->track_particles.m_particle_container,
                            amr_data->track_particles.m_rho,
                            amr_data->track_particles.m_phi,
                            amr_data->refRatio()
                        );

                        // calculate force in x,y,z
                        spacecharge::ForceFromSelfFields(
                            amr_data->track_particles.m_space_charge_field,
                            amr_data->track_particles.m_phi,
                            amr_data->Geom()
                        );

                        // gather and space-charge push in x,y,z , assuming the space-charge
                        // field is the same before/after transformation
                        // TODO: This is currently using linear order.
                        spacecharge::GatherAndPush(
                            *amr_data->track_particles.m_particle_container,
                            amr_data->track_particles.m_space_charge_field,
                            amr_data->Geom(),
                            slice_ds
                        );

                        // transform from x,y,z to x',y',t
                        transformation::CoordinateTransformation(
                            *amr_data->track_particles.m_particle_container,
                            CoordSystem::s
                        );
                    }

                    // for later: original Impact implementation as an option
                    // Redistribute particles in x',y',t
                    //   TODO: only needed if we want to gather and push space charge
                    //         in x',y',t
                    //   TODO: change geometry beforehand according to transformation
                    //m_particle_container->Redistribute();
                    //
                    // in original Impact, we gather and space-charge push in x',y',t ,
                    // assuming that the distribution did not change

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
