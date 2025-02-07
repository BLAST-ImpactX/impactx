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
#include "particles/ImpactXParticleContainer.H"
#include "particles/Push.H"
#include "diagnostics/DiagnosticOutput.H"
#include "particles/wakefields/HandleWakefield.H"

#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <memory>


namespace impactx
{
    void
    ImpactX::track_envelope ()
    {
        BL_PROFILE("ImpactX::track_envelope");

        // verbosity
        amrex::ParmParse pp_impactx("impactx");
        int verbose = 1;
        pp_impactx.queryAddWithParser("verbose", verbose);

        // a global step for diagnostics including space charge slice steps in elements
        //   before we start the evolve loop, we are in "step 0" (initial state)
        int step = 0;

        // check typos in inputs after step 1
        bool early_params_checked = false;

        // access beam data
        if (!amr_data->track_envelope.m_ref.has_value()) {
            throw std::runtime_error("track_envelope: Reference particle not set.");
        }
        if (!amr_data->track_envelope.m_cm.has_value()) {
            throw std::runtime_error("track_envelope: Envelope (covariance matrix) not set.");
        }
        auto & ref = amr_data->track_envelope.m_ref.value();
        auto & cm = amr_data->track_envelope.m_cm.value();

        // output of init state
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
            diagnostics::DiagnosticOutput(ref, "diags/ref_particle");

            // print the initial values of reduced beam characteristics
            diagnostics::DiagnosticOutput(cm, ref, "diags/reduced_beam_characteristics");

        }

        amrex::ParmParse const pp_algo("algo");
        bool space_charge = false;
        pp_algo.query("space_charge", space_charge);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!space_charge, "Space charge not yet implemented for envelope tracking.");
        if (verbose > 0) {
            amrex::Print() << " Space Charge effects: " << space_charge << "\n";
        }

        bool csr = false;
        pp_algo.query("csr", csr);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!csr, "CSR not yet implemented for envelope tracking.");
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
                ref.sedge = ref.s;

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

                    std::visit([&ref, &cm](auto&& element){
                        // push reference particle in global coordinates
                        {
                            BL_PROFILE("impactx::Push::RefPart");
                            element(ref);
                        }

                        // push Covariance Matrix
                        element(cm, ref);

                    }, element_variant);

                    // just prints an empty newline at the end of the slice_step
                    if (verbose > 0) {
                        amrex::Print() << "\n";
                    }

                    // slice-step diagnostics
                    bool slice_step_diagnostics = false;
                    pp_diag.queryAdd("slice_step_diagnostics", slice_step_diagnostics);


                    if (diag_enable && slice_step_diagnostics) {
                        // print slice step reference particle to file
                        diagnostics::DiagnosticOutput(ref, "diags/ref_particle", step, true);

                        // print slice step reduced beam characteristics to file
                        diagnostics::DiagnosticOutput(cm, ref, "diags/reduced_beam_characteristics", step, true);

                    }

                    // inputs: unused parameters (e.g. typos) check after step 1 has finished
                    if (!early_params_checked) { early_params_checked = early_param_check(); }

                } // end in-element space-charge slice-step loop

            } // end beamline element loop

        } // end periods though the lattice loop

        if (diag_enable)
        {
            // print final reference particle to file
            diagnostics::DiagnosticOutput(ref, "diags/ref_particle_final", step);

            // print the final values of the reduced beam characteristics
            diagnostics::DiagnosticOutput(cm, ref, "diags/reduced_beam_characteristics_final", step);

        }

        // loop over all beamline elements & finalize them
        finalize_elements();
    }
} // namespace impactx
