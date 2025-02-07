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

#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <memory>


namespace impactx
{
    void
    ImpactX::track_reference (RefPart & ref)
    {
        BL_PROFILE("ImpactX::track_reference");

        // verbosity
        amrex::ParmParse pp_impactx("impactx");
        int verbose = 1;
        pp_impactx.queryAddWithParser("verbose", verbose);

        // a global step for diagnostics including space charge slice steps in elements
        //   before we start the evolve loop, we are in "step 0" (initial state)
        int step = 0;

        // check typos in inputs after step 1
        bool early_params_checked = false;

        // output of init state
        amrex::ParmParse pp_diag("diag");
        bool diag_enable = true;
        pp_diag.queryAdd("enable", diag_enable);
        if (verbose > 0)
        {
            amrex::Print() << " Diagnostics: " << diag_enable << "\n";
        }

        if (diag_enable)
        {
            int file_min_digits = 6;
            pp_diag.queryAddWithParser("file_min_digits", file_min_digits);

            // print initial reference particle to file
            diagnostics::DiagnosticOutput(ref, "diags/ref_particle");

        }

        // periods through the lattice
        int num_periods = 1;
        amrex::ParmParse("lattice").queryAddWithParser("periods", num_periods);

        for (int period=0; period < num_periods; ++period)
        {
            // loop over all beamline elements
            for (auto &element_variant: m_lattice) {
                // update element edge of the reference particle
                ref.sedge = ref.s;

                // number of slices through the element
                int nslice = 1;
                amrex::ParticleReal slice_ds; // in meters
                std::visit([&nslice, &slice_ds](auto &&element)
                {
                    nslice = element.nslice();
                    slice_ds = element.ds() / nslice;
                }, element_variant);

                // sub-steps within the element
                for (int slice_step = 0; slice_step < nslice; ++slice_step)
                {
                    BL_PROFILE("ImpactX::track_reference::slice_step");
                    step++;
                    if (verbose > 0) {
                        amrex::Print() << " ++++ Starting step=" << step
                                       << " slice_step=" << slice_step << "\n";
                    }

                    std::visit([&ref](auto&& element)
                    {
                        // push reference particle in global coordinates
                        BL_PROFILE("impactx::Push::RefPart");
                        element(ref);
                    }, element_variant);

                    // just prints an empty newline at the end of the slice_step
                    if (verbose > 0)
                    {
                        amrex::Print() << "\n";
                    }

                    // slice-step diagnostics
                    bool slice_step_diagnostics = false;
                    pp_diag.queryAdd("slice_step_diagnostics", slice_step_diagnostics);


                    if (diag_enable && slice_step_diagnostics)
                    {
                        // print slice step reference particle to file
                        diagnostics::DiagnosticOutput(ref, "diags/ref_particle", step, true);
                    }

                    // inputs: unused parameters (e.g. typos) check after step 1 has finished
                    if (!early_params_checked) { early_params_checked = early_param_check(); }

                } // end in-element slice-step loop

            } // end beamline element loop

        } // end periods though the lattice loop

        if (diag_enable)
        {
            // print final reference particle to file
            diagnostics::DiagnosticOutput(ref, "diags/ref_particle_final", step);
        }

        // loop over all beamline elements & finalize them
        finalize_elements();
    }
} // namespace impactx
