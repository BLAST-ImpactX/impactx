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
#include "particles/ImpactXParticleContainer.H"
#include "particles/Push.H"
#include "envelope/spacecharge/EnvelopeSpaceChargePush.H"
#include "diagnostics/DiagnosticOutput.H"

#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>

#include <memory>
#include <stdexcept>


namespace impactx
{
    void
    ImpactX::track_envelope ()
    {
        BL_PROFILE("ImpactX::track_envelope");

        using namespace amrex::literals;

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
        if (!amr_data->track_envelope.m_ref.has_value())
        {
            throw std::runtime_error("track_envelope: Reference particle not set.");
        }
        if (!amr_data->track_envelope.m_env.has_value())
        {
            throw std::runtime_error("track_envelope: Envelope (covariance matrix) not set.");
        }
        auto & ref = amr_data->track_envelope.m_ref.value();
        auto & env = amr_data->track_envelope.m_env.value();
        auto & cm = env.m_env;
        auto & intensity = env.m_beam_intensity;

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
        auto space_charge = get_space_charge_algo();
        if (verbose > 0)
        {
            amrex::Print() << " Space Charge effects: " << to_string(space_charge) << "\n";
        }
        if (space_charge == SpaceChargeAlgo::True_3D && intensity == 0_prt) {
            ablastr::warn_manager::WMRecordWarning(
                "algo.space_charge",
                "Space charge calculations are enabled but zero bunch charge was provided. "
                "Skipping space charge calculations.",
                ablastr::warn_manager::WarnPriority::high
            );
        }
        if (space_charge == SpaceChargeAlgo::True_2D && intensity == 0_prt) {
            ablastr::warn_manager::WMRecordWarning(
                "algo.space_charge",
                "Space charge calculations are enabled but zero beam current was provided. "
                "Skipping space charge calculations.",
                ablastr::warn_manager::WarnPriority::high
            );
        }

        bool csr = false;
        pp_algo.query("csr", csr);
        if (verbose > 0)
        {
            amrex::Print() << " CSR effects: " << csr << "\n";
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!csr, "CSR effects are not yet implemented for envelope tracking.");

        // periods through the lattice
        int num_periods = 1;
        amrex::ParmParse("lattice").queryAddWithParser("periods", num_periods);

        for (int period=0; period < num_periods; ++period)
        {
            // loop over all beamline elements
            for (auto &element_variant: m_lattice)
            {
                // update element edge of the reference particle
                ref.sedge = ref.s;

                // number of slices used for the application of space charge
                int nslice = 1;
                amrex::ParticleReal slice_ds; // in meters
                std::visit([&nslice, &slice_ds](auto &&element)
                {
                    nslice = element.nslice();
                    slice_ds = element.ds() / nslice;
                }, element_variant);

                // sub-steps for space charge within the element
                int nsteps = std::floor(nslice/2.0);
                for (int slice_step = 0; slice_step < nsteps; ++slice_step)
                {
                    BL_PROFILE("ImpactX::track_envelope::slice_step");
                    step++;
                    if (verbose > 0)
                    {
                        amrex::Print() << " ++++ Starting step=" << step
                                       << " slice_step=" << slice_step << "\n";
                    }

                    std::visit([&ref, &cm](auto&& element)
                    {
                        // push reference particle in global coordinates
                        {
                            BL_PROFILE("impactx::Push::RefPart");
                            element(ref);
                        }

                        // push Covariance Matrix in external fields
                        element(cm, ref);

                    }, element_variant);

                    if (space_charge == SpaceChargeAlgo::True_2D)
                    {
                        // push Covariance Matrix in 2D space charge fields
                        envelope::spacecharge::space_charge2D_push(ref, cm, intensity, 2*slice_ds);
                    } else if (space_charge == SpaceChargeAlgo::True_3D)
                    {
                        // push Covariance Matrix in 3D space charge fields
                        envelope::spacecharge::space_charge3D_push(ref, cm, intensity, 2*slice_ds);
                    } else {
                        amrex::Print() << "Warning: Space charge is off by default." << "\n";
                    }

                    std::visit([&ref, &cm](auto&& element)
                    {
                        // push reference particle in global coordinates
                        {
                            BL_PROFILE("impactx::Push::RefPart");
                            element(ref);
                        }

                        // push Covariance Matrix in external fields
                        element(cm, ref);

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
    }
} // namespace impactx
