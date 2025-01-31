/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#include "DiagnosticOutput.H"
#include "NonlinearLensInvariants.H"
#include "particles/CovarianceMatrix.H"
#include "ReducedBeamCharacteristics.H"

#include <AMReX_BLProfiler.H> // for BL_PROFILE
#include <AMReX_ParmParse.H>  // for ParmParse
#include <AMReX_REAL.H>       // for ParticleReal
#include <AMReX_Print.H>      // for PrintToFile

#include <limits>
#include <stdexcept>
#include <utility>


namespace
{
    void
    write_column_header (
        amrex::AllPrintToFile & file_handler,
        impactx::diagnostics::OutputType otype,
        bool has_charge
    )
    {
        if (otype == impactx::diagnostics::OutputType::PrintRefParticle)
        {
            file_handler << "step s beta gamma beta_gamma x y z t px py pz pt\n";
        }
        else if (otype == impactx::diagnostics::OutputType::PrintReducedBeamCharacteristics)
        {
            // determine whether to output eigenemittances
            amrex::ParmParse pp_diag("diag");
            bool compute_eigenemittances = false;
            pp_diag.queryAdd("eigenemittances", compute_eigenemittances);

            file_handler << "step" << " " << "s" << " "
                         << "x_mean" << " " << "x_min" << " " << "x_max" << " "
                         << "y_mean" << " " << "y_min" << " " << "y_max" << " "
                         << "t_mean" << " " << "t_min" << " " << "t_max" << " "
                         << "sig_x" << " " << "sig_y" << " " << "sig_t" << " "
                         << "px_mean" << " " << "px_min" << " " << "px_max" << " "
                         << "py_mean" << " " << "py_min" << " " << "py_max" << " "
                         << "pt_mean" << " " << "pt_min" << " " << "pt_max" << " "
                         << "sig_px" << " " << "sig_py" << " " << "sig_pt" << " "
                         << "emittance_x" << " " << "emittance_y" << " " << "emittance_t" << " "
                         << "alpha_x" << " " << "alpha_y" << " " << "alpha_t" << " "
                         << "beta_x" << " " << "beta_y" << " " << "beta_t" << " "
                         << "dispersion_x" << " " << "dispersion_px" << " "
                         << "dispersion_y" << " " << "dispersion_py" << " "
                         << "emittance_xn" << " " << "emittance_yn" << " " << "emittance_tn";
            if (compute_eigenemittances)
            {
                file_handler << " "
                             << "emittance_xn" << " " << "emittance_yn" << " " << "emittance_tn";
            }
            if (has_charge)
            {
                file_handler << " " << "charge_C";
            }
            file_handler << "\n";
        }
    }

    void
    prepare_header (
        amrex::AllPrintToFile & file_handler,
        impactx::diagnostics::OutputType otype,
        bool has_charge,
        bool append
    )
    {
        file_handler.SetPrecision(std::numeric_limits<amrex::ParticleReal>::max_digits10);

        // write file header per MPI RANK
        if (!append)
        {
            write_column_header(file_handler, otype, has_charge);
        }
    }

    void
    write_ref (
        amrex::AllPrintToFile & file_handler,
        impactx::RefPart const & ref_part,
        int step
    )
    {
        amrex::ParticleReal const s = ref_part.s;
        amrex::ParticleReal const beta = ref_part.beta();
        amrex::ParticleReal const gamma = ref_part.gamma();
        amrex::ParticleReal const beta_gamma = ref_part.beta_gamma();
        amrex::ParticleReal const x = ref_part.x;
        amrex::ParticleReal const y = ref_part.y;
        amrex::ParticleReal const z = ref_part.z;
        amrex::ParticleReal const t = ref_part.t;
        amrex::ParticleReal const px = ref_part.px;
        amrex::ParticleReal const py = ref_part.py;
        amrex::ParticleReal const pz = ref_part.pz;
        amrex::ParticleReal const pt = ref_part.pt;

        // write particle data to file
        file_handler
                << step << " " << s << " "
                << beta << " " << gamma << " " << beta_gamma << " "
                << x << " " << y << " " << z << " " << t << " "
                << px << " " << py << " " << pz << " " << pt << "\n";
    }

    void
    write_rbc (
        amrex::AllPrintToFile & file_handler,
        std::unordered_map<std::string, amrex::ParticleReal> const & rbc,
        amrex::ParticleReal s,
        int step
    )
    {
        // determine whether to output eigenemittances
        amrex::ParmParse pp_diag("diag");
        bool compute_eigenemittances = false;
        pp_diag.queryAdd("eigenemittances", compute_eigenemittances);

        file_handler << step << " " << s << " "
                     << rbc.at("sig_x") << " " << rbc.at("sig_y") << " " << rbc.at("sig_t") << " "
                     << rbc.at("sig_px") << " " << rbc.at("sig_py") << " " << rbc.at("sig_pt") << " "
                     << rbc.at("emittance_x") << " " << rbc.at("emittance_y") << " " << rbc.at("emittance_t") << " "
                     << rbc.at("alpha_x") << " " << rbc.at("alpha_y") << " " << rbc.at("alpha_t") << " "
                     << rbc.at("beta_x") << " " << rbc.at("beta_y") << " " << rbc.at("beta_t") << " "
                     << rbc.at("dispersion_x") << " " << rbc.at("dispersion_px") << " "
                     << rbc.at("dispersion_y") << " " << rbc.at("dispersion_py") << " "
                     << rbc.at("emittance_xn") << " " << rbc.at("emittance_yn") << " " << rbc.at("emittance_tn");
        if (compute_eigenemittances) {
            file_handler << " "
                         << rbc.at("emittance_1") << " " << rbc.at("emittance_2") << " " << rbc.at("emittance_3");

        }
        file_handler << "\n";
    }
}


namespace impactx::diagnostics
{
    void DiagnosticOutput (
        ImpactXParticleContainer const & pc,
        OutputType const otype,
        std::string file_name,
        int step,
        bool append
    )
    {
        BL_PROFILE("impactx::diagnostics::DiagnosticOutput(pc)");

        using namespace amrex::literals; // for _rt and _prt

        // keep file open as we add more and more lines
        amrex::AllPrintToFile file_handler(std::move(file_name));
        prepare_header(file_handler, otype, true, append);

        if (otype == OutputType::PrintRefParticle)
        {
            RefPart const ref_part = pc.GetRefParticle();
            write_ref(file_handler, ref_part, step);
        }

        else if (otype == OutputType::PrintReducedBeamCharacteristics)
        {
            amrex::ParticleReal const s = pc.GetRefParticle().s;
            std::unordered_map<std::string, amrex::ParticleReal> const rbc =
                diagnostics::reduced_beam_characteristics(pc);

            write_rbc(file_handler, rbc, s, step);
        }

        if (otype == OutputType::PrintNonlinearLensInvariants)
        {
            // create a host-side particle buffer
            auto tmp = pc.make_alike<amrex::PinnedArenaAllocator>();

            // copy all particles from device to host
            bool const local = true;
            tmp.copyParticles(pc, local);

            // loop over refinement levels
            int const nLevel = tmp.finestLevel();
            for (int lev = 0; lev <= nLevel; ++lev) {
                // loop over all particle boxes
                using MyPinnedParIter = amrex::ParIterSoA<RealSoA::nattribs,IntSoA::nattribs, amrex::PinnedArenaAllocator>;
                for (MyPinnedParIter pti(tmp, lev); pti.isValid(); ++pti) {
                    const int np = pti.numParticles();

                    // preparing access to particle data: SoA of Reals
                    auto const& soa = pti.GetStructOfArrays();
                    auto const& part_x = soa.GetRealData(RealSoA::x);
                    auto const& part_y = soa.GetRealData(RealSoA::y);
                    auto const& part_px = soa.GetRealData(RealSoA::px);
                    auto const& part_py = soa.GetRealData(RealSoA::py);

                    auto const& part_idcpu = soa.GetIdCPUData().dataPtr();

                    // Parse the diagnostic parameters
                    amrex::ParmParse pp_diag("diag");

                    amrex::ParticleReal alpha = 0.0;
                    pp_diag.queryAddWithParser("alpha", alpha);

                    amrex::ParticleReal beta = 1.0;
                    pp_diag.queryAddWithParser("beta", beta);

                    amrex::ParticleReal tn = 0.4;
                    pp_diag.queryAddWithParser("tn", tn);

                    amrex::ParticleReal cn = 0.01;
                    pp_diag.queryAddWithParser("cn", cn);

                    NonlinearLensInvariants const nonlinear_lens_invariants(alpha, beta, tn, cn);

                    // print out particles (this hack works only on CPU and on GPUs with
                    // unified memory access)
                    for (int i = 0; i < np; ++i) {
                        amrex::ParticleReal const x = part_x[i];
                        amrex::ParticleReal const y = part_y[i];
                        uint64_t const global_id = part_idcpu[i];

                        amrex::ParticleReal const px = part_px[i];
                        amrex::ParticleReal const py = part_py[i];

                        // calculate invariants of motion
                        NonlinearLensInvariants::Data const HI_out =
                                nonlinear_lens_invariants(x, y, px, py);

                        // write particle invariant data to file
                        file_handler
                                << global_id << " "
                                << HI_out.H << " " << HI_out.I << "\n";

                    } // i=0...np
                } // end loop over all particle boxes
            } // end mesh-refinement level loop
        }
    }

    void DiagnosticOutput (
        Map6x6 const & cm,
        RefPart const & ref_part,
        OutputType otype,
        std::string file_name,
        int step,
        bool append
    )
    {
        BL_PROFILE("impactx::diagnostics::DiagnosticOutput(cm)");

        // keep file open as we add more and more lines
        amrex::AllPrintToFile file_handler(std::move(file_name));
        prepare_header(file_handler, otype, false, append);

        if (otype == OutputType::PrintRefParticle)
        {
            write_ref(file_handler, ref_part, step);
        }

        else if (otype == OutputType::PrintReducedBeamCharacteristics)
        {
            amrex::ParticleReal const s = ref_part.s;
            std::unordered_map<std::string, amrex::ParticleReal> const rbc =
                diagnostics::reduced_beam_characteristics(cm, ref_part);

            write_rbc(file_handler, rbc, s, step);
        }

        if (otype == OutputType::PrintNonlinearLensInvariants)
        {
            throw std::runtime_error("DiagnosticOutput(cm): PrintNonlinearLensInvariants is not implemented.");
        }
    }

} // namespace impactx::diagnostics
