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
    /** Type of beam diagnostic output
     */
    enum class OutputType
    {
        PrintRefParticle, ///< ASCII diagnostics
        PrintReducedBeamCharacteristics ///< ASCII diagnostics, for beam momenta and Twiss parameters
    };

    void
    write_column_header (
        amrex::AllPrintToFile & file_handler,
        OutputType otype
    )
    {
        if (otype == OutputType::PrintRefParticle)
        {
            file_handler << "step s beta gamma beta_gamma x y z t px py pz pt\n";
        }
        else if (otype == OutputType::PrintReducedBeamCharacteristics)
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
            file_handler << " " << "charge_C"
                         << "\n";
        }
    }

    void
    prepare_header (
        amrex::AllPrintToFile & file_handler,
        OutputType otype,
        bool append
    )
    {
        file_handler.SetPrecision(std::numeric_limits<amrex::ParticleReal>::max_digits10);

        // write file header per MPI RANK
        if (!append)
        {
            write_column_header(file_handler, otype);
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
                << rbc.at("x_mean") << " " << rbc.at("x_min") << " " << rbc.at("x_max") << " "
                << rbc.at("y_mean") << " " << rbc.at("y_min") << " " << rbc.at("y_max") << " "
                << rbc.at("t_mean") << " " << rbc.at("t_min") << " " << rbc.at("t_max") << " "
                << rbc.at("sig_x") << " " << rbc.at("sig_y") << " " << rbc.at("sig_t") << " "
                << rbc.at("px_mean") << " " << rbc.at("px_min") << " " << rbc.at("px_max") << " "
                << rbc.at("py_mean") << " " << rbc.at("py_min") << " " << rbc.at("py_max") << " "
                << rbc.at("pt_mean") << " " << rbc.at("pt_min") << " " << rbc.at("pt_max") << " "
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
        file_handler << " " << rbc.at("charge_C")
                     << "\n";
    }
}


namespace impactx::diagnostics
{
    void DiagnosticOutput (
        ImpactXParticleContainer const & pc,
        std::string file_name,
        int step,
        bool append
    )
    {
        BL_PROFILE("impactx::diagnostics::DiagnosticOutput(pc)");

        using namespace amrex::literals; // for _rt and _prt

        OutputType const otype = OutputType::PrintReducedBeamCharacteristics;

        // keep file open as we add more and more lines
        amrex::AllPrintToFile file_handler(std::move(file_name));
        prepare_header(file_handler, otype, append);

        amrex::ParticleReal const s = pc.GetRefParticle().s;
        std::unordered_map<std::string, amrex::ParticleReal> const rbc =
            diagnostics::reduced_beam_characteristics(pc);

        write_rbc(file_handler, rbc, s, step);
    }

    void DiagnosticOutput (
        Map6x6 const & cm,
        RefPart const & ref_part,
        std::string file_name,
        int step,
        bool append
    )
    {
        BL_PROFILE("impactx::diagnostics::DiagnosticOutput(cm)");

        // keep file open as we add more and more lines
        amrex::AllPrintToFile file_handler(std::move(file_name));
        prepare_header(file_handler, OutputType::PrintReducedBeamCharacteristics, append);

        amrex::ParticleReal const s = ref_part.s;
        std::unordered_map<std::string, amrex::ParticleReal> const rbc =
            diagnostics::reduced_beam_characteristics(cm, ref_part);

        write_rbc(file_handler, rbc, s, step);
    }

    void DiagnosticOutput (
        RefPart const & ref_part,
        std::string file_name,
        int step,
        bool append
    )
    {
        BL_PROFILE("impactx::diagnostics::DiagnosticOutput(pc)");

        OutputType const otype = OutputType::PrintRefParticle;

        // keep file open as we add more and more lines
        amrex::AllPrintToFile file_handler(std::move(file_name));
        prepare_header(file_handler, otype, append);
        write_ref(file_handler, ref_part, step);
    }

} // namespace impactx::diagnostics
