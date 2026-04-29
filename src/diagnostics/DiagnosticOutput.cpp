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

#include <cstdint>
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

            // determine whether to output spin moments
            amrex::ParmParse pp_algo("algo");
            bool compute_spin_moments = false;
            pp_algo.queryAdd("spin", compute_spin_moments);

            file_handler << "step" << " " << "s" << " "
                         << "mean_x" << " " << "min_x" << " " << "max_x" << " "
                         << "mean_y" << " " << "min_y" << " " << "max_y" << " "
                         << "mean_t" << " " << "min_t" << " " << "max_t" << " "
                         << "sigma_x" << " " << "sigma_y" << " " << "sigma_t" << " "
                         << "mean_px" << " " << "min_px" << " " << "max_px" << " "
                         << "mean_py" << " " << "min_py" << " " << "max_py" << " "
                         << "mean_pt" << " " << "min_pt" << " " << "max_pt" << " "
                         << "sigma_px" << " " << "sigma_py" << " " << "sigma_pt" << " "
                         << "emittance_x" << " " << "emittance_y" << " " << "emittance_t" << " "
                         << "alpha_x" << " " << "alpha_y" << " " << "alpha_t" << " "
                         << "beta_x" << " " << "beta_y" << " " << "beta_t" << " "
                         << "dispersion_x" << " " << "dispersion_px" << " "
                         << "dispersion_y" << " " << "dispersion_py" << " "
                         << "emittance_xn" << " " << "emittance_yn" << " " << "emittance_tn";
            if (compute_eigenemittances)
            {
                file_handler << " "
                             << "emittance_1" << " " << "emittance_2" << " " << "emittance_3";
            }
            if (compute_spin_moments)
            {
                file_handler << " "
                             << "mean_sx" << " " << "mean_sy" << " " << "mean_sz" << " "
                             << "sigma_sx" << " " << "sigma_sy" << " " << "sigma_sz";
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
        int64_t step
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
        int64_t step
    )
    {
        // determine whether to output eigenemittances
        amrex::ParmParse pp_diag("diag");
        bool compute_eigenemittances = false;
        pp_diag.queryAdd("eigenemittances", compute_eigenemittances);

        // determine whether to output spin moments
        amrex::ParmParse pp_algo("algo");
        bool compute_spin_moments = false;
        pp_algo.queryAdd("spin", compute_spin_moments);

        file_handler << step << " " << s << " "
                << rbc.at("mean_x") << " " << rbc.at("min_x") << " " << rbc.at("max_x") << " "
                << rbc.at("mean_y") << " " << rbc.at("min_y") << " " << rbc.at("max_y") << " "
                << rbc.at("mean_t") << " " << rbc.at("min_t") << " " << rbc.at("max_t") << " "
                << rbc.at("sigma_x") << " " << rbc.at("sigma_y") << " " << rbc.at("sigma_t") << " "
                << rbc.at("mean_px") << " " << rbc.at("mean_px") << " " << rbc.at("mean_px") << " "
                << rbc.at("mean_py") << " " << rbc.at("mean_py") << " " << rbc.at("mean_py") << " "
                << rbc.at("mean_pt") << " " << rbc.at("mean_pt") << " " << rbc.at("mean_pt") << " "
                << rbc.at("sigma_px") << " " << rbc.at("sigma_py") << " " << rbc.at("sigma_pt") << " "
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
        if (compute_spin_moments) {
            file_handler << " "
                         << rbc.at("mean_sx") << " " << rbc.at("mean_sy") << " " << rbc.at("mean_sz") << " "
                         << rbc.at("sigma_sx") << " " << rbc.at("sigma_sy") << " " << rbc.at("sigma_sz");

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
        int64_t step,
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
        int64_t step,
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
        int64_t step,
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
