/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell, Ji Qiang
 * License: BSD-3-Clause-LBNL
 */
#include "InitAMReX.H"

#include "initialization/InitParser.H"

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>

#if defined(AMREX_USE_MPI)
#   include <mpi.h>
#else
#   include <AMReX_ccse-mpi.H>
using amrex::mpidatatypes::MPI_COMM_WORLD;
#endif


namespace impactx::initialization
{
    void
    default_init_AMReX (int argc, char* argv[])
    {
        if (!amrex::Initialized())
        {
            bool const build_parm_parse = true;
            amrex::Initialize(
                    argc,
                    argv,
                    build_parm_parse,
                    MPI_COMM_WORLD,
                    impactx::initialization::overwrite_amrex_parser_defaults
            );
        }
    }

    void
    default_init_AMReX ()
    {
        // Pass a program name (argc = 1) so that AMReX runs ParmParse::Initialize.
        // With argc = 0 and build_parm_parse = true, AMReX skips ParmParse::Initialize
        // and never registers ParmParse::Finalize. ParmParse state (the global table)
        // would then leak across AMReX initialize/finalize cycles, e.g. in persistent
        // Jupyter notebook kernels that run several simulations in one process.
        int argc = 1;
        char arg0[] = "impactx";
        char* argv_storage[] = {arg0, nullptr};
        char** argv = argv_storage;
        default_init_AMReX(argc, argv);
    }
} // namespace impactx::initialization
