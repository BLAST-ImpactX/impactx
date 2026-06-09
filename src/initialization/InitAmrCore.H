/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACT_INIT_AMR_CORE_H
#define IMPACT_INIT_AMR_CORE_H

#include "AmrCoreData.H"

#include <limits>
#include <memory>


namespace impactx::initialization
{
    /** This creates AMReX boxes for each MPI-rank
     *
     * @param prec the runtime-selected beam precision to construct the data at
     * @return simulation_geometry the geometry (topology) of the simulation;
     *         amr_info contains information on mesh-refinement and box/grid blocks
     */
    std::unique_ptr<AmrCoreDataBase>
    init_amr_core (Precision prec);

    /** This lets AMReX AmrCore/AmrMesh build the boxes for each MPI-rank
     *
     * This uses the regular logic in AmrCore/AmrMesh to grid boxes and refinement levels.
     *
     * @param prec the runtime-selected beam precision to construct the data at
     * @return simulation_geometry the geometry (topology) of the simulation;
     *         amr_info contains information on mesh-refinement and box/grid blocks
     */
    std::unique_ptr<AmrCoreDataBase>
    amrex_amrcore_gridding (Precision prec);

    /** This builds one AMReX box per MPI-rank
     *
     * This is a simple decomposition for particles that we default to if
     * we do not need space charge or load balancing.
     *
     * @param prec the runtime-selected beam precision to construct the data at
     * @return simulation_geometry the geometry (topology) of the simulation;
     *         amr_info contains information on mesh-refinement and box/grid blocks
     */
    std::unique_ptr<AmrCoreDataBase>
    one_box_per_rank (Precision prec);

    //! the fake domain size (lower end), used to indicate that the user made no choice for it
    constexpr amrex::RealVect undefined_geometry_prob_lo(-1.0);
    //! the fake domain size (higher end), used to indicate that the user made no choice for it
    constexpr amrex::RealVect undefined_geometry_prob_hi(1.0);

} // namespace impactx::initialization

#endif // IMPACT_INIT_AMR_CORE_H
