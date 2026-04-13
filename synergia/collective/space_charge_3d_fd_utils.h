#ifndef SPACE_CHARGE_3D_FD_UTILS_H_
#define SPACE_CHARGE_3D_FD_UTILS_H_

#include "logger.h"
#include "space_charge_3d_fd_impl.h"

PetscErrorCode init_solver_subcomms(SubcommCtx& sctx, GlobalCtx& gctx);

PetscErrorCode init_local_vecs(LocalCtx& lctx, GlobalCtx& gctx);
PetscErrorCode init_subcomm_vecs(SubcommCtx& sctx, GlobalCtx& gctx);

PetscErrorCode init_subcomm_mat(LocalCtx& lctx,
                                SubcommCtx& sctx,
                                GlobalCtx& gctx);

PetscErrorCode init_solver(LocalCtx& lctx,
                           SubcommCtx& sctx,
                           GlobalCtx& gctx,
                           Logger& logger,
                           bool fixed_domain);

PetscErrorCode compute_mat(LocalCtx& lctx, SubcommCtx& sctx, GlobalCtx& gctx);

PetscErrorCode solve(SubcommCtx& sctx, GlobalCtx& gctx);

PetscErrorCode init_global_subcomm_scatters(SubcommCtx& sctx, GlobalCtx& gctx);
PetscErrorCode init_subcomm_local_scatters(LocalCtx& lctx,
                                           SubcommCtx& sctx,
                                           GlobalCtx& gctx);

PetscErrorCode finalize(LocalCtx& lctx, SubcommCtx& sctx, GlobalCtx& gctx);

PetscErrorCode MyMonitor(KSP ksp, PetscInt it, PetscReal rnorm, void* mctx);

#endif
