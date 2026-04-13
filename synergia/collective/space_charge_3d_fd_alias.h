#ifndef SPACE_CHARGE_3D_FD_ALIAS_H_
#define SPACE_CHARGE_3D_FD_ALIAS_H_

#include "space_charge_3d_fd_impl.h"

PetscErrorCode init_global_local_aliases(LocalCtx& lctx, GlobalCtx& gctx);

PetscErrorCode init_global_subcomm_aliases(SubcommCtx& sctx, GlobalCtx& gctx);

PetscErrorCode init_subcomm_local_aliases(LocalCtx& lctx,
                                          SubcommCtx& sctx,
                                          GlobalCtx& gctx);

#endif
