#include "space_charge_3d_fd_alias.h"

/* --------------------------------------------------------------------- */
/*!
  Initialize global aliases of local vectors
  \param   lctx - local context
  \param   gctx - global context
  \return  ierr - PetscErrorCode
  */
PetscErrorCode
init_global_local_aliases(LocalCtx& lctx, GlobalCtx& gctx)
{

    PetscInt phi_localsize; /* size of vector on this MPI rank */
    PetscInt rho_localsize; /* size of vector on this MPI rank */

    PetscScalar const*
        d_rho_val; /* read-only pointer to vector's contents on device */
    PetscScalar const*
        d_phi_val; /* read-only pointer to vector's contents on device */

    PetscMemType mtype_phi; /* memory type for phi vector */
    PetscMemType mtype_rho; /* memory type for rho vector */

    PetscFunctionBeginUser;

    /* Get size of local vector */
    PetscCall(VecGetLocalSize(lctx.seqphi, &phi_localsize));
    PetscCall(VecGetLocalSize(lctx.seqrho, &rho_localsize));

    /* Get access to local vector */
    PetscCall(VecGetArrayReadAndMemType(lctx.seqphi, &d_phi_val, &mtype_phi));
    PetscCall(VecGetArrayReadAndMemType(lctx.seqrho, &d_rho_val, &mtype_rho));

    /* consistency check */
    if (gctx.debug) {
        if (mtype_phi != mtype_rho)
            SETERRQ(gctx.bunch_comm,
                    PETSC_ERR_ARG_NOTSAMETYPE,
                    "phi and rho vectors don't have the memory type!");
        if (!gctx.VecCreate_type_WithArray)
            SETERRQ(gctx.bunch_comm,
                    PETSC_ERR_ARG_BADPTR,
                    "the function pointer in gctx is invalid!");
    }

    PetscCall(gctx.VecCreate_type_WithArray(gctx.bunch_comm,
                                            1,
                                            phi_localsize,
                                            PETSC_DECIDE,
                                            d_phi_val,
                                            &gctx.phi_global_local));
    PetscCall(gctx.VecCreate_type_WithArray(gctx.bunch_comm,
                                            1,
                                            rho_localsize,
                                            PETSC_DECIDE,
                                            d_rho_val,
                                            &gctx.rho_global_local));

    PetscCall(PetscObjectSetName((PetscObject)(gctx.phi_global_local),
                                 "phi_global_local_on_gctx"));
    PetscCall(PetscObjectSetName((PetscObject)(gctx.rho_global_local),
                                 "rho_global_local_on_gctx"));

    /* Restore local vector arrays */
    PetscCall(VecRestoreArrayReadAndMemType(lctx.seqphi, &d_phi_val));
    PetscCall(VecRestoreArrayReadAndMemType(lctx.seqrho, &d_rho_val));

    PetscFunctionReturn(0);
}

/* --------------------------------------------------------------------- */
/*!
  Initialize global aliases of subcomm vectors
  \param   sctx - subcomm context
  \param   gctx - global context
  \return  ierr - PetscErrorCode
  */
PetscErrorCode
init_global_subcomm_aliases(SubcommCtx& sctx, GlobalCtx& gctx)
{

    PetscInt phi_localsize; /* size of vector on this MPI rank */
    PetscInt rho_localsize; /* size of vector on this MPI rank */
    PetscScalar const*
        d_phi_val; /* read-only pointer to vector's contents on device */
    PetscScalar const*
        d_rho_val; /* read-only pointer to vector's contents on device */
    PetscMemType mtype_phi;
    PetscMemType mtype_rho;
    PetscInt size;

    PetscFunctionBeginUser;

    /* Get size of local vector */
    PetscCall(VecGetLocalSize(sctx.phi_subcomm, &phi_localsize));
    PetscCall(VecGetLocalSize(sctx.rho_subcomm, &rho_localsize));

    /* Get access to local vector */
    PetscCall(
        VecGetArrayReadAndMemType(sctx.phi_subcomm, &d_phi_val, &mtype_phi));
    PetscCall(
        VecGetArrayReadAndMemType(sctx.rho_subcomm, &d_rho_val, &mtype_rho));

    /* resize to ensure we only create as many aliases
       as the number of subcomms ! */
    gctx.phi_global_subcomm.resize(gctx.nsubcomms);
    gctx.rho_global_subcomm.resize(gctx.nsubcomms);

    PetscCall(PetscBarrier(NULL));
    for (PetscInt i = 0; i < gctx.nsubcomms; i++) {

        if (gctx.sids[i] == sctx.solversubcommid) {

            if (gctx.debug) {
                PetscCall(PetscPrintf(gctx.bunch_comm, "\n"));
                PetscCall(
                    PetscPrintf(PetscObjectComm((PetscObject)sctx.phi_subcomm),
                                "Hi there from solver-subcomm number %d, with "
                                "subcommid of %d\n",
                                i,
                                gctx.sids[i]));
                PetscCall(PetscPrintf(gctx.bunch_comm, "\n"));
                if (mtype_phi != mtype_rho)
                    SETERRQ(gctx.bunch_comm,
                            PETSC_ERR_ARG_NOTSAMETYPE,
                            "phi and rho vectors don't have the memory type!");
                if (!gctx.VecCreate_type_WithArray)
                    SETERRQ(gctx.bunch_comm,
                            PETSC_ERR_ARG_BADPTR,
                            "the function pointer in gctx is invalid!");
            }

            /* Create global aliases of subcomm vectors */
            PetscCall(
                gctx.VecCreate_type_WithArray(gctx.bunch_comm,
                                              1,
                                              phi_localsize,
                                              PETSC_DECIDE,
                                              d_phi_val,
                                              &(gctx.phi_global_subcomm[i])));
            PetscCall(
                gctx.VecCreate_type_WithArray(gctx.bunch_comm,
                                              1,
                                              rho_localsize,
                                              PETSC_DECIDE,
                                              d_rho_val,
                                              &(gctx.rho_global_subcomm[i])));

        } else {

            /* ranks outside the subcomm whose vector is being aliased
               do not contribute any values to the global alias,
               yet they must participate in collective calls */
            PetscCall(
                gctx.VecCreate_type_WithArray(gctx.bunch_comm,
                                              1,
                                              0,
                                              PETSC_DECIDE,
                                              NULL,
                                              &(gctx.phi_global_subcomm[i])));
            PetscCall(
                gctx.VecCreate_type_WithArray(gctx.bunch_comm,
                                              1,
                                              0,
                                              PETSC_DECIDE,
                                              NULL,
                                              &(gctx.rho_global_subcomm[i])));
        }

        PetscCall(PetscObjectSetName(
            (PetscObject)(gctx.phi_global_subcomm[i]),
            std::string("phi_global_subcomm_" + std::to_string(i) + "_on_gctx")
                .c_str()));
        PetscCall(PetscObjectSetName(
            (PetscObject)(gctx.rho_global_subcomm[i]),
            std::string("rho_global_subcomm_" + std::to_string(i) + "_on_gctx")
                .c_str()));
    }
    PetscCall(PetscBarrier(NULL));

    /* Restore local vector arrays */
    PetscCall(VecRestoreArrayReadAndMemType(sctx.phi_subcomm, &d_phi_val));
    PetscCall(VecRestoreArrayReadAndMemType(sctx.rho_subcomm, &d_rho_val));

    if (gctx.debug) {
        for (PetscInt i = 0; i < gctx.nsubcomms; i++) {
            PetscCall(VecGetSize(gctx.rho_global_subcomm[i], &size));
            PetscCall(PetscPrintf(
                gctx.bunch_comm,
                "Hello! the size of the global alias of subcomm vec "
                "with index %d is %d\n",
                i,
                size));
        }
        PetscCall(PetscPrintf(gctx.bunch_comm, "\n"));
    }

    PetscFunctionReturn(0);
}

/* --------------------------------------------------------------------- */
/*!
  Initialize subcomm aliases of local vectors
  \param  lctx - local context
  \param  sctx - subcomm context
  \param  gctx - global context
  \return ierr - PetscErrorCode
  */
PetscErrorCode
init_subcomm_local_aliases(LocalCtx& lctx, SubcommCtx& sctx, GlobalCtx& gctx)
{

    PetscInt phi_localsize; /* size of vector on this MPI rank */
    PetscInt rho_localsize; /* size of vector on this MPI rank */
    PetscScalar const*
        d_phi_val; /* read-only pointer to vector's contents on device */
    PetscScalar const*
        d_rho_val; /* read-only pointer to vector's contents on device */
    PetscMemType mtype_phi;
    PetscMemType mtype_rho;

    PetscFunctionBeginUser;

    /* Get size of local vector */
    PetscCall(VecGetLocalSize(lctx.seqphi, &phi_localsize));
    PetscCall(VecGetLocalSize(lctx.seqrho, &rho_localsize));

    /* Get access to local vector */
    PetscCall(VecGetArrayReadAndMemType(lctx.seqphi, &d_phi_val, &mtype_phi));
    PetscCall(VecGetArrayReadAndMemType(lctx.seqrho, &d_rho_val, &mtype_rho));

    /* consistency check */
    if (gctx.debug) {
        if (mtype_phi != mtype_rho)
            SETERRQ(gctx.bunch_comm,
                    PETSC_ERR_ARG_NOTSAMETYPE,
                    "phi and rho vectors don't have the memory type!");
        if (!gctx.VecCreate_type_WithArray)
            SETERRQ(gctx.bunch_comm,
                    PETSC_ERR_ARG_BADPTR,
                    "the function pointer in gctx is invalid!");
    }

    /* Create subcomm aliases of local vectors */
    PetscCall(gctx.VecCreate_type_WithArray(
        PetscObjectComm(PetscObject(sctx.phi_subcomm)),
        1,
        phi_localsize,
        PETSC_DECIDE,
        d_phi_val,
        &sctx.phi_subcomm_local));
    PetscCall(gctx.VecCreate_type_WithArray(
        PetscObjectComm(PetscObject(sctx.phi_subcomm)),
        1,
        rho_localsize,
        PETSC_DECIDE,
        d_rho_val,
        &sctx.rho_subcomm_local));
    PetscCall(PetscObjectSetName((PetscObject)(sctx.phi_subcomm_local),
                                 "phi_subcomm_local_on_sctx"));
    PetscCall(PetscObjectSetName((PetscObject)(sctx.rho_subcomm_local),
                                 "rho_subcomm_local_on_sctx"));

    /* Restore local vector arrays */
    PetscCall(VecRestoreArrayReadAndMemType(lctx.seqphi, &d_phi_val));
    PetscCall(VecRestoreArrayReadAndMemType(lctx.seqrho, &d_rho_val));

    PetscFunctionReturn(0);
}
