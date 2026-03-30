#ifndef SPACE_CHARGE_3D_FD_IMPL_H_
#define SPACE_CHARGE_3D_FD_IMPL_H_

#include <functional>
#include <vector>

#include <petscdmda.h>
#include <petscksp.h>
#include <petscmat.h>
#include <petscviewerhdf5.h>

#include "synergia/utils/kokkos_views.h"
#include "synergia/utils/simple_timer.h"

/*! Local (per MPI rank) context */
struct LocalCtx {

    Vec seqphi;  /*! local seq vector */
    Vec seqrho;  /*! local seq vector */
    DM da3d_seq; /*! sequential DMDA to manage grid and vecs */

    karray1d_dev seqphi_view; /*! kokkos view for seqphi */
    karray1d_dev seqrho_view; /*! kokkos view for seqrho */

    karray1d_dev enx; /*! kokkos view for electric field along x */
    karray1d_dev eny; /*! kokkos view for electric field along y */
    karray1d_dev enz; /*! kokkos view for electric field along z */

    karray1d_dev coo_v; /*! kokkos view for matrix entries in coo format */
};

/*! Subcomm context */
struct SubcommCtx {

    Vec phi_subcomm_local; /*! subcomm alias of the local vector on each MPI
                             rank */
    Vec rho_subcomm_local; /*! subcomm alias of the local vector on each MPI
                              rank */
    Vec phi_subcomm;       /*! vector on the subcomm */
    Vec rho_subcomm;       /*! vector on the subcomm */

    DM da;                         /* DMDA to manage grid and vecs */
    Mat A;                         /* discretization matrix */
    KSP ksp;                       /* krylov solver */
    PC pc;                         /* preconditioner */
    PetscBool reuse = PETSC_FALSE; /* state of ksp/pc re-use */

    VecScatter scat_subcomm_to_local; /*! VecScatter from subcomm vector to
                                        constituent local vectors */
    IS ix_scat_subcomms_to_local; /*! IndexSet for scatter from subcomm vector
                                    to constituent local vectors */
    IS iy_scat_subcomms_to_local; /*! IndexSet for scatter from subcomm vector
                                    to constituent local vectors */

    PetscSubcomm solverpsubcomm; /*! PETSc subcommunicator for this concurrent
                                   solver instance */
    MPI_Comm solversubcomm; /*! MPI subcommunicator for this concurrent solver
                              instance */

    PetscMPIInt solversubcommid; /*! a unique Int on each solver subcomm (max
                                   value of global ranks in subcomm) */
    PetscMPIInt solver_rank;     /*! solver subcomm MPI communicator rank */
    PetscMPIInt solver_size;     /*! solver subcomm MPI communicator size */
};

/*! Global (over all MPI ranks encompassing the bunch) context */
struct GlobalCtx {

    PetscInt nsubcomms = 1; /*! total number of subcomms */
    PetscInt nsize;         /*! the size of the problem */

    PetscInt nsize_x; /*! the size of the grid along x axis */
    PetscInt nsize_y; /*! the size of the grid along y axis */
    PetscInt nsize_z; /*! the size of the grid along z axis */

    MPI_Comm bunch_comm; /*! MPI communicator over which the bunch has been
                            defined */
    bool ksp_converged_reason =
        false;             /*! enable KSPConvergedReason to stdout */
    bool ksp_view = false; /*! enable KSPView to stdout */
    bool ksp_monitor_residual =
        false;          /*! enable KSPMonitorResidual to stdout */
    bool debug = false; /*! enable verbose outputs */
    bool dumps = false; /*! enable dumping states to HDF5 files */

    PetscMPIInt global_rank = -1; /*! global MPI communicator rank */
    PetscMPIInt global_size = -1; /*! global MPI communicator size */

    double Lx; /*! length along x */
    double Ly; /*! length along x */
    double Lz; /*! length along x */

    double Lx_ref = -1; /*! length along x, when preconditioner was created */
    double Ly_ref = -1; /*! length along y, when preconditioner was created */
    double Lz_ref = -1; /*! length along z, when preconditioner was created */

    Vec phi_global_local; /*! global alias of the local vector on each MPI rank
                           */
    Vec rho_global_local; /*! global alias of the local vector on each MPI rank
                           */

    std::vector<Vec> phi_global_subcomm; /*! array of global alias of the vector
                                           on each subcomm */
    std::vector<Vec> rho_global_subcomm; /*! array of global alias of the vector
                                           on each subcomm */

    std::vector<VecScatter>
        scat_glocal_to_subcomms; /*! array of VecScatter(s) from global alias of
                                   local vectors to all subcomm vectors */
    IS ix_scat_glocal_to_subcomms; /*! IndexSet for scatters from global alias
                                     of local vectors to all subcomm vectors */
    IS iy_scat_glocal_to_subcomms; /*! IndexSet for scatters from global alias
                                     of local vectors to all subcomm vectors */
    std::vector<PetscMPIInt>
        sids; /*! holds the solversubcommmid's from each solver */

    /*! type of matrix created */
    static constexpr MatType
#if defined KOKKOS_ENABLE_CUDA
        mattype = "aijcusparse";
#elif defined KOKKOS_ENABLE_OPENMP
        mattype = "aij";
#endif

    /*! type of all vectors on subcomms */
    static constexpr VecType
#if defined KOKKOS_ENABLE_CUDA
        vectype = "cuda";
#elif defined KOKKOS_ENABLE_OPENMP
        vectype = "standard";
#endif

    /*! function pointer to hold the appropriate function for creating vectors
      with array */
    std::function<PetscErrorCode(MPI_Comm,
                                 PetscInt,
                                 PetscInt,
                                 PetscInt,
                                 const PetscScalar*,
                                 Vec*)>
#if defined KOKKOS_ENABLE_CUDA
        VecCreate_type_WithArray = VecCreateMPICUDAWithArray;
#elif defined KOKKOS_ENABLE_OPENMP
        VecCreate_type_WithArray = VecCreateMPIWithArray;
#endif
};

#endif
