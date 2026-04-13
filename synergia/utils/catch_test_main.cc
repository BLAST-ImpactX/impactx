

#include <Kokkos_Core.hpp>
#include <catch2/catch_session.hpp>
#include <mpi.h>

#if defined BUILD_FD_SPACE_CHARGE_SOLVER
#include <petscsys.h>
static std::string help = "This is a synergia2 test program!\n\n";
#endif

int
main(int argc, char* argv[])
{

#if defined BUILD_FD_SPACE_CHARGE_SOLVER
    PETSC_MPI_THREAD_REQUIRED = MPI_THREAD_MULTIPLE;
    [[maybe_unused]] PetscErrorCode ierr;
    ierr = PetscInitialize(&argc, &argv, (char*)0, help.c_str());
    auto initargs =
        Kokkos::InitializationSettings{}; /* use default constructor */
    Kokkos::initialize(initargs);
#else
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);
#endif

    int result = Catch::Session().run(argc, argv);

    Kokkos::finalize();
#if defined BUILD_FD_SPACE_CHARGE_SOLVER
    ierr = PetscFinalize();
#else
    MPI_Finalize();
#endif

    return result;
}
