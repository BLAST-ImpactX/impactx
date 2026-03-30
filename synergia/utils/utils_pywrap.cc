
#include <Kokkos_Core.hpp>
#include <pybind11/pybind11.h>

#if defined BUILD_FD_SPACE_CHARGE_SOLVER
#include <petsc.h>
#include <string>
#endif

namespace py = pybind11;
using namespace py::literals;

PYBIND11_MODULE(utils, m)
{
  m.def("init", []() {

#if defined BUILD_FD_SPACE_CHARGE_SOLVER
    PetscErrorCode ierr;
    ierr = PetscInitialize(nullptr,
                           nullptr,
                           (char*)0,
                           std::string("synergia2-v3 program!\n").c_str());
#endif
    auto settings =
      Kokkos::InitializationSettings(); /* use default constructor */
    auto num_threads_chars = std::getenv("OMP_NUM_THREADS");
    if (num_threads_chars != nullptr) {
      settings.set_num_threads(std::stoi(num_threads_chars));
    }
    Kokkos::initialize(settings);
  });

  m.def("finalize", []() {
    Kokkos::finalize();
#if defined BUILD_FD_SPACE_CHARGE_SOLVER
    PetscErrorCode ierr;
    ierr = PetscFinalize();
#endif
  });
}
