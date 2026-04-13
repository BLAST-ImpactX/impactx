add_library(synergia_foundation four_momentum.cc reference_particle.cc
                                distribution.cc)
target_compile_definitions(synergia_foundation PUBLIC ${eigen_defs})
target_link_libraries(
  synergia_foundation
  synergia_parallel_utils
  lsexpr
  GSL::gsl
  GSL::gslcblas
  Eigen3::Eigen
  ${kokkos_libs})

target_link_options(synergia_foundation PRIVATE ${LINKER_OPTIONS})

if(BUILD_PYTHON_BINDINGS)
  pybind11_add_module(foundation MODULE NO_EXTRAS foundation_pywrap.cc)
  target_link_libraries(foundation PRIVATE synergia_foundation ${kokkos_libs})
  copy_file(__init__.py foundation)
endif()

install(TARGETS synergia_foundation DESTINATION ${SYNERGIA3_LIB_INSTALL_DIR})
install(FILES distribution.h four_momentum.h physical_constants.h
              reference_particle.h
        DESTINATION ${SYNERGIA3_INCLUDE_INSTALL_DIR}/synergia/foundation)

if(BUILD_PYTHON_BINDINGS)
  install(FILES __init__.py
          DESTINATION ${PYTHON_INSTALL_DIR}/synergia/foundation)
  install(TARGETS foundation
          DESTINATION ${PYTHON_INSTALL_DIR}/synergia/foundation)
endif()

add_subdirectory(tests)
