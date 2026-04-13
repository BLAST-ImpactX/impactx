# This library is built separately to avoid linking to kokkos directly and avoid
# calling host only functions from host-device code.
add_library(synergia_bunch_hostonly diagnostics_full2_host.cc populate_host.cc)
target_link_libraries(synergia_bunch_hostonly PUBLIC Eigen3::Eigen)
target_compile_definitions(synergia_bunch_hostonly PUBLIC ${eigen_defs})
target_link_options(synergia_bunch_hostonly PRIVATE ${LINKER_OPTIONS})

add_library(
  synergia_bunch
  bunch_particles.cc
  bunch.cc
  bunch_train.cc
  core_diagnostics.cc
  diagnostics_worker.cc
  diagnostics_io.cc
  diagnostics_full2.cc
  diagnostics_full2_host.cc
  diagnostics_particles.cc
  diagnostics_loss.cc
  diagnostics_bulk_track.cc
  populate.cc
  populate_global.cc
  populate_host.cc
  period.cc)
target_link_libraries(
  synergia_bunch
  PUBLIC synergia_bunch_hostonly
         synergia_foundation
         synergia_utils
         # $<$<STREQUAL:${USE_OPENPMD_IO},OFF>:synergia_hdf5_utils>
         synergia_hdf5_utils
         $<$<STREQUAL:${USE_OPENPMD_IO},ON>:openPMD::openPMD>
  PRIVATE ${kokkos_libs})
target_link_options(synergia_bunch PRIVATE ${LINKER_OPTIONS})

if(BUILD_PYTHON_BINDINGS)
  add_library(synergia_pydiag diagnostics_py.cc)

  # the -fvisibility flag suppress the warning, but it also prevents the pydiag
  # from exposing to Python as well set_target_properties(synergia_pydiag
  # PROPERTIES CXX_VISIBILITY_PRESET "hidden")

  target_link_libraries(synergia_pydiag PRIVATE synergia_bunch ${kokkos_libs}
                                                pybind11::module)

  pybind11_add_module(bunch MODULE NO_EXTRAS bunch_pywrap.cc)
  target_link_libraries(bunch PRIVATE synergia_bunch synergia_pydiag)

  copy_file(__init__.py bunch)

  install(FILES __init__.py DESTINATION ${PYTHON_INSTALL_DIR}/synergia/bunch)
  install(TARGETS synergia_pydiag DESTINATION ${SYNERGIA3_LIB_INSTALL_DIR})
  install(TARGETS bunch DESTINATION ${PYTHON_INSTALL_DIR}/synergia/bunch)

endif()

install(TARGETS synergia_bunch synergia_bunch_hostonly
        DESTINATION ${SYNERGIA3_LIB_INSTALL_DIR})

install(
  FILES bunch.h
        bunch_train.h
        core_diagnostics.h
        diagnostics.h
        diagnostics_full2.h
        diagnostics_track.h
        diagnostics_bulk_track.h
        diagnostics_particles.h
        fixed_t_z_converter.h
        populate.h
        period.h
  DESTINATION ${SYNERGIA3_INCLUDE_INSTALL_DIR}/synergia/bunch)

add_subdirectory(tests)
