add_library(synergia_bunchsim bunch_simulator.cc independent_operation.cc
                              independent_operator.cc operation_extractor.cc)
target_link_libraries(synergia_bunchsim synergia_foundation synergia_bunch
                      synergia_lattice)
target_link_options(synergia_bunchsim PRIVATE ${LINKER_OPTIONS})

add_library(synergia_simulation_hostonly checkpoint_json.cc
                                         lattice_simulator_host.cc)
target_compile_definitions(synergia_simulation_hostonly PUBLIC ${eigen_defs})
target_link_libraries(synergia_simulation_hostonly Eigen3::Eigen ${kokkos_libs})
target_link_options(synergia_simulation_hostonly PRIVATE ${LINKER_OPTIONS})

add_library(
  synergia_simulation
  step.cc
  independent_stepper_elements.cc
  split_operator_stepper.cc
  split_operator_stepper_elements.cc
  propagator.cc
  lattice_simulator.cc
  checkpoint.cc)
target_link_libraries(
  synergia_simulation
  PUBLIC synergia_simulation_hostonly synergia_foundation synergia_bunch
         synergia_lattice synergia_bunchsim synergia_collective)
target_link_options(synergia_simulation PRIVATE ${LINKER_OPTIONS})

if(BUILD_PYTHON_BINDINGS)

  pybind11_add_module(simulation MODULE NO_EXTRAS simulation_pywrap.cc)
  set_target_properties(simulation PROPERTIES CXX_VISIBILITY_PRESET "default")
  target_link_libraries(simulation PRIVATE synergia_simulation synergia_bunch
                                           synergia_pydiag)

  copy_file(__init__.py simulation)

endif()

install(TARGETS synergia_simulation_hostonly synergia_simulation
                synergia_bunchsim DESTINATION ${SYNERGIA3_LIB_INSTALL_DIR})

install(
  FILES bunch_simulator.h
        independent_operation.h
        aperture_operation.h
        lattice_simulator.h
        operation_extractor.h
        operator.h
        propagate_actions.h
        propagator.h
        step.h
        stepper.h
        independent_stepper_elements.h
        split_operator_stepper.h
        split_operator_stepper_elements.h
        populate_stationary.h
  DESTINATION ${SYNERGIA3_INCLUDE_INSTALL_DIR}/synergia/simulation)

if(BUILD_PYTHON_BINDINGS)
  install(FILES __init__.py
          DESTINATION ${PYTHON_INSTALL_DIR}/synergia/simulation)
  install(TARGETS simulation
          DESTINATION ${PYTHON_INSTALL_DIR}/synergia/simulation)
endif()

add_subdirectory(tests)
