add_executable(test_bunch_traits test_bunch_traits.cc ${test_main})
target_link_libraries(test_bunch_traits PRIVATE synergia_bunch ${testing_libs})
add_mpi_test(test_bunch_traits 1)

add_executable(test_bunch test_bunch.cc ${test_main})
target_link_libraries(test_bunch PRIVATE synergia_bunch ${testing_libs})
add_mpi_test(test_bunch 1)

add_executable(test_adjust_particles_ref_energy test_adjust_particles_ref_energy.cc ${test_main})
target_link_libraries(test_adjust_particles_ref_energy PRIVATE synergia_bunch ${testing_libs})
add_mpi_test(test_adjust_particles_ref_energy 1)

add_executable(test_offset_bunch_cdt test_offset_bunch_cdt.cc ${test_main})
target_link_libraries(test_offset_bunch_cdt PRIVATE synergia_bunch ${testing_libs})
add_mpi_test(test_offset_bunch_cdt 1)

add_executable(test_bunch_particles test_bunch_particles.cc ${test_main})
target_link_libraries(test_bunch_particles PRIVATE synergia_bunch
                                                   ${testing_libs})
add_mpi_test(test_bunch_particles 1)

if(BUILD_PYTHON_BINDINGS)
  add_py_test(test_bunch.py)
  add_py_test(test_bunch_reference_particles.py)
  add_py_test(test_adjust_particles_ref_energy.py)
  add_py_test(test_offset_bunch_cdt.py)
  add_py_test(test_particle_masks.py)
  endif()
