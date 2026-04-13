add_executable(test_deposit test_deposit.cc ${test_main})
target_link_libraries(test_deposit PRIVATE synergia_collective ${testing_libs})
add_mpi_test(test_deposit 1)

add_executable(test_space_charge_3d_open_hockney_mpi
               test_space_charge_3d_open_hockney_mpi.cc ${test_main})
target_link_libraries(test_space_charge_3d_open_hockney_mpi
                      PRIVATE synergia_collective ${testing_libs})
add_mpi_test(test_space_charge_3d_open_hockney_mpi 1)

add_executable(test_space_charge_3d_rectangular_mpi
               test_space_charge_3d_rectangular_mpi.cc ${test_main})
target_link_libraries(test_space_charge_3d_rectangular_mpi
                      PRIVATE synergia_collective ${testing_libs})
add_mpi_test(test_space_charge_3d_rectangular_mpi 1)

if(BUILD_PYTHON_BINDINGS)
  add_py_test(test_impedance_options.py)
endif()

if(${BUILD_FD_SPACE_CHARGE_SOLVER})
  add_executable(test_space_charge_3d_fd_mpi test_space_charge_3d_fd_mpi.cc
                                             ${test_main})
  target_link_libraries(test_space_charge_3d_fd_mpi PRIVATE synergia_collective
                                                            ${testing_libs})
  add_mpi_test(test_space_charge_3d_fd_mpi 1)

endif()
