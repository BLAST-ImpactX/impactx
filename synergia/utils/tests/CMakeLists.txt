if(BUILD_PYTHON_BINDINGS)
  add_py_test(test_parallel_utils.py)
endif()

add_executable(test_command_line_arg test_command_line_arg.cc ${test_main})
target_link_libraries(test_command_line_arg PRIVATE synergia_command_line
                                                    ${testing_libs})
add_mpi_test(test_command_line_arg 1)

add_executable(test_commxx_mpi test_commxx_mpi.cc ${test_main})
target_link_libraries(test_commxx_mpi PRIVATE synergia_parallel_utils
                                              ${testing_libs})
add_mpi_test(test_commxx_mpi 1)
add_mpi_test(test_commxx_mpi 2)
add_mpi_test(test_commxx_mpi 3)
add_mpi_test(test_commxx_mpi 4)

add_executable(test_commxx_serdes test_commxx_serdes.cc ${test_main})
target_link_libraries(test_commxx_serdes PRIVATE synergia_parallel_utils
                                                 ${testing_libs})
add_mpi_test(test_commxx_serdes 1)
add_mpi_test(test_commxx_serdes 2)
add_mpi_test(test_commxx_serdes 3)
add_mpi_test(test_commxx_serdes 4)

add_executable(test_distributed_fft2d test_distributed_fft2d.cc ${test_main})
target_link_libraries(test_distributed_fft2d PRIVATE synergia_distributed_fft
                                                     ${testing_libs})
add_mpi_test(test_distributed_fft2d 1)
add_mpi_test(test_distributed_fft2d 2)
add_mpi_test(test_distributed_fft2d 3)
add_mpi_test(test_distributed_fft2d 4)

add_executable(test_distributed_fft3d test_distributed_fft3d.cc ${test_main})
target_link_libraries(test_distributed_fft3d PRIVATE synergia_distributed_fft
                                                     ${testing_libs})
add_mpi_test(test_distributed_fft3d 1)
add_mpi_test(test_distributed_fft3d 2)
add_mpi_test(test_distributed_fft3d 3)
add_mpi_test(test_distributed_fft3d 4)

if(NOT ${USE_OPENPMD_IO})
  add_executable(test_hdf5_file_mpi test_hdf5_file_mpi.cc ${test_main})
  target_link_libraries(test_hdf5_file_mpi PRIVATE synergia_hdf5_utils
                                                   ${testing_libs})

  add_mpi_test(test_hdf5_file_mpi 1)
  add_mpi_test(test_hdf5_file_mpi 2)
  add_mpi_test(test_hdf5_file_mpi 3)
  add_mpi_test(test_hdf5_file_mpi 4)

  add_executable(test_hdf5_write_mpi test_hdf5_write_mpi.cc ${test_main})
  target_link_libraries(test_hdf5_write_mpi PRIVATE synergia_hdf5_utils
                                                    ${testing_libs})

  add_mpi_test(test_hdf5_write_mpi 1)
  add_mpi_test(test_hdf5_write_mpi 2)
  add_mpi_test(test_hdf5_write_mpi 3)
  add_mpi_test(test_hdf5_write_mpi 4)

  add_executable(test_hdf5_append_mpi test_hdf5_append_mpi.cc ${test_main})
  target_link_libraries(test_hdf5_append_mpi PRIVATE synergia_hdf5_utils
                                                     ${testing_libs})

  add_mpi_test(test_hdf5_append_mpi 1)
  add_mpi_test(test_hdf5_append_mpi 2)
  add_mpi_test(test_hdf5_append_mpi 3)
  add_mpi_test(test_hdf5_append_mpi 4)

  add_executable(test_hdf5_read_mpi test_hdf5_read_mpi.cc ${test_main})
  target_link_libraries(test_hdf5_read_mpi PRIVATE synergia_hdf5_utils
                                                   ${testing_libs})

  add_mpi_test(test_hdf5_read_mpi 1)
  add_mpi_test(test_hdf5_read_mpi 2)
  add_mpi_test(test_hdf5_read_mpi 3)
  add_mpi_test(test_hdf5_read_mpi 4)
endif()

if("${ENABLE_KOKKOS_BACKEND}" STREQUAL "CUDA")
  add_executable(test_kokkos test_kokkos.cc ${test_main})
  target_link_libraries(test_kokkos PRIVATE ${testing_libs})

  add_mpi_test(test_kokkos 1)
endif()

if(0)
  add_executable(test_loading_kokkos test_loading_kokkos.cc)
  target_link_libraries(test_loading_kokkos PRIVATE dl)
  target_include_directories(
    test_loading_kokkos
    PRIVATE ${SYNERGIA2_BINARY_DIR}/src/synergia/utils/kokkos
            ${SYNERGIA2_BINARY_DIR}/src/synergia/utils/kokkos/core/src
            ${SYNERGIA2_BINARY_DIR}/src/synergia/utils/kokkos/core/src/desul
            ${SYNERGIA2_SOURCE_DIR}/src/synergia/utils/kokkos/tpls/desul/include
            ${SYNERGIA2_SOURCE_DIR}/src/synergia/utils/kokkos/core/src)
  add_test(test_loading_kokkos test_loading_kokkos)
endif()
