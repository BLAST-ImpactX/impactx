if(BUILD_PYTHON_BINDINGS)
  add_py_test(test_verify_mpi4py.py)
  add_test(
    NAME test_mpi4py_works
    COMMAND
      ${CMAKE_COMMAND} -E env OMP_NUM_THREADS=1 OMP_PROC_BIN=spread
      OMP_PLACES=cores OMPI_MCA_rmaps_base_oversubscribe=true
      PRTE_MCA_rmaps_default_mapping_policy=:oversubscribe
      ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
      ${Python_EXECUTABLE} -m mpi4py.bench helloworld ${MPIEXEC_POSTFLAGS})
  add_py_test(test_version.py)
  add_py_test(test_four_momentum.py)
  add_py_test(test_reference_particle.py)
  add_py_test(test_kg_to_GeV.py)

endif()

add_executable(test_trigon test_trigon.cc ${test_main})
target_link_libraries(test_trigon PRIVATE synergia_foundation ${testing_libs})
add_mpi_test(test_trigon 1)

add_executable(test_four_momentum test_four_momentum.cc ${test_main})
target_link_libraries(test_four_momentum PRIVATE synergia_foundation
                                                 ${testing_libs})
add_mpi_test(test_four_momentum 1)

add_executable(test_reference_particle test_reference_particle.cc ${test_main})
target_link_libraries(test_reference_particle PRIVATE synergia_foundation
                                                      ${testing_libs})
add_mpi_test(test_reference_particle 1)
