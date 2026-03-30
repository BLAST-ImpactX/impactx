add_test(
  test_one_turn_map_py
  ${MPIEXEC}
  ${MPIEXEC_NUMPROC_FLAG}
  1
  ${SYNERGIA2_BINARY_DIR}/synergia-local
  -c
  "import nose; nose.main()"
  test_one_turn_map.py)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_one_turn_map.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/lattices
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
