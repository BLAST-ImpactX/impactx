if(BUILD_PYTHON_BINDINGS)
  add_py_test(test_lattice.py)
  add_py_test(test_lattice_element.py)
  add_py_test(test_dynamic_lattice.py)
  add_py_test(test_parse_matrix.py)
endif()

add_executable(test_lattice test_lattice.cc ${test_main})
target_link_libraries(test_lattice PRIVATE synergia_lattice ${testing_libs})
add_mpi_test(test_lattice 1)

add_executable(test_mx_expr test_mx_expr.cc ${test_main})
target_link_libraries(test_mx_expr PRIVATE synergia_lattice ${testing_libs})
add_mpi_test(test_mx_expr 1)

add_executable(test_madx_parser test_madx_parser.cc ${test_main})
target_link_libraries(test_madx_parser PRIVATE synergia_lattice ${testing_libs})
add_mpi_test(test_madx_parser 1)

copy_file(foo.dbx test_madx_parser)

add_executable(test_dynamic_lattice test_dynamic_lattice.cc ${test_main})
target_link_libraries(test_dynamic_lattice PRIVATE synergia_lattice
                                                   ${testing_libs})
add_mpi_test(test_dynamic_lattice 1)
