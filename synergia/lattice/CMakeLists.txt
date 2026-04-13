# This library is built separately to avoid linking to kokkos directly and avoid
# calling host only functions from host-device code.
add_library(synergia_lattice_hostonly mx_expr.cc mx_parse.cc mx_tree.cc madx.cc
                                      lattice_tree.cc)
target_link_libraries(synergia_lattice_hostonly synergia_foundation
                      synergia_boost)
target_compile_definitions(synergia_lattice_hostonly PUBLIC ${eigen_defs})
target_link_options(synergia_lattice_hostonly PRIVATE ${LINKER_OPTIONS})

add_library(
  synergia_lattice
  dynamic_lattice.cc madx_reader.cc lattice_element.cc lattice_element_slice.cc
  lattice_element_processor.cc lattice.cc)
target_link_libraries(synergia_lattice PUBLIC synergia_lattice_hostonly
                                              synergia_hdf5_utils)
target_link_options(synergia_lattice PRIVATE ${LINKER_OPTIONS})

if(BUILD_PYTHON_BINDINGS)

  pybind11_add_module(lattice MODULE NO_EXTRAS lattice_pywrap.cc)
  target_link_libraries(lattice PRIVATE synergia_lattice)

  copy_file(__init__.py lattice)
  copy_file(mad8_parser.py lattice)
  copy_file(mad8_reader.py lattice)
  copy_file(simplify.py lattice)

endif()

install(TARGETS synergia_lattice synergia_lattice_hostonly
        DESTINATION ${SYNERGIA3_LIB_INSTALL_DIR})

install(
  FILES lattice_element.h
        lattice_element_slice.h
        lattice.h
        madx.h
        madx_reader.h
        mx_expr.h
        mx_parse.h
        mx_tree.h
  DESTINATION ${SYNERGIA3_INCLUDE_INSTALL_DIR}/synergia/lattice)

if(BUILD_PYTHON_BINDINGS)
  install(FILES __init__.py mad8_parser.py mad8_reader.py simplify.py
          DESTINATION ${PYTHON_INSTALL_DIR}/synergia/lattice)

  install(TARGETS lattice DESTINATION ${PYTHON_INSTALL_DIR}/synergia/lattice)

endif()

add_subdirectory(tests)
