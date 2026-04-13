# synergia components
add_subdirectory(utils)
add_subdirectory(foundation)
add_subdirectory(bunch)
add_subdirectory(lattice)
add_subdirectory(libFF)
add_subdirectory(simulation)
add_subdirectory(collective)
add_subdirectory(tools)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set(SYNERGIA_MPI_LIBRARY_VERSION "")
string(REGEX REPLACE "[\n\t]+" " " SYNERGIA_MPI_LIBRARY_VERSION
                     "${MPI_C_LIBRARY_VERSION_STRING}")

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/version.py.in"
               "${CMAKE_CURRENT_BINARY_DIR}/version.py" IMMEDIATE)

add_library(synergia INTERFACE)
target_link_libraries(
  synergia INTERFACE synergia_utils synergia_foundation synergia_bunch
                     synergia_lattice synergia_collective synergia_simulation)

if(BUILD_PYTHON_BINDINGS)
  # add_subdirectory(convertors)
  copy_file(__init__.py synergia)
endif()

install(TARGETS synergia DESTINATION ${SYNERGIA3_LIB_INSTALL_DIR})
install(FILES __init__.py "${CMAKE_CURRENT_BINARY_DIR}/version.py"
        DESTINATION ${PYTHON_INSTALL_DIR}/synergia)
