# synmadx: standalone MadX lattice parser imported from Synergia
# No external dependencies beyond C++ standard library and bundled Boost headers

add_library(synmadx SHARED
    synergia/lattice/mx_expr.cc
    synergia/lattice/mx_parse.cc
    synergia/lattice/mx_tree.cc
    synergia/lattice/madx.cc
    synergia/lattice/lattice_tree.cc
    synergia/lattice/madx_reader.cc
    synergia/lattice/lattice.cc
    synergia/lattice/lattice_element.cc
    synergia/lattice/lattice_element_slice.cc
    synergia/lattice/lattice_element_processor.cc
    synergia/lattice/dynamic_lattice.cc
    synergia/foundation/four_momentum.cc
    synergia/foundation/reference_particle.cc
    synergia/utils/lsexpr.cc
)

set_target_properties(synmadx PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(synmadx PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/synergia/utils
)

target_compile_features(synmadx PUBLIC cxx_std_17)

# install shared library
install(TARGETS synmadx
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# Python bindings
if(ImpactX_PYTHON)
    pybind11_add_module(pysynmadx MODULE python/synmadx_pywrap.cc)
    target_link_libraries(pysynmadx PRIVATE synmadx)
    set_target_properties(pysynmadx PROPERTIES
        OUTPUT_NAME synmadx_pybind
        CXX_VISIBILITY_PRESET "hidden"
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_PYTHON_OUTPUT_DIRECTORY}/impactx/synmadx
    )
    pybind11_extension(pysynmadx)
    if(NOT MSVC AND NOT ${CMAKE_BUILD_TYPE} MATCHES Debug|RelWithDebInfo)
        pybind11_strip(pysynmadx)
    endif()
endif()

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
