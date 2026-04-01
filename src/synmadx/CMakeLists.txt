# synmadx: standalone MadX lattice parser imported from Synergia
# No external dependencies beyond C++ standard library and bundled Boost headers

add_library(synmadx
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

target_include_directories(synmadx PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/synergia/utils
)

target_compile_features(synmadx PUBLIC cxx_std_17)

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
