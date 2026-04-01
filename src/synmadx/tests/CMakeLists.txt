# Fetch Catch2 v3 for unit tests
include(FetchContent)
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.5.2
)
FetchContent_MakeAvailable(Catch2)

set(SYNMADX_TEST_SOURCES
    test_lattice.cc
    test_dynamic_lattice.cc
    test_madx_parser.cc
    test_mx_expr.cc
)

foreach(test_src ${SYNMADX_TEST_SOURCES})
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name} ${test_src})
    target_link_libraries(${test_name} PRIVATE synmadx Catch2::Catch2WithMain)
    add_test(NAME synmadx.${test_name} COMMAND ${test_name})
endforeach()

# Copy test data files
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/foo.dbx
               ${CMAKE_CURRENT_BINARY_DIR}/foo.dbx COPYONLY)
