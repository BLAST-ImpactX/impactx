add_executable(test_libff_elements test_libff_elements.cc ${test_main})
target_link_libraries(test_libff_elements PRIVATE synergia_simulation
                                                  ${testing_libs})
add_mpi_test(test_libff_elements 1)
copy_file(fodo.madx test_libff_elements)

add_executable(test_madx_multipoles test_madx_multipoles.cc ${test_main})
target_link_libraries(test_madx_multipoles PRIVATE synergia_simulation
                                                   ${testing_libs})
add_mpi_test(test_madx_multipoles 1)

add_executable(test_madx_elements test_madx_elements.cc ${test_main})
target_link_libraries(test_madx_elements PRIVATE synergia_simulation
                                                 ${testing_libs})
add_mpi_test(test_madx_elements 1)

add_executable(test_foils test_foils.cc ${test_main})

target_link_libraries(test_foils PRIVATE synergia_simulation ${testing_libs})
add_mpi_test(test_foils 1)

add_subdirectory(lattices)

add_executable(test_zeroparticle_propagate test_zeroparticle_propagate.cc
                                           ${test_main})

target_link_libraries(test_zeroparticle_propagate PRIVATE synergia_simulation
                                                          ${testing_libs})

add_mpi_test(test_zeroparticle_propagate 1)

add_executable(test_zeroparticle2_propagate test_zeroparticle2_propagate.cc
                                            ${test_main})

target_link_libraries(test_zeroparticle2_propagate PRIVATE synergia_simulation
                                                           ${testing_libs})

add_mpi_test(test_zeroparticle2_propagate 1)

add_executable(test_booster_propagate test_booster_propagate.cc ${test_main})

target_link_libraries(test_booster_propagate PRIVATE synergia_simulation
                                                     ${testing_libs})

add_mpi_test(test_booster_propagate 1)

add_executable(test_booster2_propagate test_booster2_propagate.cc ${test_main})

target_link_libraries(test_booster2_propagate PRIVATE synergia_simulation
                                                      ${testing_libs})
add_mpi_test(test_booster2_propagate 1)
