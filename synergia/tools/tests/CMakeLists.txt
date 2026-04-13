add_definitions(-DBOOST_TEST_DYN_LINK)

# add_test(NAME test_three_bump_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
# ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()"
# test_three_bump.py)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_three_bump.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/lattices
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
                                    "bump_basic.h5;step_basic.h5")
