# add pytest tests
#
set(pytest_name pytest.ImpactX)
set(pytest_rundir ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${pytest_name})

#   make a unique run directory
file(MAKE_DIRECTORY ${pytest_rundir})
file(MAKE_DIRECTORY ${pytest_rundir}.dashboard)

#   copy input files
file(COPY ${ImpactX_SOURCE_DIR}/examples
     DESTINATION ${pytest_rundir})

#   run
add_test(NAME ${pytest_name}
    COMMAND ${Python_EXECUTABLE} -m pytest -s -vvvv
        ${ImpactX_SOURCE_DIR}/tests/python
        --ignore=${ImpactX_SOURCE_DIR}/tests/python/dashboard
    WORKING_DIRECTORY ${pytest_rundir}
)

#    run dashboard-specific test
set(dashboard_test_name ${pytest_name}.dashboard)
add_test(NAME ${dashboard_test_name}
    COMMAND ${Python_EXECUTABLE} -m pytest -s -vvvv
        ${ImpactX_SOURCE_DIR}/tests/python/dashboard
    WORKING_DIRECTORY ${pytest_rundir}.dashboard
)
set_property(TEST ${dashboard_test_name} PROPERTY LABELS "dashboard")

#   limit threads
set_property(TEST ${pytest_name} APPEND PROPERTY ENVIRONMENT "OMP_NUM_THREADS=2")
set_property(TEST ${pytest_name}.dashboard APPEND PROPERTY ENVIRONMENT "OMP_NUM_THREADS=2")


#   set PYTHONPATH and PATH (for .dll files)
impactx_test_set_pythonpath(${pytest_name})
impactx_test_set_pythonpath(${pytest_name}.dashboard)
