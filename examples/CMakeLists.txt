# Configuration ###############################################################
#
if(ImpactX_MPI)
    # OpenMPI root guard: https://github.com/open-mpi/ompi/issues/4451
    if("$ENV{USER}" STREQUAL "root")
        # calling even --help as root will abort and warn on stderr
        execute_process(COMMAND ${MPIEXEC_EXECUTABLE} --help
                ERROR_VARIABLE MPIEXEC_HELP_TEXT
                OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(${MPIEXEC_HELP_TEXT} MATCHES "^.*allow-run-as-root.*$")
            set(MPI_ALLOW_ROOT --allow-run-as-root)
        endif()
    endif()
    set(MPI_TEST_EXE
        ${MPIEXEC_EXECUTABLE}
        ${MPI_ALLOW_ROOT}
        ${MPIEXEC_NUMPROC_FLAG} 2
        ${MPIEXEC_POSTFLAGS}
        ${MPIEXEC_PREFLAGS}
    )
endif()


# Add an ImpactX test set (with sub-tests).
#
# The test will be run and dependent follow-up tests will be added. The latter
# can optionally perform an analysis and/or a plot test on the generated data.
#
# Parameters
# ----------
# name: unique ImpactX test name
# input: inputs file or Python script
# is_mpi: run with mpirun?
# analysis_script: an analysis script that validates test output
# plot_script: a plot script that visualizes test output
#
function(add_impactx_test name input is_mpi analysis_script plot_script)
    # is this a Python script (ends in .py)
    if(input MATCHES ".*\.py$")
        set(is_python ON)
    else()
        set(is_python OFF)
    endif()

    # cannot run Python tests w/o Python support
    if(is_python AND NOT ImpactX_PYTHON)
        return()
    endif()

    # cannot run MPI tests w/o MPI support
    if(is_mpi AND NOT ImpactX_MPI)
        return()
    endif()

    # make a unique run directory
    file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${name})

    # get input file/script and optional command-line arguments
    separate_arguments(INPUTS_LIST UNIX_COMMAND "${input}")
    list(GET INPUTS_LIST 0 INPUTS_FILE)
    list(LENGTH INPUTS_LIST INPUTS_LIST_LENGTH)
    if(INPUTS_LIST_LENGTH GREATER 1)
        list(SUBLIST INPUTS_LIST 1 -1 INPUTS_ARGS)
        list(JOIN INPUTS_ARGS " " INPUTS_ARGS)
    else()
        set(INPUTS_ARGS "")
    endif()
    cmake_path(SET INPUTS_FILE "${ImpactX_SOURCE_DIR}/${INPUTS_FILE}")

    # get analysis script and optional command-line arguments
    separate_arguments(ANALYSIS_LIST UNIX_COMMAND "${analysis_script}")
    list(GET ANALYSIS_LIST 0 ANALYSIS_FILE)
    cmake_path(SET ANALYSIS_FILE "${CMAKE_CURRENT_SOURCE_DIR}/${ANALYSIS_FILE}")
    list(LENGTH ANALYSIS_LIST ANALYSIS_LIST_LENGTH)
    if(ANALYSIS_LIST_LENGTH GREATER 1)
        list(SUBLIST ANALYSIS_LIST 1 -1 ANALYSIS_ARGS)
        list(JOIN ANALYSIS_ARGS " " ANALYSIS_ARGS)
    else()
        set(ANALYSIS_ARGS "")
    endif()

    # test run
    set(THIS_WORKING_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${name})
    set(THIS_MPI_TEST_EXE)
    if(is_mpi)
        set(THIS_MPI_TEST_EXE ${MPI_TEST_EXE})
    endif()
    set(THIS_Python_EXE)
    if(is_python)
        set(THIS_Python_EXE ${Python_EXECUTABLE})
        # for argparse, do not pass command-line arguments as one quoted string
        separate_arguments(INPUTS_ARGS UNIX_COMMAND "${INPUTS_ARGS}")
        add_test(NAME ${name}.run
                 COMMAND ${THIS_MPI_TEST_EXE} ${THIS_Python_EXE} ${INPUTS_FILE} ${INPUTS_ARGS}
                 WORKING_DIRECTORY ${THIS_WORKING_DIR}
        )
        # TODO:
        # amrex.throw_exception=1
        # amrex.signal_handling=0
        # amrex.the_arena_init_size=0
        impactx_test_set_pythonpath(${name}.run)
    else()
        add_test(NAME ${name}.run
                 COMMAND
                    ${THIS_MPI_TEST_EXE} $<TARGET_FILE:app> ${INPUTS_FILE}
                        amrex.abort_on_unused_inputs=1
                        amrex.throw_exception=1
                        amrex.signal_handling=0
                        # allocate GPU memory on-demand instead of pre-allocating 3/4th
                        # to enable parallel test runs on the same GPU
                        amrex.the_arena_init_size=0
                        impactx.always_warn_immediately=1
                        impactx.abort_on_warning_threshold=low
                        ${INPUTS_ARGS}
                 WORKING_DIRECTORY ${THIS_WORKING_DIR}
        )
    endif()
    if(is_mpi)
        set_property(TEST ${name}.run APPEND PROPERTY ENVIRONMENT "OMP_NUM_THREADS=1")
    else()
        set_property(TEST ${name}.run APPEND PROPERTY ENVIRONMENT "OMP_NUM_THREADS=2")
    endif()
    # special return code for skipped tests (e.g., runtime prerequisite fails)
    set_tests_properties(${name}.run PROPERTIES SKIP_RETURN_CODE 42)

    # analysis and plots
    set(THIS_Python_SCRIPT_EXE)
    if(is_python OR WIN32)
        set(THIS_Python_SCRIPT_EXE ${Python_EXECUTABLE})
    endif()
    if(analysis_script)
        add_test(NAME ${name}.analysis
                 COMMAND ${THIS_Python_SCRIPT_EXE} ${ImpactX_SOURCE_DIR}/${analysis_script}
                 WORKING_DIRECTORY ${THIS_WORKING_DIR}
        )
        set_property(TEST ${name}.analysis APPEND PROPERTY DEPENDS "${name}.run")

        # make HDF5 I/O more robust on various filesystems
        set_property(TEST ${name}.analysis APPEND PROPERTY ENVIRONMENT "HDF5_USE_FILE_LOCKING=FALSE")

        # run test failed? Mark this as skipped
        set_property(TEST ${name}.analysis PROPERTY SKIP_REGULAR_EXPRESSION
            "Supplied directory is not valid: diags"
        )
    endif()
    if(plot_script)
        add_test(NAME ${name}.plot
                 COMMAND ${THIS_Python_SCRIPT_EXE} ${ImpactX_SOURCE_DIR}/${plot_script} --save-png
                 WORKING_DIRECTORY ${THIS_WORKING_DIR}
        )
        set_property(TEST ${name}.plot APPEND PROPERTY DEPENDS "${name}.run")

        # make HDF5 I/O more robust on various filesystems
        set_property(TEST ${name}.plot APPEND PROPERTY ENVIRONMENT "HDF5_USE_FILE_LOCKING=FALSE")

        # run test failed? Mark this as skipped
        set_property(TEST ${name}.plot PROPERTY SKIP_REGULAR_EXPRESSION
            "ValueError: No objects to concatenate"
        )
    endif()

    # CI: remove test directory after run
    if(ImpactX_TEST_CLEANUP)
        add_test(
            NAME ${name}.cleanup
            COMMAND ${CMAKE_COMMAND} -P ${ImpactX_SOURCE_DIR}/examples/test_cleanup.cmake ${THIS_WORKING_DIR}
        )
        # test cleanup depends on test run
        set_property(TEST ${name}.cleanup APPEND PROPERTY DEPENDS "${name}.run")
        if(analysis_script)
            # test cleanup depends on test analysis
            set_property(TEST ${name}.cleanup APPEND PROPERTY DEPENDS "${name}.analysis")
        endif()
        if(plot_script)
            # test cleanup depends on test analysis
            set_property(TEST ${name}.cleanup APPEND PROPERTY DEPENDS "${name}.plot")
        endif()
    endif()
endfunction()


# Add a CTest label to an ImpactX test set.
#
# The test set must already been added with add_impactx_test. Labeling it here
# will add the label to the run test and its optional (dependent) analysis
# and/or plot follow-up tests.
#
# Parameters
# ----------
# name: unique ImpactX test name
# label: ctest LABELS property value to add
#
function(label_impactx_test name label)
    # see add_test calls in add_impactx_test for values
    set(_test_names "${name}.run;${name}.analysis;${name}.plot")

    foreach(_test_name IN LISTS _test_names)
        if(TEST ${_test_name})
            set_property(TEST ${_test_name} APPEND PROPERTY LABELS "${label}")
        endif()
    endforeach()
endfunction()


# FODO Cell ###################################################################
#
add_impactx_test(FODO
    examples/fodo/input_fodo.in
    OFF  # ImpactX MPI-parallel
    examples/fodo/analysis_fodo.py
    examples/fodo/plot_fodo.py
)
add_impactx_test(FODO.py
    examples/fodo/run_fodo.py
    OFF  # ImpactX MPI-parallel
    examples/fodo/analysis_fodo.py
    examples/fodo/plot_fodo.py
)


# FODO Cell w/ Twiss input ####################################################
#
add_impactx_test(FODO.Twiss
    examples/fodo/input_fodo_twiss.in
    OFF  # ImpactX MPI-parallel
    examples/fodo/analysis_fodo.py
    examples/fodo/plot_fodo.py
)
add_impactx_test(FODO.Twiss.py
    examples/fodo/run_fodo_twiss.py
    OFF  # ImpactX MPI-parallel
    examples/fodo/analysis_fodo.py
    examples/fodo/plot_fodo.py
)


# MPI-Parallel FODO Cell ######################################################
#
add_impactx_test(FODO.MPI
    examples/fodo/input_fodo.in
    ON   # ImpactX MPI-parallel
    examples/fodo/analysis_fodo.py
    examples/fodo/plot_fodo.py
)
add_impactx_test(FODO.MPI.py
    examples/fodo/run_fodo.py
    ON  # ImpactX MPI-parallel
    examples/fodo/analysis_fodo.py
    examples/fodo/plot_fodo.py
)


# FODO channel using envelope tracking ##############################
#
# w/o space charge
add_impactx_test(FODO.envelope
    examples/fodo/input_fodo_envelope.in
    ON   # ImpactX MPI-parallel
    examples/fodo/analysis_fodo_envelope.py
    OFF  # no plot script yet
)
add_impactx_test(FODO.envelope.py
    examples/fodo/run_fodo_envelope.py
    OFF   # ImpactX MPI-parallel
    examples/fodo/analysis_fodo_envelope.py
    OFF  # no plot script yet
)

# FODO Cell using nonlinear tracking ################################################
#
add_impactx_test(FODO-exact
    examples/fodo/input_fodo_exact.in
    ON  # ImpactX MPI-parallel
    examples/fodo/analysis_fodo_exact.py
    OFF  # no plot script yet
)
add_impactx_test(FODO-exact.py
    examples/fodo/run_fodo_exact.py
    ON  # ImpactX MPI-parallel
    examples/fodo/analysis_fodo_exact.py
    OFF  # no plot script yet
)

# Python: FODO Cell w/ programmed element for the Drifts ######################
#
add_impactx_test(FODO.programmable.py
    examples/fodo_programmable/run_fodo_programmable.py
    OFF  # ImpactX MPI-parallel
    examples/fodo_programmable/analysis_fodo.py
    examples/fodo_programmable/plot_fodo.py
)


# Python MADX: FODO Cell ######################################################
#
# copy MAD-X lattice file
file(COPY ${ImpactX_SOURCE_DIR}/examples/fodo/fodo.madx
     DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/FODO.MADX.py)

add_impactx_test(FODO.MADX.py
    examples/fodo/run_fodo_madx.py
    OFF  # ImpactX MPI-parallel
    examples/fodo/analysis_fodo.py
    examples/fodo/plot_fodo.py
)


# Python: FODO Cell w/ custom linear element  #################################
#
add_impactx_test(FODO.userdef
    examples/fodo_userdef/input_fodo_userdef.in
    ON   # ImpactX MPI-parallel
    examples/fodo_userdef/analysis_fodo.py
    examples/fodo_userdef/plot_fodo.py
)
add_impactx_test(FODO.userdef.py
    examples/fodo_userdef/run_fodo_userdef.py
    OFF  # ImpactX MPI-parallel
    examples/fodo_userdef/analysis_fodo.py
    examples/fodo_userdef/plot_fodo.py
)


# FODO Channel ################################################################
#
add_impactx_test(FODO.channel
    examples/fodo_channel/input_fodo.in
    OFF  # ImpactX MPI-parallel
    examples/fodo_channel/analysis_fodo.py
    examples/fodo_channel/plot_fodo.py
)
add_impactx_test(FODO.channel.py
    examples/fodo_channel/run_fodo.py
    OFF  # ImpactX MPI-parallel
    examples/fodo_channel/analysis_fodo.py
    examples/fodo_channel/plot_fodo.py
)
label_impactx_test(FODO_channel slow)
label_impactx_test(FODO_channel.py slow)


# Chicane #####################################################################
#
add_impactx_test(chicane
    examples/chicane/input_chicane.in
    OFF  # ImpactX MPI-parallel
    examples/chicane/analysis_chicane.py
    examples/chicane/plot_chicane.py
)


# Python: Chicane #############################################################
#
add_impactx_test(chicane.py
    examples/chicane/run_chicane.py
    OFF  # ImpactX MPI-parallel
    examples/chicane/analysis_chicane.py
    examples/chicane/plot_chicane.py
)


# Python MADX: Chicane ########################################################
#
# copy MAD-X lattice file
file(COPY ${ImpactX_SOURCE_DIR}/examples/chicane/chicane.madx DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/chicane.MADX.py)

add_impactx_test(chicane.MADX.py
    examples/chicane/run_chicane_madx.py
    OFF  # ImpactX MPI-parallel
    examples/chicane/analysis_chicane.py
    examples/chicane/plot_chicane.py
)


# Constant Focusing Channel ###################################################
#
add_impactx_test(cfchannel
    examples/cfchannel/input_cfchannel.in
    OFF  # ImpactX MPI-parallel
    examples/cfchannel/analysis_cfchannel.py
    OFF  # no plot script yet
)
add_impactx_test(cfchannel.py
    examples/cfchannel/run_cfchannel.py
    OFF  # ImpactX MPI-parallel
    examples/cfchannel/analysis_cfchannel.py
    OFF  # no plot script yet
)


# Constant Focusing Channel with Space Charge #################################
#
add_impactx_test(cfchannel_spacecharge_mlmg
    examples/cfchannel/input_cfchannel_10nC_mlmg.in
    OFF  # ImpactX MPI-parallel
    examples/cfchannel/analysis_cfchannel_10nC.py
    OFF  # no plot script yet
)
if(ImpactX_FFT)
    add_impactx_test(cfchannel_spacecharge_fft
        examples/cfchannel/input_cfchannel_10nC_fft.in
        OFF  # ImpactX MPI-parallel
        examples/cfchannel/analysis_cfchannel_10nC.py
        OFF  # no plot script yet
    )
endif()


# Python: Constant Focusing Channel with Space Charge #########################
#
add_impactx_test(cfchannel_spacecharge_mlmg.py
    examples/cfchannel/run_cfchannel_10nC_mlmg.py
    OFF  # ImpactX MPI-parallel
    examples/cfchannel/analysis_cfchannel_10nC.py
    OFF  # no plot script yet
)
if(ImpactX_FFT)
    add_impactx_test(cfchannel_spacecharge_fft.py
        examples/cfchannel/run_cfchannel_10nC_fft.py
        OFF  # ImpactX MPI-parallel
        examples/cfchannel/analysis_cfchannel_10nC.py
        OFF  # no plot script yet
    )
endif()


# 6D Gaussian Distribution from Twiss Test ####################################
#
add_impactx_test(gaussian_twiss
    examples/distgen/input_gaussian_twiss.in
    OFF  # ImpactX MPI-parallel
    examples/distgen/analysis_gaussian.py
    OFF  # no plot script yet
)


# Python: 6D Gaussian Distribution  w/ Twiss input ############################
#
add_impactx_test(gaussian_twiss.py
    examples/distgen/run_gaussian_twiss.py
    OFF  # ImpactX MPI-parallel
    examples/distgen/analysis_gaussian.py
    OFF  # no plot script yet
)


# K-V Distribution from Twiss Test ############################################
#
add_impactx_test(kvdist_twiss
    examples/distgen/input_kvdist_twiss.in
    OFF  # ImpactX MPI-parallel
    examples/distgen/analysis_kvdist.py
    OFF  # no plot script yet
)
add_impactx_test(kvdist_twiss.py
    examples/distgen/run_kvdist_twiss.py
    OFF  # ImpactX MPI-parallel
    examples/distgen/analysis_kvdist.py
    OFF  # no plot script yet
)


# FODO Cell with RF ###########################################################
#
add_impactx_test(FODO.RF
    examples/fodo_rf/input_fodo_rf.in
    OFF  # ImpactX MPI-parallel
    examples/fodo_rf/analysis_fodo_rf.py
    OFF  # no plot script yet
)
add_impactx_test(FODO.RF.py
    examples/fodo_rf/run_fodo_rf.py
    OFF  # ImpactX MPI-parallel
    examples/fodo_rf/analysis_fodo_rf.py
    OFF  # no plot script yet
)


# 4D Kurth Distribution Test ##################################################
#
add_impactx_test(kurth4d
    examples/distgen/input_kurth4d.in
    OFF  # ImpactX MPI-parallel
    examples/distgen/analysis_kurth4d.py
    OFF  # no plot script yet
)
add_impactx_test(kurth4d.py
    examples/distgen/run_kurth4d.py
    OFF  # ImpactX MPI-parallel
    examples/distgen/analysis_kurth4d.py
    OFF  # no plot script yet
)


# Semi-Gaussian Distribution Test #############################################
#
add_impactx_test(semigaussian
    examples/distgen/input_semigaussian.in
    OFF  # ImpactX MPI-parallel
    examples/distgen/analysis_semigaussian.py
    OFF  # no plot script yet
)
add_impactx_test(semigaussian.py
    examples/distgen/run_semigaussian.py
    OFF  # ImpactX MPI-parallel
    examples/distgen/analysis_semigaussian.py
    OFF  # no plot script yet
)


# Chain of Multipoles Test ####################################################
#
add_impactx_test(multipole
    examples/multipole/input_multipole.in
    OFF  # ImpactX MPI-parallel
    examples/multipole/analysis_multipole.py
    OFF  # no plot script yet
)
add_impactx_test(multipole.py
    examples/multipole/run_multipole.py
    OFF  # ImpactX MPI-parallel
    examples/multipole/analysis_multipole.py
    OFF  # not plotting script yet
)


# Expanding Beam Test #########################################################
#
add_impactx_test(expanding_beam_mlmg
    examples/expanding_beam/input_expanding_mlmg.in
    OFF  # ImpactX MPI-parallel
    examples/expanding_beam/analysis_expanding.py
    OFF  # no plot script yet
)
add_impactx_test(expanding_beam_mlmg.py
    examples/expanding_beam/run_expanding_mlmg.py
    OFF  # ImpactX MPI-parallel
    examples/expanding_beam/analysis_expanding.py
    OFF  # no plot script yet
)
if(ImpactX_FFT)
    add_impactx_test(expanding_beam_fft
        examples/expanding_beam/input_expanding_fft.in
        OFF  # ImpactX MPI-parallel
        examples/expanding_beam/analysis_expanding.py
        OFF  # no plot script yet
    )
    add_impactx_test(expanding_beam_fft.py
        examples/expanding_beam/run_expanding_fft.py
        OFF  # ImpactX MPI-parallel
        examples/expanding_beam/analysis_expanding.py
        OFF  # no plot script yet
    )
endif()

# Expanding Beam Test Using Envelope Tracking ##################################
#
add_impactx_test(expanding_beam_envelope
    examples/expanding_beam/input_expanding_envelope.in
    OFF  # ImpactX MPI-parallel
    examples/expanding_beam/analysis_expanding_envelope.py
    OFF  # no plot script yet
)
add_impactx_test(expanding_beam_envelope.py
    examples/expanding_beam/run_expanding_envelope.py
    OFF  # ImpactX MPI-parallel
    examples/expanding_beam/analysis_expanding_envelope.py
    OFF  # no plot script yet
)

# IOTA Nonlinear Focusing Channel Test ########################################
#
add_impactx_test(iotalens
    examples/iota_lens/input_iotalens.in
    OFF  # ImpactX MPI-parallel
    examples/iota_lens/analysis_iotalens.py
    OFF  # no plot script yet
)
add_impactx_test(iotalens.py
    examples/iota_lens/run_iotalens.py
    OFF  # ImpactX MPI-parallel
    examples/iota_lens/analysis_iotalens.py
    OFF  # not plotting script yet
)


# IOTA Linear Lattice Test ####################################################
#
add_impactx_test(iotalattice.MPI
    examples/iota_lattice/input_iotalattice.in
    ON   # ImpactX MPI-parallel
    examples/iota_lattice/analysis_iotalattice.py
    OFF  # no plot script yet
)
add_impactx_test(iotalattice.MPI.py
    examples/iota_lattice/run_iotalattice.py
    ON   # ImpactX MPI-parallel
    examples/iota_lattice/analysis_iotalattice.py
    OFF  # no plot script yet
)


# Kurth Beam in a Periodic Channel Test ###################################################
#
# w/o space charge
add_impactx_test(kurth_periodic
    examples/kurth/input_kurth_periodic.in
    OFF  # ImpactX MPI-parallel
    examples/kurth/analysis_kurth_periodic.py
    OFF  # no plot script yet
)
add_impactx_test(kurth_periodic.py
    examples/kurth/run_kurth_periodic.py
    OFF  # ImpactX MPI-parallel
    examples/kurth/analysis_kurth_periodic.py
    OFF  # no plot script yet
)

# w/ space charge
add_impactx_test(kurth_10nC_periodic
    examples/kurth/input_kurth_10nC_periodic.in
    OFF  # ImpactX MPI-parallel
    examples/kurth/analysis_kurth_10nC_periodic.py
    OFF  # no plot script yet
)
add_impactx_test(kurth_10nC_periodic.py
    examples/kurth/run_kurth_10nC_periodic.py
    OFF  # ImpactX MPI-parallel
    examples/kurth/analysis_kurth_10nC_periodic.py
    OFF  # no plot script yet
)

# Kurth Beam in a Periodic Channel with Envelope Tracking ######################
#
# with space charge
add_impactx_test(kurth_10nC_periodic_envelope
    examples/kurth/input_kurth_10nC_periodic_envelope.in
    OFF  # ImpactX MPI-parallel
    examples/kurth/analysis_kurth_10nC_periodic_envelope.py
    OFF  # no plot script yet
)
add_impactx_test(kurth_10nC_periodic_envelope.py
    examples/kurth/run_kurth_10nC_periodic_envelope.py
    OFF  # ImpactX MPI-parallel
    examples/kurth/analysis_kurth_10nC_periodic_envelope.py
    OFF  # no plot script yet
)

# Acceleration by RF Cavities #################################################
#
add_impactx_test(rfcavity.MPI
    examples/rfcavity/input_rfcavity.in
    ON   # ImpactX MPI-parallel
    examples/rfcavity/analysis_rfcavity.py
    OFF  # no plot script yet
)
add_impactx_test(rfcavity.MPI.py
    examples/rfcavity/run_rfcavity.py
    ON   # ImpactX MPI-parallel
    examples/rfcavity/analysis_rfcavity.py
    OFF  # no plot script yet
)

# Acceleration by RF Cavities (Reference Particle) ####################################
#
add_impactx_test(rfcavity-ref-part
    examples/rfcavity/input_rfcavity_ref_part.in
    OFF   # ImpactX MPI-parallel
    examples/rfcavity/analysis_rfcavity_ref_part.py
    OFF  # no plot script yet
)
add_impactx_test(rfcavity-ref-part.py
    examples/rfcavity/run_rfcavity_ref_part.py
    OFF   # ImpactX MPI-parallel
    examples/rfcavity/analysis_rfcavity_ref_part.py
    OFF  # no plot script yet
)

# Ideal, Hard-Edge Solenoid ###################################################
#
# w/o space charge
add_impactx_test(solenoid
    examples/solenoid/input_solenoid.in
    OFF  # ImpactX MPI-parallel
    examples/solenoid/analysis_solenoid.py
    OFF  # no plot script yet
)
add_impactx_test(solenoid.py
    examples/solenoid/run_solenoid.py
    OFF  # ImpactX MPI-parallel
    examples/solenoid/analysis_solenoid.py
    OFF  # no plot script yet
)

# copy MAD-X lattice file
file(COPY ${ImpactX_SOURCE_DIR}/examples/solenoid/solenoid.madx
     DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/solenoid.MADX.py)
add_impactx_test(solenoid.MADX.py
    examples/solenoid/run_solenoid_madx.py
    OFF  # ImpactX MPI-parallel
    examples/solenoid/analysis_solenoid.py
    OFF  # no plot script yet
)


# Ideal, Hard-Edge Solenoid (Restart + 270 degrees) ###########################
#
# w/o space charge
add_impactx_test(solenoid.restart
    examples/solenoid_restart/input_solenoid.in
    OFF  # ImpactX MPI-parallel
    examples/solenoid_restart/analysis_solenoid.py
    OFF  # no plot script yet
)
set_property(TEST solenoid.restart.run APPEND PROPERTY DEPENDS "solenoid.run")
if(ImpactX_TEST_CLEANUP)
    # do not clean up dependency test before current test is completed
    set_property(TEST solenoid.cleanup APPEND PROPERTY DEPENDS "solenoid.restart.cleanup")
endif()

add_impactx_test(solenoid.restart.py
    examples/solenoid_restart/run_solenoid.py
    OFF  # ImpactX MPI-parallel
    examples/solenoid_restart/analysis_solenoid.py
    OFF  # no plot script yet
)
if(ImpactX_PYTHON)
    set_property(TEST solenoid.restart.py.run APPEND PROPERTY DEPENDS "solenoid.py.run")
endif()
if(ImpactX_TEST_CLEANUP)
    # do not clean up dependency test before current test is completed
    set_property(TEST solenoid.py.cleanup APPEND PROPERTY DEPENDS "solenoid.restart.py.cleanup")
endif()


# Pole-face rotation ##########################################################
#
# w/o space charge
add_impactx_test(rotation
    examples/rotation/input_rotation.in
    ON   # ImpactX MPI-parallel
    examples/rotation/analysis_rotation.py
    OFF  # no plot script yet
)
add_impactx_test(rotation.py
    examples/rotation/run_rotation.py
    ON   # ImpactX MPI-parallel
    examples/rotation/analysis_rotation.py
    OFF  # no plot script yet
)


# Soft-Edge Solenoid ##########################################################
#
# w/o space charge
add_impactx_test(solenoid_softedge
    examples/solenoid_softedge/input_solenoid_softedge.in
    OFF  # ImpactX MPI-parallel
    examples/solenoid_softedge/analysis_solenoid_softedge.py
    OFF  # no plot script yet
)
add_impactx_test(solenoid_softedge.py
    examples/solenoid_softedge/run_solenoid_softedge.py
    OFF  # ImpactX MPI-parallel
    examples/solenoid_softedge/analysis_solenoid_softedge.py
    OFF  # no plot script yet
)


# Soft-Edge Quadrupole ########################################################
#
# w/o space charge
add_impactx_test(quadrupole_softedge
    examples/quadrupole_softedge/input_quadrupole_softedge.in
    OFF  # ImpactX MPI-parallel
    examples/quadrupole_softedge/analysis_quadrupole_softedge.py
    OFF  # no plot script yet
)
add_impactx_test(quadrupole_softedge.py
    examples/quadrupole_softedge/run_quadrupole_softedge.py
    OFF  # ImpactX MPI-parallel
    examples/quadrupole_softedge/analysis_quadrupole_softedge.py
    OFF  # no plot script yet
)


# FODO Cell with Chromatic Elements ###########################################
#
# w/o space charge
add_impactx_test(FODO.chromatic
    examples/fodo_chromatic/input_fodo_chr.in
    ON  # ImpactX MPI-parallel
    examples/fodo_chromatic/analysis_fodo_chr.py
    OFF  # no plot script yet
)
add_impactx_test(FODO.chromatic.py
    examples/fodo_chromatic/run_fodo_chr.py
    OFF  # ImpactX MPI-parallel
    examples/fodo_chromatic/analysis_fodo_chr.py
    OFF  # no plot script yet
)


# Positron Channel ############################################################
#
# w/o space charge
add_impactx_test(positron_channel
    examples/positron_channel/input_positron.in
    ON  # ImpactX MPI-parallel
    examples/positron_channel/analysis_positron.py
    OFF  # no plot script yet
)
add_impactx_test(positron_channel.py
    examples/positron_channel/run_positron.py
    OFF  # ImpactX MPI-parallel
    examples/positron_channel/analysis_positron.py
    OFF  # no plot script yet
)


# Cyclotron ###################################################################
#
# w/o space charge
add_impactx_test(cyclotron
    examples/cyclotron/input_cyclotron.in
    ON  # ImpactX MPI-parallel
    examples/cyclotron/analysis_cyclotron.py
    OFF  # no plot script yet
)
add_impactx_test(cyclotron.py
    examples/cyclotron/run_cyclotron.py
    OFF  # ImpactX MPI-parallel
    examples/cyclotron/analysis_cyclotron.py
    OFF  # no plot script yet
)


# Combined-function bend ######################################################
#
# w/o space charge
add_impactx_test(cfbend
    examples/cfbend/input_cfbend.in
    ON  # ImpactX MPI-parallel
    examples/cfbend/analysis_cfbend.py
    OFF  # no plot script yet
)
add_impactx_test(cfbend.py
    examples/cfbend/run_cfbend.py
    OFF  # ImpactX MPI-parallel
    examples/cfbend/analysis_cfbend.py
    OFF  # no plot script yet
)


# Ballistic compression #######################################################
#
# w/o space charge
add_impactx_test(compression
    examples/compression/input_compression.in
    ON  # ImpactX MPI-parallel
    examples/compression/analysis_compression.py
    OFF  # no plot script yet
)
add_impactx_test(compression.py
    examples/compression/run_compression.py
    OFF  # ImpactX MPI-parallel
    examples/compression/analysis_compression.py
    OFF  # no plot script yet
)


# Kicker test #################################################################
#
# w/o space charge
add_impactx_test(kicker
    examples/kicker/input_kicker.in
    ON   # ImpactX MPI-parallel
    examples/kicker/analysis_kicker.py
    OFF  # no plot script yet
)
add_impactx_test(kicker.py
    examples/kicker/run_kicker.py
    OFF   # ImpactX MPI-parallel
    examples/kicker/analysis_kicker.py
    OFF  # no plot script yet
)
# copy MAD-X lattice file
file(COPY ${ImpactX_SOURCE_DIR}/examples/kicker/kicker.madx
        DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/kicker_madx.py)
add_impactx_test(kicker_madx.py
    examples/kicker/run_kicker_madx.py
    OFF   # ImpactX MPI-parallel
    examples/kicker/analysis_kicker.py
    OFF  # no plot script yet
)
# copy MAD-X lattice file
file(COPY ${ImpactX_SOURCE_DIR}/examples/kicker/hvkicker.madx
        DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/hvkicker_madx.py)
add_impactx_test(hvkicker_madx.py
    examples/kicker/run_hvkicker_madx.py
    OFF   # ImpactX MPI-parallel
    examples/kicker/analysis_kicker.py
    OFF  # no plot script yet
)


# FODO + RF EPAC2004 ##########################################################
#
# with space charge
add_impactx_test(FODO.rf.sc
    examples/epac2004_benchmarks/input_fodo_rf_SC.in
    ON  # ImpactX MPI-parallel
    examples/epac2004_benchmarks/analysis_fodo_rf_SC.py
    OFF  # no plot script yet
)
add_impactx_test(FODO.rf.sc.py
    examples/epac2004_benchmarks/run_fodo_rf_SC.py
    OFF  # ImpactX MPI-parallel
    examples/epac2004_benchmarks/analysis_fodo_rf_SC.py
    OFF  # no plot script yet
)

# FODO + RF EPAC2004 Using Envelope Tracking ##################################
#
# with space charge
add_impactx_test(FODO.rf.sc.envelope
    examples/epac2004_benchmarks/input_fodo_rf_SC_envelope.in
    ON  # ImpactX MPI-parallel
    examples/epac2004_benchmarks/analysis_fodo_rf_SC_envelope.py
    OFF  # no plot script yet
)
add_impactx_test(FODO.rf.sc.envelope.py
    examples/epac2004_benchmarks/run_fodo_rf_SC_envelope.py
    OFF  # ImpactX MPI-parallel
    examples/epac2004_benchmarks/analysis_fodo_rf_SC_envelope.py
    OFF  # no plot script yet
)

# Thermal Beam EPAC2004 #######################################################
#
# with space charge
add_impactx_test(thermal
    examples/epac2004_benchmarks/input_thermal.in
    ON  # ImpactX MPI-parallel
    examples/epac2004_benchmarks/analysis_thermal.py
    OFF  # no plot script yet
)
label_impactx_test(thermal slow)


# Bithermal Beam EPAC2004 #####################################################
#
# with space charge
add_impactx_test(bithermal
    examples/epac2004_benchmarks/input_bithermal.in
    ON  # ImpactX MPI-parallel
    examples/epac2004_benchmarks/analysis_bithermal.py
    examples/epac2004_benchmarks/plot_bithermal.py
)
add_impactx_test(bithermal.py
    examples/epac2004_benchmarks/run_bithermal.py
    ON  # ImpactX MPI-parallel
    examples/epac2004_benchmarks/analysis_bithermal.py
    examples/epac2004_benchmarks/plot_bithermal.py
)


# IOTA s-dependent nonlinear lens test ########################################
#
# w/o space charge
add_impactx_test(IOTA_nll
    examples/iota_lens/input_iotalens_sdep.in
    ON   # ImpactX MPI-parallel
    examples/iota_lens/analysis_iotalens_sdep.py
    OFF  # no plot script yet
)
add_impactx_test(IOTA_nll.py
    examples/iota_lens/run_iotalens_sdep.py
    OFF   # ImpactX MPI-parallel
    examples/iota_lens/analysis_iotalens_sdep.py
    OFF  # no plot script yet
)


# IOTA nonlinear lattice test #################################################
#
# w/o space charge
add_impactx_test(IOTA_lattice
    examples/iota_lattice/input_iotalattice_sdep.in
    ON   # ImpactX MPI-parallel
    examples/iota_lattice/analysis_iotalattice_sdep.py
    OFF  # no plot script yet
)
add_impactx_test(IOTA_lattice.py
    examples/iota_lattice/run_iotalattice_sdep.py
    OFF   # ImpactX MPI-parallel
    examples/iota_lattice/analysis_iotalattice_sdep.py
    OFF  # no plot script yet
)


# Thin dipole #################################################################
#
# w/o space charge
add_impactx_test(thin_dipole
    examples/thin_dipole/input_thin_dipole.in
    ON  # ImpactX MPI-parallel
    examples/thin_dipole/analysis_thin_dipole.py
    OFF  # no plot script yet
)
add_impactx_test(thin_dipole.py
    examples/thin_dipole/run_thin_dipole.py
    OFF  # ImpactX MPI-parallel
    examples/thin_dipole/analysis_thin_dipole.py
    OFF  # no plot script yet
)


# Aperture collimation ########################################################
#
# w/o space charge
add_impactx_test(aperture
    examples/aperture/input_aperture.in
    ON  # ImpactX MPI-parallel
    examples/aperture/analysis_aperture.py
    OFF  # no plot script yet
)
add_impactx_test(aperture.py
    examples/aperture/run_aperture.py
    OFF  # ImpactX MPI-parallel
    examples/aperture/analysis_aperture.py
    OFF  # no plot script yet
)


# Absorber collimation ########################################################
#
# w/o space charge
add_impactx_test(absorber
    examples/aperture/input_absorber.in
    ON  # ImpactX MPI-parallel
    examples/aperture/analysis_absorber.py
    OFF  # no plot script yet
)
add_impactx_test(absorber.py
    examples/aperture/run_absorber.py
    OFF  # ImpactX MPI-parallel
    examples/aperture/analysis_absorber.py
    OFF  # no plot script yet
)


# Apochromat drift-quad example ##########################################################
#
# w/o space charge
add_impactx_test(apochromat
    examples/apochromatic/input_apochromatic.in
    ON  # ImpactX MPI-parallel
    examples/apochromatic/analysis_apochromatic.py
    OFF  # no plot script yet
)
add_impactx_test(apochromat.py
    examples/apochromatic/run_apochromatic.py
    OFF  # ImpactX MPI-parallel
    examples/apochromatic/analysis_apochromatic.py
    OFF  # no plot script yet
)


# Apochromat drift-plasma lens example ###################################################
#
# w/o space charge
add_impactx_test(apochromat_pl
    examples/apochromatic/input_apochromatic_pl.in
    ON  # ImpactX MPI-parallel
    examples/apochromatic/analysis_apochromatic_pl.py
    OFF  # no plot script yet
)
add_impactx_test(apochromat_pl.py
    examples/apochromatic/run_apochromatic_pl.py
    OFF  # ImpactX MPI-parallel
    examples/apochromatic/analysis_apochromatic_pl.py
    OFF  # no plot script yet
)


# Misalignment Quad ###########################################################
#
# w/o space charge
add_impactx_test(alignment
    examples/alignment/input_alignment.in
    ON   # ImpactX MPI-parallel
    examples/alignment/analysis_alignment.py
    OFF  # no plot script yet
)
add_impactx_test(alignment.py
    examples/alignment/run_alignment.py
    OFF  # ImpactX MPI-parallel
    examples/alignment/analysis_alignment.py
    OFF  # no plot script yet
)


# Tune Calculation ###########################################################
#
# w/o space charge
add_impactx_test(FODO.tune
    examples/fodo_tune/input_fodo_tune.in
    ON   # ImpactX MPI-parallel
    examples/fodo_tune/analysis_fodo_tune.py
    OFF  # no plot script yet
)
add_impactx_test(FODO.tune.py
    examples/fodo_tune/run_fodo_tune.py
    OFF  # ImpactX MPI-parallel
    examples/fodo_tune/analysis_fodo_tune.py
    OFF  # no plot script yet
)


# Optimize Triplet ###########################################################
#
# w/o space charge
add_impactx_test(triplet.py
    examples/optimize_triplet/run_triplet.py
    OFF  # ImpactX MPI-parallel
    examples/optimize_triplet/analysis_triplet.py
    examples/optimize_triplet/plot_triplet.py
)
label_impactx_test(triplet.py slow)


# Initialize a Beam from Arrays ##############################################
#
# w/o space charge
add_impactx_test(initialize_from_array.py
    examples/initialize_from_array/run_from_array.py
    OFF  # ImpactX MPI-parallel
    examples/initialize_from_array/analyze_from_array.py
    examples/initialize_from_array/visualize_from_array.py
)


# Achromatic Spectrometer ####################################################
#
# w/o space charge
add_impactx_test(spectrometer
    examples/achromatic_spectrometer/input_spectrometer.in
    ON   # ImpactX MPI-parallel
    examples/achromatic_spectrometer/analysis_spectrometer.py
    OFF  # no plot script yet
)
add_impactx_test(spectrometer.py
    examples/achromatic_spectrometer/run_spectrometer.py
    OFF  # ImpactX MPI-parallel
    examples/achromatic_spectrometer/analysis_spectrometer.py
    OFF  # no plot script yet
)


# Chicane with CSR ###########################################################
#
if(ImpactX_FFT)
    add_impactx_test(chicane_csr
        examples/chicane/input_chicane_csr.in
        OFF  # ImpactX MPI-parallel
        examples/chicane/analysis_chicane_csr.py
        OFF  # no plot script yet
    )
    add_impactx_test(chicane_csr.py
        examples/chicane/run_chicane_csr.py
        OFF  # ImpactX MPI-parallel
        examples/chicane/analysis_chicane_csr.py
        OFF  # no plot script yet
    )
endif()


# Dogleg #####################################################################
#
# w/o space charge
add_impactx_test(dogleg
    examples/dogleg/input_dogleg.in
    ON   # ImpactX MPI-parallel
    examples/dogleg/analysis_dogleg.py
    OFF  # no plot script yet
)
add_impactx_test(dogleg.py
    examples/dogleg/run_dogleg.py
    OFF  # ImpactX MPI-parallel
    examples/dogleg/analysis_dogleg.py
    OFF  # no plot script yet
)


# Coupled Optics #############################################################
#
# w/o space charge
add_impactx_test(coupled-optics
    examples/coupled_optics/input_coupled_optics.in
    ON   # ImpactX MPI-parallel
    examples/coupled_optics/analysis_coupled_optics.py
    OFF  # no plot script yet
)
add_impactx_test(coupled-optics.py
    examples/coupled_optics/run_coupled_optics.py
    OFF  # ImpactX MPI-parallel
    examples/coupled_optics/analysis_coupled_optics.py
    OFF  # no plot script yet
)


# Aperture with Periodic Masking #########################################################
#
# w/o space charge
add_impactx_test(aperture-periodic
    examples/aperture/input_aperture_periodic.in
    ON   # ImpactX MPI-parallel
    examples/aperture/analysis_aperture_periodic.py
    OFF  # no plot script yet
)
add_impactx_test(aperture-periodic.py
    examples/aperture/run_aperture_periodic.py
    OFF  # ImpactX MPI-parallel
    examples/aperture/analysis_aperture_periodic.py
    OFF  # no plot script yet
)


# Rotation in the transverse plane #########################################################
#
# w/o space charge
add_impactx_test(rotation-xy
    examples/rotation/input_rotation_xy.in
    ON   # ImpactX MPI-parallel
    examples/rotation/analysis_rotation_xy.py
    OFF  # no plot script yet
)
add_impactx_test(rotation-xy.py
    examples/rotation/run_rotation_xy.py
    OFF  # ImpactX MPI-parallel
    examples/rotation/analysis_rotation_xy.py
    OFF  # no plot script yet
)


# PIP-II linac segment #########################################################
#
# with space charge
add_impactx_test(linac-segment
    examples/linac_segment/input_linac_segment.in
    ON   # ImpactX MPI-parallel
    examples/linac_segment/analysis_linac_segment.py
    OFF  # no plot script yet
)
add_impactx_test(linac-segment.py
    examples/linac_segment/run_linac_segment.py
    OFF  # ImpactX MPI-parallel
    examples/linac_segment/analysis_linac_segment.py
    OFF  # no plot script yet
)


# Cyclotron with dynamical losses ###############################################
#
# with space charge
add_impactx_test(cyclotron-loss
    examples/cyclotron/input_cyclotron_loss.in
    ON   # ImpactX MPI-parallel
    examples/cyclotron/analysis_cyclotron_loss.py
    OFF  # no plot script yet
)
add_impactx_test(cyclotron-loss.py
    examples/cyclotron/run_cyclotron_loss.py
    OFF  # ImpactX MPI-parallel
    examples/cyclotron/analysis_cyclotron_loss.py
    OFF  # no plot script yet
)


# Drift with transverse aperture ################################################
#
# w/o space charge
add_impactx_test(aperture-thick
    examples/aperture/input_aperture_thick.in
    ON   # ImpactX MPI-parallel
    examples/aperture/analysis_aperture_thick.py
    OFF  # no plot script yet
)
add_impactx_test(aperture-thick.py
    examples/aperture/run_aperture_thick.py
    OFF  # ImpactX MPI-parallel
    examples/aperture/analysis_aperture_thick.py
    OFF  # no plot script yet
)


# Iteration of a linear one-turn map #########################################
#
# w/o space charge
add_impactx_test(linear-map
    examples/linear_map/input_map.in
    ON   # ImpactX MPI-parallel
    examples/linear_map/analysis_map.py
    OFF  # no plot script yet
)
add_impactx_test(linear-map.py
    examples/linear_map/run_map.py
    ON  # ImpactX MPI-parallel
    examples/linear_map/analysis_map.py
    OFF  # no plot script yet
)


# PyTorch Surrogate: Staged LPA ##############################################
#
add_impactx_test(pytorch_surrogate_model
    examples/pytorch_surrogate_model/run_ml_surrogate_15_stage.py
    OFF  # ImpactX MPI-parallel
    examples/pytorch_surrogate_model/analyze_ml_surrogate_15_stage.py
    examples/pytorch_surrogate_model/visualize_ml_surrogate_15_stage.py
)
label_impactx_test(pytorch_surrogate_model slow)


# IOTA s-dependent nonlinear lens with aperture ##############################
#
# w/o space charge
add_impactx_test(IOTA_nll_aperture
    examples/iota_lens/input_iotalens_sdep_aperture.in
    ON   # ImpactX MPI-parallel
    examples/iota_lens/analysis_iotalens_sdep_aperture.py
    OFF  # no plot script yet
)
add_impactx_test(IOTA_nll_aperture.py
    examples/iota_lens/run_iotalens_sdep_aperture.py
    OFF   # ImpactX MPI-parallel
    examples/iota_lens/analysis_iotalens_sdep_aperture.py
    OFF  # no plot script yet
)


# IOTA linear lattice using envelope tracking ############################
#
add_impactx_test(iotalattice.envelope
    examples/iota_lattice/input_iotalattice_envelope.in
    ON   # ImpactX MPI-parallel
    examples/iota_lattice/analysis_iotalattice_envelope.py
    OFF  # no plot script yet
)
add_impactx_test(iotalattice.envelope.py
    examples/iota_lattice/run_iotalattice_envelope.py
    ON   # ImpactX MPI-parallel
    examples/iota_lattice/analysis_iotalattice_envelope.py
    OFF  # no plot script yet
)


# Expanding beam scraping against a vacuum pipe ##############################
#
# w/o space charge
add_impactx_test(scraping
    examples/scraping_beam/input_scraping.in
    ON   # ImpactX MPI-parallel
    examples/scraping_beam/analysis_scraping.py
    OFF  # no plot script yet
)
add_impactx_test(scraping.py
    examples/scraping_beam/run_scraping.py
    OFF   # ImpactX MPI-parallel
    examples/scraping_beam/analysis_scraping.py
    OFF  # no plot script yet
)

# FODO cell with 2D space charge using envelope tracking ######################
#
# with space charge
add_impactx_test(FODO.envelope.sc
    examples/fodo_space_charge/input_fodo_envelope_sc.in
    ON   # ImpactX MPI-parallel
    examples/fodo_space_charge/analysis_fodo_envelope_sc.py
    OFF  # no plot script yet
)
add_impactx_test(FODO.envelope.sc.py
    examples/fodo_space_charge/run_fodo_envelope_sc.py
    OFF   # ImpactX MPI-parallel
    examples/fodo_space_charge/analysis_fodo_envelope_sc.py
    OFF  # no plot script yet
)

# FODO cell with Gaussian 3D space charge using particle tracking #############
#
add_impactx_test(FODO.Gauss3d.sc
    examples/fodo_space_charge/input_fodo_Gauss3D_sc.in
    ON   # ImpactX MPI-parallel
    examples/fodo_space_charge/analysis_fodo_Gauss3D_sc.py
    OFF  # no plot script yet
)
add_impactx_test(FODO.Gauss3d.sc.py
    examples/fodo_space_charge/run_fodo_Gauss3D_sc.py
    ON   # ImpactX MPI-parallel
    examples/fodo_space_charge/analysis_fodo_Gauss3D_sc.py
    OFF  # no plot script yet
)

# A single bend with incoherent synchrotron radiation ######################
#
# no space charge
add_impactx_test(bend-isr
    examples/incoherent_synchrotron/input_bend_isr.in
    ON   # ImpactX MPI-parallel
    examples/incoherent_synchrotron/analysis_bend_isr.py
    OFF  # no plot script yet
)
add_impactx_test(bend-isr.py
    examples/incoherent_synchrotron/run_bend_isr.py
    OFF   # ImpactX MPI-parallel
    examples/incoherent_synchrotron/analysis_bend_isr.py
    OFF  # no plot script yet
)

# Bending Dipole with Zero Field Strength ######################
#
# no space charge
add_impactx_test(zero-bend-inverse
    examples/reversibility/input_zero_bend_inverse.in
    ON   # ImpactX MPI-parallel
    examples/reversibility/analysis_zero_bend_inverse.py
    OFF  # no plot script yet
)
add_impactx_test(zero-bend-inverse.py
    examples/reversibility/run_zero_bend_inverse.py
    OFF   # ImpactX MPI-parallel
    examples/reversibility/analysis_zero_bend_inverse.py
    OFF  # no plot script yet
)

# Testing Charge and Field Sign Consistency ######################
#
# no space charge
add_impactx_test(charge-sign
    examples/charge_sign/input_charge_sign.in
    ON   # ImpactX MPI-parallel
    examples/charge_sign/analysis_charge_sign.py
    OFF  # no plot script yet
)
add_impactx_test(charge-sign.py
    examples/charge_sign/run_charge_sign.py
    ON   # ImpactX MPI-parallel
    examples/charge_sign/analysis_charge_sign.py
    OFF  # no plot script yet
)

# Chicane with ISR ###########################################################
#
add_impactx_test(chicane_isr
    examples/chicane/input_chicane_isr.in
    OFF  # ImpactX MPI-parallel
    examples/chicane/analysis_chicane_isr.py
    OFF  # no plot script yet
)
add_impactx_test(chicane_isr.py
    examples/chicane/run_chicane_isr.py
    OFF  # ImpactX MPI-parallel
    examples/chicane/analysis_chicane_isr.py
    OFF  # no plot script yet
)

# Symplectic Integration in an Exact Nonlinear Quadrupole ##############
#
# no space charge
add_impactx_test(exact-quad
    examples/symplectic_integration/input_exact_quad.in
    ON   # ImpactX MPI-parallel
    examples/symplectic_integration/analysis_exact_quad.py
    OFF  # no plot script yet
)
add_impactx_test(exact-quad.py
    examples/symplectic_integration/run_exact_quad.py
    ON   # ImpactX MPI-parallel
    examples/symplectic_integration/analysis_exact_quad.py
    OFF  # no plot script yet
)

# FODO Channel with Quads Treated as Exact Multipoles ##############
#
# no space charge
add_impactx_test(fodo-multipole
    examples/symplectic_integration/input_fodo_multipole.in
    ON   # ImpactX MPI-parallel
    examples/symplectic_integration/analysis_fodo_multipole.py
    OFF  # no plot script yet
)
add_impactx_test(fodo_multipole.py
    examples/symplectic_integration/run_fodo_multipole.py
    ON   # ImpactX MPI-parallel
    examples/symplectic_integration/analysis_fodo_multipole.py
    OFF  # no plot script yet
)

# Symplectic Integration in a Long Sextupole ##############
#
# no space charge
file(COPY ${ImpactX_SOURCE_DIR}/examples/symplectic_integration/initial_coords.csv
        DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/long-sextupole.py)
file(COPY ${ImpactX_SOURCE_DIR}/examples/symplectic_integration/final_coords.csv
        DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/long-sextupole.py)
add_impactx_test(long-sextupole.py
    examples/symplectic_integration/run_sextupole.py
    ON   # ImpactX MPI-parallel
    examples/symplectic_integration/analysis_sextupole.py
    OFF  # no plot script yet
)

# Dogleg in Reverse  ##############
#
# no space charge
add_impactx_test(dogleg-reverse
    examples/dogleg/input_dogleg_reverse.in
    ON   # ImpactX MPI-parallel
    examples/dogleg/analysis_dogleg_reverse.py
    OFF  # no plot script yet
)
add_impactx_test(dogleg-reverse.py
    examples/dogleg/run_dogleg_reverse.py
    ON   # ImpactX MPI-parallel
    examples/dogleg/analysis_dogleg_reverse.py
    OFF  # no plot script yet
)

# Dogleg with Energy Jitter ##############
#
# no space charge
add_impactx_test(dogleg-jitter
    examples/dogleg/input_dogleg_jitter.in
    ON   # ImpactX MPI-parallel
    examples/dogleg/analysis_dogleg_jitter.py
    OFF  # no plot script yet
)
add_impactx_test(dogleg-jitter.py
    examples/dogleg/run_dogleg_jitter.py
    ON   # ImpactX MPI-parallel
    examples/dogleg/analysis_dogleg_jitter.py
    OFF  # no plot script yet
)

# HTU Beamline Model ##############
#
# no space charge
file(COPY ${ImpactX_SOURCE_DIR}/examples/htu_beamline/htu_lattice.py
        DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/htu-beamline.py)
add_impactx_test(htu-beamline.py
    examples/htu_beamline/run_impactx.py
    ON   # ImpactX MPI-parallel
    examples/htu_beamline/analysis_htu_beamline.py
    OFF  # no plot script yet
)

# HTU Beamline Model with Off-Energy Transport ##############
#
file(COPY ${ImpactX_SOURCE_DIR}/examples/htu_beamline/htu_lattice.py
        DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/htu-beamline-offenergy1.py)
file(COPY ${ImpactX_SOURCE_DIR}/examples/htu_beamline/htu_lattice.py
        DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/htu-beamline-offenergy2.py)
add_impactx_test(htu-beamline-offenergy1.py
    examples/htu_beamline/run_impactx_offenergy1.py
    ON   # ImpactX MPI-parallel
    examples/htu_beamline/analysis_htu_beamline_offenergy.py
    OFF  # no plot script yet
)
add_impactx_test(htu-beamline-offenergy2.py
    examples/htu_beamline/run_impactx_offenergy2.py
    ON   # ImpactX MPI-parallel
    examples/htu_beamline/analysis_htu_beamline_offenergy.py
    OFF  # no plot script yet
)

# Symplectic Integration in an Exact Combined-Function Bend ##############
#
# no space charge
add_impactx_test(exact-cfbend
    examples/symplectic_integration/input_exact_cfbend.in
    ON   # ImpactX MPI-parallel
    examples/symplectic_integration/analysis_exact_cfbend.py
    OFF  # no plot script yet
)
add_impactx_test(exact-cfbend.py
    examples/symplectic_integration/run_exact_cfbend.py
    ON   # ImpactX MPI-parallel
    examples/symplectic_integration/analysis_exact_cfbend.py
    OFF  # no plot script yet
)

# Symplectic Integration in an Exact Combined-Function Bend (Multipole test) ##############
#
# no space charge
add_impactx_test(exact-cfbend-multipole
    examples/symplectic_integration/input_exact_cfbend_multipole.in
    ON   # ImpactX MPI-parallel
    examples/symplectic_integration/analysis_exact_cfbend_multipole.py
    OFF  # no plot script yet
)
add_impactx_test(exact-cfbend_multipole.py
    examples/symplectic_integration/run_exact_cfbend_multipole.py
    ON   # ImpactX MPI-parallel
    examples/symplectic_integration/analysis_exact_cfbend_multipole.py
    OFF  # no plot script yet
)

# FODO Cell Based on PALS-Compliant Lattice Input ##############
#
# copy PALS lattice file
file(COPY ${ImpactX_SOURCE_DIR}/examples/pals/fodo.pals.yaml
     DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/fodo-pals.py)
add_impactx_test(fodo-pals.py
    examples/pals/run_fodo_pals.py
    ON   # ImpactX MPI-parallel
    examples/pals/analysis_fodo_pals.py
    OFF  # no plot script yet
)

# Active Plasma Lens tests ####################################
#
# ChrPlasmaLens, k = 0 (acting like a drift)
file(COPY ${ImpactX_SOURCE_DIR}/examples/active_plasma_lens/run_APL.py
     DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/APL_ChrPlasmaLens_zero.py)
file(COPY ${ImpactX_SOURCE_DIR}/examples/active_plasma_lens/analysis_APL_ChrPlasmaLens.py
  DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/APL_ChrPlasmaLens_zero.py)
file(COPY ${ImpactX_SOURCE_DIR}/examples/active_plasma_lens/plot_APL_ChrPlasmaLens.py
  DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/APL_ChrPlasmaLens_zero.py)

add_impactx_test(APL_ChrPlasmaLens_zero.py
  examples/active_plasma_lens/run_APL_ChrPlasmaLens_zero.py
  OFF # No MPI needed
  examples/active_plasma_lens/analysis_APL_ChrPlasmaLens_zero.py
  examples/active_plasma_lens/plot_APL_ChrPlasmaLens_zero.py
)

# ChrPlasmaLens, mg = -1000T/m (focusing)
file(COPY ${ImpactX_SOURCE_DIR}/examples/active_plasma_lens/run_APL.py
     DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/APL_ChrPlasmaLens_focusing.py)
file(COPY ${ImpactX_SOURCE_DIR}/examples/active_plasma_lens/analysis_APL_ChrPlasmaLens.py
  DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/APL_ChrPlasmaLens_focusing.py)
file(COPY ${ImpactX_SOURCE_DIR}/examples/active_plasma_lens/plot_APL_ChrPlasmaLens.py
  DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/APL_ChrPlasmaLens_focusing.py)

add_impactx_test(APL_ChrPlasmaLens_focusing.py
  examples/active_plasma_lens/run_APL_ChrPlasmaLens_focusing.py
  OFF # No MPI needed
  examples/active_plasma_lens/analysis_APL_ChrPlasmaLens_focusing.py
  examples/active_plasma_lens/plot_APL_ChrPlasmaLens_focusing.py
)

# ChrPlasmaLens, mg = 1000T/m (defocusing)
file(COPY ${ImpactX_SOURCE_DIR}/examples/active_plasma_lens/run_APL.py
     DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/APL_ChrPlasmaLens_defocusing.py)
file(COPY ${ImpactX_SOURCE_DIR}/examples/active_plasma_lens/analysis_APL_ChrPlasmaLens.py
  DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/APL_ChrPlasmaLens_defocusing.py)
file(COPY ${ImpactX_SOURCE_DIR}/examples/active_plasma_lens/plot_APL_ChrPlasmaLens.py
  DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/APL_ChrPlasmaLens_defocusing.py)

add_impactx_test(APL_ChrPlasmaLens_defocusing.py
  examples/active_plasma_lens/run_APL_ChrPlasmaLens_defocusing.py
  OFF # No MPI needed
  examples/active_plasma_lens/analysis_APL_ChrPlasmaLens_defocusing.py
  examples/active_plasma_lens/plot_APL_ChrPlasmaLens_defocusing.py
  )

# Exactly-solvable (non-uniform) soft-edge solenoid ##############
#
# copy PALS lattice file
add_impactx_test(solenoid-softedge-solvable.py
    examples/solenoid_softedge/run_solenoid_softedge_solvable.py
    ON   # ImpactX MPI-parallel
    examples/solenoid_softedge/analysis_solenoid_softedge_solvable.py
    OFF  # no plot script yet
)
