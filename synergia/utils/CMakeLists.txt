set(synergia_parallel_utils_src parallel_utils.cc commxx.cc logger.cc
                                simple_timer.cc base64.cpp)
add_library(synergia_parallel_utils SHARED ${synergia_parallel_utils_src})
add_library(synergia_parallel_utils_static STATIC
            ${synergia_parallel_utils_src})
target_link_libraries(
  synergia_parallel_utils cereal::cereal MPI::MPI_C
  $<$<STREQUAL:${BUILD_FD_SPACE_CHARGE_SOLVER},ON>:PkgConfig::PETSC>)
target_link_libraries(
  synergia_parallel_utils_static cereal::cereal MPI::MPI_C
  $<$<STREQUAL:${BUILD_FD_SPACE_CHARGE_SOLVER},ON>:PkgConfig::PETSC>)
target_link_options(synergia_parallel_utils PRIVATE ${LINKER_OPTIONS})
target_link_options(synergia_parallel_utils_static PRIVATE ${LINKER_OPTIONS})

add_library(synergia_serialization cereal_files.cc)
target_link_libraries(synergia_serialization synergia_parallel_utils)
target_link_options(synergia_serialization PRIVATE ${LINKER_OPTIONS})

add_library(synergia_test_main STATIC catch_test_main.cc)
target_link_libraries(
  synergia_test_main
  PRIVATE Catch2::Catch2 ${kokkos_libs} MPI::MPI_C
          $<$<STREQUAL:${BUILD_FD_SPACE_CHARGE_SOLVER},ON>:PkgConfig::PETSC>)

add_library(synergia_hdf5_utils hdf5_file.cc hdf5_misc.cc)
target_link_libraries(synergia_hdf5_utils ${HDF5_LIBRARIES}
                      synergia_parallel_utils ${kokkos_libs})
target_include_directories(synergia_hdf5_utils PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_options(synergia_hdf5_utils PRIVATE ${LINKER_OPTIONS})

if("${ENABLE_KOKKOS_BACKEND}" STREQUAL "CUDA")
  find_package(CUDAToolkit REQUIRED)
  set(FFT_SRC distributed_fft2d_cuda.cc distributed_fft3d_cuda.cc
              distributed_fft3d_rect_cuda.cc)
  set(FFT_LIB CUDA::cufft)
else()
  set(FFT_SRC distributed_fft2d_fftw.cc distributed_fft3d_fftw.cc
              distributed_fft3d_rect_fftw.cc)
  set(FFT_LIB ${PARALLEL_FFTW_LIBRARIES})
endif()
add_library(synergia_distributed_fft ${FFT_SRC})
if("${ENABLE_KOKKOS_BACKEND}" STREQUAL "OpenMP")
  target_include_directories(synergia_distributed_fft
                             PUBLIC ${FFTW3_INCLUDE_DIR})
endif()
target_link_libraries(synergia_distributed_fft ${FFT_LIB}
                      synergia_parallel_utils ${kokkos_libs})
target_link_options(synergia_distributed_fft PRIVATE ${LINKER_OPTIONS})

add_library(synergia_command_line command_line_arg.cc)
target_link_libraries(synergia_command_line)
target_link_options(synergia_command_line PRIVATE ${LINKER_OPTIONS})

add_library(lsexpr SHARED lsexpr.cc)
add_library(lsexpr_static STATIC lsexpr.cc)
target_link_options(lsexpr PRIVATE ${LINKER_OPTIONS})
target_link_options(lsexpr_static PRIVATE ${LINKER_OPTIONS})

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/synergia_version.h.in"
               "${CMAKE_CURRENT_BINARY_DIR}/synergia_version.h" IMMEDIATE @ONLY)

add_library(synergia_utils INTERFACE)
target_include_directories(synergia_utils
                           INTERFACE "${CMAKE_CURRENT_BINARY_DIR}")

add_library(synergia_boost INTERFACE)
target_include_directories(synergia_boost
                           INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

if(BUILD_PYTHON_BINDINGS)
  pybind11_add_module(parallel_utils MODULE NO_EXTRAS parallel_utils_pywrap.cc)
  target_link_libraries(parallel_utils PRIVATE synergia_parallel_utils)

  pybind11_add_module(pylsexpr MODULE NO_EXTRAS lsexpr_pywrap.cc)
  target_link_libraries(pylsexpr PRIVATE lsexpr)

  pybind11_add_module(utils MODULE NO_EXTRAS utils_pywrap.cc)
  target_link_libraries(
    utils
    PRIVATE ${kokkos_libs}
            $<$<STREQUAL:${BUILD_FD_SPACE_CHARGE_SOLVER},ON>:PkgConfig::PETSC>
            MPI::MPI_C)

  copy_file(__init__.py utils)

  install(FILES __init__.py DESTINATION ${PYTHON_INSTALL_DIR}/synergia/utils)
  install(TARGETS parallel_utils
          DESTINATION ${PYTHON_INSTALL_DIR}/synergia/utils)
endif()

install(TARGETS synergia_parallel_utils synergia_distributed_fft
                synergia_serialization synergia_command_line lsexpr
        DESTINATION ${SYNERGIA3_LIB_INSTALL_DIR})
install(TARGETS synergia_hdf5_utils DESTINATION ${SYNERGIA3_LIB_INSTALL_DIR})

install(
  FILES command_line_arg.h
        commxx.h
        container_conversions.h
        distributed_fft3d.h
        distributed_fft2d.h
        fast_int_floor.h
        floating_point.h
        gsvector.h
        hdf5_misc.h
        hdf5_file.h
        hdf5_serial_writer.h
        hdf5_writer.h
        complex_error_function.h
        kokkos_views.h
        kokkos_utils.h
        parallel_utils.h
        simple_timer.h
        cereal.h
        cereal_files.h
        digits.h
        logger.h
        lsexpr.h
        synergia_config.h
        utils.h
  DESTINATION ${SYNERGIA3_INCLUDE_INSTALL_DIR}/synergia/utils)
if(BUILD_PYTHON_BINDINGS)
  install(FILES __init__.py DESTINATION ${PYTHON_INSTALL_DIR}/synergia/utils)
  install(TARGETS parallel_utils pylsexpr utils
          DESTINATION ${PYTHON_INSTALL_DIR}/synergia/utils)
endif()

add_subdirectory(tests)
