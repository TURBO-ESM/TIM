#***********************************************************************
#*                   GNU Lesser General Public License
#*
#* This file is part of the GFDL Flexible Modeling System (FMS).
#*
#* FMS is free software: you can redistribute it and/or modify it under
#* the terms of the GNU Lesser General Public License as published by
#* the Free Software Foundation, either version 3 of the License, or (at
#* your option) any later version.
#*
#* FMS is distributed in the hope that it will be useful, but WITHOUT
#* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#* for more details.
#*
#* You should have received a copy of the GNU Lesser General Public
#* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************

# Copyright (c) GFDL, @underwoo

cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

# add build type for debug, overrides default flags (set with $FCFLAGS, $CFLAGS)
set(CMAKE_Fortran_FLAGS_DEBUG)

# Define the CMake project
project(TIM
  VERSION 2023.01.0
  DESCRIPTION  "TURBO Infrastructure for MOM (TIM) is an AMReX alternative to the TURBO project. It aims to provide the infrastructure that FMS2 did."
  LANGUAGES C CXX Fortran)

include(GNUInstallDirs)


if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release|RelWithDebInfo|MinSizeRel)$")
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

if(NOT CMAKE_C_COMPILER_ID MATCHES "^(Intel|GNU|Clang|IntelLLVM)$")
  message(
    WARNING "Compiler not officially supported: ${CMAKE_C_COMPILER_ID}")
endif()

if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|GNU|IntelLLVM)$")
  message(
    WARNING "Compiler not officially supported: ${CMAKE_Fortran_COMPILER_ID}")
endif()

# Append directory that contains CMake Modules for building TIM
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Build options
option(OPENMP "Build TIM with OpenMP support"        OFF)
option(32BIT  "Build 32-bit (r4) TIM library"        OFF)
option(64BIT  "Build 64-bit (r8) TIM library"         ON)
option(FPIC   "Build with position independent code" OFF)
option(TIM_ADD_TEST_TARGETS "Build C++ unit test targets (test_tim)"          OFF)
option(TIM_FETCH_DEPS       "Fetch missing test dependencies (GoogleTest)"     ON)

# Options for compiler definitions
option(INTERNAL_FILE_NML     "Enable compiler definition -DINTERNAL_FILE_NML"      ON)
option(ENABLE_QUAD_PRECISION "Enable compiler definition -DENABLE_QUAD_PRECISION"  ON)
option(GFS_PHYS              "Enable compiler definition -DGFS_PHYS"              OFF)
option(LARGEFILE             "Enable compiler definition -Duse_LARGEFILE"         OFF)

if(32BIT)
  list(APPEND kinds "r4")
endif()
if(64BIT)
  list(APPEND kinds "r8")
endif()
if(NOT kinds)
  message(STATUS "Single Precision 32BIT: ${32BIT}")
  message(STATUS "Double Precision 64BIT: ${64BIT}")
  message(FATAL_ERROR "Either 32BIT or 64BIT should be ON")
endif()

# Find dependencies
find_package(MPI REQUIRED COMPONENTS C CXX Fortran)
find_package(NetCDF REQUIRED COMPONENTS C Fortran)
find_package(AMReX REQUIRED)
find_package(PIO REQUIRED COMPONENTS C)

# Check for the OpenMP library and set the required compile flags
if (OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS C Fortran)
endif()

# Enables position independent code (i.e., -fPIC)
if (FPIC)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif ()

# Collect FMS Fortran source files
list(APPEND tim_fortran_src_files
  affinity/fms_affinity.F90
  astronomy/astronomy.F90
  axis_utils/axis_utils2.F90
  block_control/block_control.F90
  constants/constants.F90
  constants/fmsconstants.F90
  coupler/atmos_ocean_fluxes.F90
  coupler/coupler_types.F90
  coupler/ensemble_manager.F90
  data_override/get_grid_version.F90
  data_override/data_override.F90
  diag_manager/diag_axis.F90
  diag_manager/diag_data.F90
  diag_manager/diag_manager.F90
  diag_manager/diag_output.F90
  diag_manager/diag_table.F90
  diag_manager/diag_util.F90
  diag_manager/fms_diag_time_reduction.F90
  diag_manager/fms_diag_outfield.F90
  diag_manager/fms_diag_elem_weight_procs.F90
  diag_manager/fms_diag_fieldbuff_update.F90
  diag_manager/fms_diag_bbox.F90
  field_manager/field_manager.F90
  field_manager/fm_util.F90
  fms/fms_io.F90
  fms/fms.F90
  fms2_io/fms_io_utils.F90
  fms2_io/fms_netcdf_domain_io.F90
  fms2_io/fms2_io.F90
  fms2_io/netcdf_io.F90
  horiz_interp/horiz_interp_bicubic.F90
  horiz_interp/horiz_interp_bilinear.F90
  horiz_interp/horiz_interp_conserve.F90
  horiz_interp/horiz_interp_spherical.F90
  horiz_interp/horiz_interp_type.F90
  horiz_interp/horiz_interp.F90
  memutils/memutils.F90
  mosaic2/grid2.F90
  mosaic2/mosaic2.F90
  mpp/mpp.F90
  mpp/mpp_data.F90
  mpp/mpp_domains.F90
  mpp/mpp_io.F90
  mpp/mpp_memutils.F90
  mpp/mpp_parameter.F90
  platform/platform.F90
  random_numbers/mersennetwister.F90
  random_numbers/random_numbers.F90
  string_utils/fms_string_utils.F90
  time_interp/time_interp_external2.F90
  time_interp/time_interp.F90
  time_manager/get_cal_time.F90
  time_manager/time_manager.F90
  libFMS.F90
  tim/fortran/tim_coms_infra_interface.F90
)

# Collect FMS C source files
list(APPEND tim_c_src_files
  diag_manager/mppnccombine.c
  affinity/affinity.c
  mosaic2/create_xgrid.c
  mosaic2/gradient_c2l.c
  mosaic2/interp.c
  mosaic2/mosaic_util.c
  mosaic2/read_mosaic.c
  mpp/mpp_memuse.c
  string_utils/fms_string_utils_binding.c
)

# Collect FMS header files
list(APPEND tim_header_files
  include/file_version.h
  include/fms_platform.h
)

# Standard FMS compiler definitions
list(APPEND tim_defs
  use_libMPI
  use_netCDF)

# check gettid
include(CheckFunctionExists)
check_function_exists(gettid HAVE_GETTID)
if(HAVE_GETTID)
  list(APPEND tim_defs HAVE_GETTID)
endif()

# Additional (optional) compiler definitions
if(NOT CONSTANTS)
    set(CONSTANTS GFDL)
endif()
if(CONSTANTS STREQUAL "GFS")
  list(APPEND tim_defs GFS_CONSTANTS)
elseif(CONSTANTS STREQUAL "GEOS")
  list(APPEND tim_defs GEOS_CONSTANTS)
elseif(CONSTANTS STREQUAL "GFDL")
  list(APPEND tim_defs GFDL_CONSTANTS)
else()
  message(FATAL_ERROR "CONSTANTS=${CONSTANTS} option not supported")
endif()

if(GFS_PHYS)
  list(APPEND tim_defs GFS_PHYS)
endif()

if(INTERNAL_FILE_NML)
  list(APPEND tim_defs INTERNAL_FILE_NML)
endif()

if(ENABLE_QUAD_PRECISION)
  list(APPEND tim_defs ENABLE_QUAD_PRECISION)
endif()

if(LARGEFILE)
  list(APPEND tim_defs use_LARGEFILE)
endif()

# Precision-based compiler definitions
if(32BIT)
  list(APPEND r4_defs OVERLOAD_R4 OVERLOAD_R8)
endif()

# Add platform specific compiler definitions
if(APPLE)
  list(APPEND tim_defs __APPLE__)
endif()

# Obtain compiler-specific flags
include(fms_compiler_flags)

# TIM C++ infrastructure layer
add_library(tim STATIC
    tim/cpp/tim_coms_infra.cpp
    tim/cpp/tim_coms_infra_C_API.cpp
    tim/cpp/tim_profile.cpp
)
target_include_directories(tim PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/tim/cpp>
    $<INSTALL_INTERFACE:include/tim>)
# tim's C++ object files and AMReX::amrex reference the C++ standard
# library (e.g. __cxa_call_terminate). Position the runtime *after*
# AMReX::amrex so it lands after libamrex_3d.a in the consumer's link
# line -- otherwise GNU ld with --as-needed (default on some Linux
# distros, incl. NCAR Derecho) drops libstdc++ before amrex's undefined
# refs are seen.
target_link_libraries(tim PUBLIC AMReX::amrex MPI::MPI_CXX
                                 $<$<CXX_COMPILER_ID:GNU,IntelLLVM>:stdc++>
                                 $<$<CXX_COMPILER_ID:Clang,AppleClang>:c++>)
# PUBLIC so consumers including tim_profile.hpp / tim_coms_infra_C_API.h
# (installed below) get the C++20 requirement propagated.
target_compile_features(tim PUBLIC cxx_std_20)
set_target_properties(tim PROPERTIES CXX_EXTENSIONS OFF)
install(FILES
    tim/cpp/tim_coms_infra_C_API.h
    tim/cpp/tim_profile.hpp
    DESTINATION include/tim)

foreach(kind ${kinds})

  set(libTgt tim_${kind})
  set(includeDir "include_${kind}")
  set(moduleDir "${CMAKE_CURRENT_BINARY_DIR}/${includeDir}")

  # C
  add_library(${libTgt}_c OBJECT ${tim_c_src_files})

  target_include_directories(${libTgt}_c PRIVATE include)
  target_compile_definitions(${libTgt}_c PRIVATE "${tim_defs}")

  target_link_libraries(${libTgt}_c PRIVATE NetCDF::NetCDF_C
                                            MPI::MPI_C)

  if(OpenMP_C_FOUND)
    target_link_libraries(${libTgt}_c PRIVATE OpenMP::OpenMP_C)
  endif()

  # Fortran
  add_library(${libTgt}_f OBJECT ${tim_fortran_src_files})

  target_include_directories(${libTgt}_f PRIVATE include
                                                 fms
                                                 fms2_io/include
                                                 string_utils/include
                                                 mpp/include
                                                 diag_manager/include
                                                 constants)
  target_compile_definitions(${libTgt}_f PRIVATE "${tim_defs}")
  target_compile_definitions(${libTgt}_f PRIVATE "${${kind}_defs}")

  if(${kind}_flags)
    set_target_properties(${libTgt}_f PROPERTIES COMPILE_FLAGS
                                               "${${kind}_flags}")
  endif()
  set_target_properties(${libTgt}_f PROPERTIES Fortran_MODULE_DIRECTORY
                                               ${moduleDir})

  target_link_libraries(${libTgt}_f PRIVATE NetCDF::NetCDF_Fortran
                                            MPI::MPI_Fortran)

  if(OpenMP_Fortran_FOUND)
    target_link_libraries(${libTgt}_f PRIVATE OpenMP::OpenMP_Fortran)
  endif()

  # Check if gnu 10 or higher with mpich
  if ( CMAKE_Fortran_COMPILER_VERSION MATCHES "1[0-9]\.[0-9]*\.[0-9]*" AND CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    if(MPI_C_COMPILER MATCHES ".*mpich.*" )
      message(STATUS "Adding -fallow-argument-mismatch flag to compile with GCC >=10 and MPICH")
      set_target_properties(${libTgt}_f PROPERTIES COMPILE_FLAGS "-fallow-argument-mismatch -w")
    endif()
  endif()

  # FMS (C + Fortran)
  add_library(${libTgt} STATIC $<TARGET_OBJECTS:${libTgt}_c>
                               $<TARGET_OBJECTS:${libTgt}_f>)

  target_include_directories(${libTgt} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fms>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fms2_io/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/string_utils/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mpp/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/diag_manager/include>)

  target_include_directories(${libTgt} INTERFACE
    $<BUILD_INTERFACE:${moduleDir}>
    $<INSTALL_INTERFACE:${includeDir}>)

  target_compile_definitions(${libTgt} PRIVATE "${tim_defs}")
  target_compile_definitions(${libTgt} PRIVATE "${${kind}_defs}")

  target_link_libraries(${libTgt} PUBLIC NetCDF::NetCDF_C
                                         NetCDF::NetCDF_Fortran
                                         MPI::MPI_Fortran)

  if(OpenMP_Fortran_FOUND)
    target_link_libraries(${libTgt} PRIVATE OpenMP::OpenMP_Fortran)
  endif()

  target_link_libraries(${libTgt} PUBLIC tim)

  add_library(TIM::${libTgt} ALIAS ${libTgt})

  list(APPEND LIB_TARGETS ${libTgt})
  install(DIRECTORY ${moduleDir}    DESTINATION .)
  install(FILES ${tim_header_files} DESTINATION ${includeDir})

endforeach()


install(
  TARGETS ${LIB_TARGETS} tim
  EXPORT TIMExports
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

### Package config
include(CMakePackageConfigHelpers)
set(CONFIG_INSTALL_DESTINATION lib/cmake/tim)

export(EXPORT TIMExports
  NAMESPACE TIM::
  FILE ${CMAKE_CURRENT_BINARY_DIR}/${CONFIG_INSTALL_DESTINATION}/tim-targets.cmake)

# Versioned find_package(TIM <version>) derives the filename from the
# config stem; with TIMConfig.cmake it must be TIMConfigVersion.cmake.
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/${CONFIG_INSTALL_DESTINATION}/TIMConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${CONFIG_INSTALL_DESTINATION}/TIMConfigVersion.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/TIMConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${CONFIG_INSTALL_DESTINATION}/TIMConfig.cmake
  INSTALL_DESTINATION ${CONFIG_INSTALL_DESTINATION})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${CONFIG_INSTALL_DESTINATION}/TIMConfig.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})

# FindNetCDF.cmake must sit next to TIMConfig.cmake so the module path
# adjustment in TIMConfig.cmake.in can find it.
configure_file(cmake/FindNetCDF.cmake ${CMAKE_CURRENT_BINARY_DIR}/${CONFIG_INSTALL_DESTINATION}/FindNetCDF.cmake COPYONLY)
install(FILES cmake/FindNetCDF.cmake DESTINATION ${CONFIG_INSTALL_DESTINATION})

install(EXPORT TIMExports
  NAMESPACE TIM::
  FILE tim-targets.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})

# C++ unit tests (reuses the AMReX/PIO/tim targets defined above).
if(TIM_ADD_TEST_TARGETS)
  enable_testing()
  add_subdirectory(test_tim)
endif()
