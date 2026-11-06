# TIM's C++ infrastructure layer (namespace TIM)
#
# Built two ways:
#   - As part of the full TIM build, via the top-level CMakeLists.txt
#   - Or, standalone by C++ consumers (e.g. protoMOMxx):
#     add_subdirectory() or `cmake -S tim/cpp`
#
# Targets:
#   TIM::tim_base  abort, Time, Stagger. MPI only; never links AMReX or PIO.
#   TIM::tim_io    the I/O layer: tim_base + PIO. Never links AMReX. Only when
#                  TIM_ENABLE_IO is on.
#   TIM::tim       the AMReX core (Runtime, Domain, checksums). Links tim_base
#                  only, so a consumer of tim alone carries no PIO requirement.

cmake_minimum_required(VERSION 3.24)
project(TIM_cpp LANGUAGES C CXX)

# ON for full TIM builds and the test suite (both set it before add_subdirectory)
# and for direct configures of this directory; OFF by default when a consumer
# embeds this directory, so embedding tim alone never resolves or fetches PIO.
option(TIM_ENABLE_IO "Build the AMReX-free I/O targets (requires PIO)" ${PROJECT_IS_TOP_LEVEL})
option(TIM_FETCH_DEPS "Fetch missing dependencies (GoogleTest, PIO)" ON)

if(NOT TARGET MPI::MPI_CXX)
  find_package(MPI REQUIRED COMPONENTS CXX)
endif()
if(NOT TARGET AMReX::amrex)
  # An AMReX built with AMReX_FORTRAN=ON requires the importing project to enable Fortran
  include(CheckLanguage)
  check_language(Fortran)
  if(CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
  endif()
  find_package(AMReX REQUIRED)
endif()

# Position the C++ runtime last on consumers' link lines: tim's object files and
# AMReX::amrex reference it (e.g. __cxa_call_terminate), and GNU ld with
# --as-needed (default on some Linux distros, incl. NCAR Derecho) drops it if it
# comes before the archives that need it. Fortran-linked consumers rely on this.
set(_tim_cxx_runtime
    $<$<CXX_COMPILER_ID:GNU,IntelLLVM>:stdc++>
    $<$<CXX_COMPILER_ID:Clang,AppleClang>:c++>)

# --- tim_base ------------------------------------------------------------------
add_library(tim_base STATIC
    core/tim_abort.cpp
    core/tim_time.cpp
)
add_library(TIM::tim_base ALIAS tim_base)

target_include_directories(tim_base PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/tim>)
target_link_libraries(tim_base PUBLIC MPI::MPI_CXX ${_tim_cxx_runtime})
target_compile_features(tim_base PUBLIC cxx_std_20)
set_target_properties(tim_base PROPERTIES CXX_EXTENSIONS OFF)

# --- tim_io --------------------------------------------------------------------
if(TIM_ENABLE_IO)
  # find_package(PIO) or FetchContent pio2_6_8 -> PIO::PIO_C
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindOrFetchPIO.cmake)

  add_library(tim_io INTERFACE)
  add_library(TIM::tim_io ALIAS tim_io)
  target_link_libraries(tim_io INTERFACE tim_base PIO::PIO_C)
endif()

# --- tim -----------------------------------------------------------------------
add_library(tim STATIC
    core/tim_domain.cpp
    core/tim_runtime.cpp
    tim_coms_infra.cpp
    tim_coms_infra_C_API.cpp
    tim_profile.cpp
)
add_library(TIM::tim ALIAS tim)

# tim_base supplies MPI, the include path and the C++20 requirement.
target_link_libraries(tim PUBLIC tim_base AMReX::amrex ${_tim_cxx_runtime})
set_target_properties(tim PROPERTIES CXX_EXTENSIONS OFF)

# Install/export rules: ON for full TIM builds (the top level sets
# TIM_CPP_INSTALL before add_subdirectory) and for direct configures of
# this directory; OFF by default when a consumer embeds this directory
option(TIM_CPP_INSTALL "Install the tim library and its headers" ${PROJECT_IS_TOP_LEVEL})
if(TIM_CPP_INSTALL)
  set(_tim_install_targets tim_base tim)
  if(TIM_ENABLE_IO)
    list(APPEND _tim_install_targets tim_io)
  endif()
  install(TARGETS ${_tim_install_targets}
    EXPORT TIMExports
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)
  # turbotmp_bridge_c_types.h is included by tim_coms_infra_C_API.h, so it
  # installs alongside it -- consumers of an installed TIM::tim need both.
  install(FILES
      tim_coms_infra.hpp
      tim_coms_infra_C_API.h
      tim_profile.hpp
      turbotmp_bridge_c_types.h
      DESTINATION include/tim)
  install(FILES
      core/tim_abort.hpp
      core/tim_domain.hpp
      core/tim_runtime.hpp
      core/tim_stagger.hpp
      core/tim_time.hpp
      DESTINATION include/tim/core)
endif()
