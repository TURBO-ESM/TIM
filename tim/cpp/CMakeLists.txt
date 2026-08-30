# TIM's C++ infrastructure layer (namespace TIM)
#
# Built two ways:
#   - As part of the full TIM build, via the top-level CMakeLists.txt
#   - Or, standalone by C++ consumers (e.g. protoMOMxx):
#     add_subdirectory() or `cmake -S tim/cpp`

cmake_minimum_required(VERSION 3.24)
project(TIM_cpp LANGUAGES C CXX)

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

add_library(tim STATIC
    core/tim_abort.cpp
    core/tim_domain.cpp
    core/tim_runtime.cpp
    tim_coms_infra.cpp
    tim_coms_infra_C_API.cpp
    tim_profile.cpp
)
add_library(TIM::tim ALIAS tim)

target_include_directories(tim PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
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

# Install/export rules: ON for full TIM builds (the top level sets
# TIM_CPP_INSTALL before add_subdirectory) and for direct configures of
# this directory; OFF by default when a consumer embeds this directory
option(TIM_CPP_INSTALL "Install the tim library and its headers" ${PROJECT_IS_TOP_LEVEL})
if(TIM_CPP_INSTALL)
  install(TARGETS tim
    EXPORT TIMExports
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)
  install(FILES
      tim_coms_infra.hpp
      tim_coms_infra_C_API.h
      tim_profile.hpp
      DESTINATION include/tim)
  install(FILES
      core/tim_abort.hpp
      core/tim_domain.hpp
      core/tim_runtime.hpp
      core/tim_stagger.hpp
      DESTINATION include/tim/core)
endif()
