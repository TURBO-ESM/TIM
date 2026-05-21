# Standalone CMake project for MOM C++ unit tests.
#
# Configure & run with, e.g.:
#   cmake -S test_mom -B build/test_mom
#   cmake --build build/test_mom
#   ctest --test-dir build/test_mom --output-on-failure
#
# AMReX is brought in as a source subdirectory from the turbo-stack submodule
# (../../../amrex relative to this file). Override with -DAMREX_SOURCE_DIR=...
# or point CMAKE_PREFIX_PATH at an installed AMReX to use find_package instead.
#
# Captured Fortran I/O data is consumed via the `--data-dir=PATH` CLI flag on
# the test executable. ctest invocations get it injected from
# TEST_MOM_DATA_DIR (CMake cache variable) -- override on the configure line
# to point at a different capture directory.

cmake_minimum_required(VERSION 3.20)

project(test_mom CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(AMReX_GPU_BACKEND STREQUAL "CUDA")
    enable_language(CUDA)
endif()

# --- AMReX ----------------------------------------------------------------
# Prefer an installed AMReX if available; otherwise build the in-tree source.
set(AMREX_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../amrex"
    CACHE PATH "Path to the AMReX source tree")

find_package(AMReX QUIET)
if(NOT AMReX_FOUND)
    if(NOT EXISTS "${AMREX_SOURCE_DIR}/CMakeLists.txt")
        message(FATAL_ERROR
            "AMReX not found via find_package and AMREX_SOURCE_DIR "
            "('${AMREX_SOURCE_DIR}') does not contain a CMakeLists.txt.")
    endif()
    message(STATUS "Building AMReX from source at ${AMREX_SOURCE_DIR}")
    set(AMReX_SPACEDIM 3 CACHE STRING "" FORCE)
    add_subdirectory(${AMREX_SOURCE_DIR}
                     ${CMAKE_CURRENT_BINARY_DIR}/_amrex)
endif()

# --- GTest ----------------------------------------------------------------
option(TEST_MOM_FETCH_GTEST
       "Download GoogleTest via FetchContent if find_package fails" ON)

find_package(GTest QUIET)
if(NOT GTest_FOUND)
    if(TEST_MOM_FETCH_GTEST)
        include(FetchContent)
        FetchContent_Declare(
            googletest
            URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
            DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        )
        FetchContent_MakeAvailable(googletest)
    else()
        message(FATAL_ERROR
            "GTest not found. Install GoogleTest or set CMAKE_PREFIX_PATH "
            "to include its install location, or configure with "
            "-DTEST_MOM_FETCH_GTEST=ON to download it automatically.")
    endif()
endif()

enable_testing()
include(GoogleTest)

# --- Captured-data directory (used by --data-dir CLI flag) ---------------
set(TEST_MOM_DATA_DIR ""
    CACHE PATH "Directory containing the captured Fortran I/O .bin/.meta files")
if(NOT TEST_MOM_DATA_DIR)
    message(FATAL_ERROR
        "TEST_MOM_DATA_DIR is not set. Pass\n"
        "  -DTEST_MOM_DATA_DIR=/path/to/captured/fixtures\n"
        "on the cmake configure line.")
endif()

# --- MOM C++ kernels library (target: mom_continuity) ---------------------
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../mom/cpp
                 ${CMAKE_CURRENT_BINARY_DIR}/mom_cpp)

# --- mom_continuity_ppm test ----------------------------------------------
add_executable(test_mom_continuity_ppm
    test_mom_main.cpp
    common/captured_io.cpp
    common/amrex_assertions.cpp
    mom_continuity_ppm/test_mom_continuity_ppm.cpp
)

# Every source in this target includes AMReX headers. With the CUDA backend
# those headers use __host__/__device__ qualifiers and CUDA kernel-launch
# syntax, so nvcc must compile them. CMake won't infer this from .cpp, so we
# set the language explicitly.
if(AMReX_GPU_BACKEND STREQUAL "CUDA")
    set_source_files_properties(
        test_mom_main.cpp
        common/captured_io.cpp
        common/amrex_assertions.cpp
        mom_continuity_ppm/test_mom_continuity_ppm.cpp
        PROPERTIES LANGUAGE CUDA)
endif()

target_include_directories(test_mom_continuity_ppm
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/common
)
target_link_libraries(test_mom_continuity_ppm
    PRIVATE
        mom_continuity
        GTest::gtest
)
gtest_discover_tests(test_mom_continuity_ppm
    EXTRA_ARGS --data-dir=${TEST_MOM_DATA_DIR}
    # Discover at test-run time rather than build time. With the CUDA backend
    # the binary calls amrex::Initialize() which aborts on hosts without a GPU
    # (e.g. Derecho login nodes), so running it at POST_BUILD would break the
    # build.
    DISCOVERY_MODE PRE_TEST
)
