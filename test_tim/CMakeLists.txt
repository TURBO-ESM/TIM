# TIM C++ unit tests
#
# Run from the TIM root, with parallelio (PIO) and AMReX installs on CMAKE_PREFIX_PATH:
#   ./build.sh --tests                              # build + install + ctest
# or with raw CMake:
#   cmake -S . -B build -DTIM_ADD_TEST_TARGETS=ON \
#         -DCMAKE_PREFIX_PATH="<amrex-install>;$NCAR_ROOT_PARALLELIO"
#   cmake --build build --target test_tim_smoke
#   ctest --test-dir build -L unit --output-on-failure
# or standalone -- just tim/cpp and its tests, without the full TIM/FMS build
#   cmake -S test_tim -B build -DCMAKE_PREFIX_PATH=<amrex-install>
#   cmake --build build && ctest --test-dir build -L unit

# Standalone configuration: provide what the top-level build provides
# otherwise. AMReX comes from CMAKE_PREFIX_PATH. PIO comes from an install
# on CMAKE_PREFIX_PATH / $PIO, or a libpioc source build (FindOrFetchPIO).
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  cmake_minimum_required(VERSION 3.24)
  project(TIM_cpp_tests LANGUAGES C CXX)
  option(TIM_FETCH_DEPS "Fetch missing dependencies (GoogleTest, PIO)" ON)
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")
  include(FindOrFetchPIO)
  add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../tim/cpp" tim_cpp)
endif()

find_package(GTest QUIET)
if(NOT GTest_FOUND)
  if(TIM_FETCH_DEPS)
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
      "-DTIM_FETCH_DEPS=ON to download it automatically.")
  endif()
endif()
enable_testing()

if(AMReX_GPU_BACKEND STREQUAL "CUDA")
  enable_language(CUDA)
  set_source_files_properties(
    test_tim_main.cpp
    smoke/test_smoke.cpp
    runtime/test_runtime_owner.cpp
    runtime/test_runtime_guest.cpp
    domain/test_tim_domain.cpp
    PROPERTIES LANGUAGE CUDA)
endif()

add_executable(test_tim_smoke
  test_tim_main.cpp
  smoke/test_smoke.cpp
)

add_executable(test_tim_runtime_owner
  runtime/test_runtime_owner.cpp
)

add_executable(test_tim_runtime_guest
  runtime/test_runtime_guest.cpp
)

add_executable(test_tim_domain
  domain/test_tim_domain.cpp
)

target_link_libraries(test_tim_smoke
  GTest::gtest
  tim
  PIO::PIO_C
)

foreach(runtime_test test_tim_runtime_owner test_tim_runtime_guest test_tim_domain)
  target_link_libraries(${runtime_test}
    GTest::gtest
    tim
  )
endforeach()

include(GoogleTest)

gtest_discover_tests(test_tim_smoke
  DISCOVERY_MODE PRE_TEST
  PROPERTIES LABELS "unit")

gtest_discover_tests(test_tim_runtime_owner
  DISCOVERY_MODE PRE_TEST
  PROPERTIES LABELS "unit")

gtest_discover_tests(test_tim_runtime_guest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES LABELS "unit")

gtest_discover_tests(test_tim_domain
  DISCOVERY_MODE PRE_TEST
  PROPERTIES LABELS "unit")
