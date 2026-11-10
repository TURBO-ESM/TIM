# TIM C++ unit tests
#
# Run from the TIM root, with parallelio (PIO) and AMReX installs on CMAKE_PREFIX_PATH:
#   ./build.sh --tests                              # build + install + ctest
# or with raw CMake:
#   cmake -S . -B build -DTIM_ADD_TEST_TARGETS=ON \
#         -DCMAKE_PREFIX_PATH="<amrex-install>;$NCAR_ROOT_PARALLELIO"
#   cmake --build build
#   ctest --test-dir build -L unit --output-on-failure
# or standalone -- just tim/cpp and its tests, without the full TIM/FMS build
#   cmake -S test_tim -B build -DCMAKE_PREFIX_PATH=<amrex-install>
#   cmake --build build && ctest --test-dir build -L unit
#
# Labels: `unit` runs every test case serially, one ctest entry per case;
# `mpi<N>` runs whole test executables on N MPI ranks (e.g. ctest -L mpi2).

# Standalone configuration: provide what the top-level build provides
# otherwise. AMReX comes from CMAKE_PREFIX_PATH. tim/cpp resolves PIO for the
# io targets: an install on CMAKE_PREFIX_PATH / $PIO, or a libpioc source build.
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  cmake_minimum_required(VERSION 3.24)
  project(TIM_cpp_tests LANGUAGES C CXX)
  set(TIM_ENABLE_IO ON)
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

# tim_add_mpi_test(<target> RANKS <n>...)
#   Run the whole test executable under the MPI launcher once per rank count,
#   as ctest entry <target>_np<n> with label mpi<n>.
function(tim_add_mpi_test target)
  cmake_parse_arguments(PARSE_ARGV 1 arg "" "" "RANKS")
  foreach(n IN LISTS arg_RANKS)
    add_test(NAME ${target}_np${n}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${n} ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:${target}> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(${target}_np${n} PROPERTIES LABELS "mpi${n}" PROCESSORS ${n})
  endforeach()
endfunction()

if(AMReX_GPU_BACKEND STREQUAL "CUDA")
  enable_language(CUDA)
  set_source_files_properties(
    test_tim_main.cpp
    smoke/test_smoke.cpp
    runtime/test_runtime_owner.cpp
    runtime/test_runtime_guest.cpp
    domain/test_tim_domain.cpp
    time/test_tim_time.cpp
    coms_infra/test_tim_checksum.cpp
    PROPERTIES LANGUAGE CUDA)
endif()

add_executable(test_tim_smoke
  test_tim_main.cpp
  smoke/test_smoke.cpp
)

add_executable(test_tim_io_smoke
  io/test_tim_io_smoke.cpp
)

add_executable(test_tim_runtime_owner
  runtime/test_runtime_owner.cpp
)

add_executable(test_tim_runtime_guest
  runtime/test_runtime_guest.cpp
)

add_executable(test_tim_domain
  test_tim_main.cpp
  domain/test_tim_domain.cpp
)

add_executable(test_tim_time
  time/test_tim_time.cpp
)

add_executable(test_tim_checksum
  test_tim_main.cpp
  coms_infra/test_tim_checksum.cpp
)

target_link_libraries(test_tim_smoke
  GTest::gtest
  tim
)

# The AMReX-free layers link nothing above them: test_tim_io_smoke and
# test_tim_time are built without AMReX on their link or include lines.
target_link_libraries(test_tim_io_smoke
  GTest::gtest
  tim_io
)

foreach(tim_test test_tim_runtime_owner test_tim_runtime_guest test_tim_domain)
  target_link_libraries(${tim_test}
    GTest::gtest
    tim
  )
endforeach()

target_link_libraries(test_tim_time
  GTest::gtest_main
  tim_base
)

target_link_libraries(test_tim_checksum
  GTest::gtest
  tim
)

include(GoogleTest)

gtest_discover_tests(test_tim_smoke
  DISCOVERY_MODE PRE_TEST
  PROPERTIES LABELS "unit")

gtest_discover_tests(test_tim_io_smoke
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

gtest_discover_tests(test_tim_time
  DISCOVERY_MODE PRE_TEST
  PROPERTIES LABELS "unit")

gtest_discover_tests(test_tim_checksum
  DISCOVERY_MODE PRE_TEST
  PROPERTIES LABELS "unit")

foreach(tim_test test_tim_io_smoke test_tim_runtime_owner test_tim_runtime_guest
                 test_tim_domain)
  tim_add_mpi_test(${tim_test} RANKS 2)
endforeach()
