# TIM C++ unit tests
#
# Run from the TIM root, with parallelio (PIO) and AMReX installs on CMAKE_PREFIX_PATH:
#   ./build.sh --tests                              # build + install + ctest
# or with raw CMake:
#   cmake -S . -B build -DTIM_ADD_TEST_TARGETS=ON \
#         -DCMAKE_PREFIX_PATH="<amrex-install>;$NCAR_ROOT_PARALLELIO"
#   cmake --build build --target test_tim_smoke
#   ctest --test-dir build -L unit --output-on-failure

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

add_executable(test_tim_smoke
  test_tim_main.cpp
  smoke/test_smoke.cpp
)

target_link_libraries(test_tim_smoke
  GTest::gtest
  tim
  PIO::PIO_C
)

target_link_options(test_tim_smoke PRIVATE LINKER:--allow-shlib-undefined)

if(AMReX_GPU_BACKEND STREQUAL "CUDA")
  set_source_files_properties(
    test_tim_main.cpp
    smoke/test_smoke.cpp
    PROPERTIES LANGUAGE CUDA)
endif()

include(GoogleTest)

gtest_discover_tests(test_tim_smoke
  DISCOVERY_MODE PRE_TEST
  PROPERTIES LABELS "unit")
