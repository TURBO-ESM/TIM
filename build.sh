#!/bin/bash
# Usage: ./build.sh [options]
#
# Builds TIM using CMake. Assumes dependencies (MPI, NetCDF, AMReX) are
# already on CMAKE_PREFIX_PATH (e.g. via an active Spack environment).
#
# Optional environment variables:
#   TIM_ROOT    Path to the TIM repository clone (defaults to this script's
#               directory, so it normally need not be set)
#
# Options:
#   --clean         Full clean rebuild (deletes build dir and rebuilds from scratch)
#   --debug         Build in Debug mode (default: Release)
#   --ninja         Use Ninja generator instead of the default (Unix Makefiles)
#   --build-dir DIR Override build directory (default: $TIM_ROOT/build)
#   --prefix DIR    Install prefix (default: $build_dir/install)
#   --parallel N    Build with N parallel jobs
#   --tests         Build the C++ unit tests (test_tim) and run ctest. Requires
#                   AMReX and PIO to be on CMAKE_PREFIX_PATH.

set -eo pipefail

TIM_ROOT="${TIM_ROOT:-$(cd -P "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

# Sanity-check the deduced/overridden root really is a TIM tree
if [[ ! -f "$TIM_ROOT/CMakeLists.txt" ]]; then
    echo "Error: TIM_ROOT ('$TIM_ROOT') doesn't look like a TIM checkout." >&2
    exit 1
fi

# Defaults
build_dir="$TIM_ROOT/build"
build_type="Release"
debug=false
clean=false
ninja=false
parallel_jobs=""
prefix=""
run_tests=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --debug)        debug=true; shift ;;
        --clean)        clean=true; shift ;;
        --ninja)        ninja=true; shift ;;
        --build-dir)    build_dir="$2"; shift 2 ;;
        --prefix)       prefix="$2"; shift 2 ;;
        --parallel)     parallel_jobs="$2"; shift 2 ;;
        --tests)        run_tests=true; shift ;;
        *)              echo "Error: unknown argument '$1'" >&2; exit 1 ;;
    esac
done

# Derive the install prefix from the (possibly overridden) build dir unless
# the user set one explicitly.
prefix="${prefix:-$build_dir/install}"

if [[ "$debug" == true ]]; then
    build_type="Debug"
fi

if [[ "$clean" == true ]]; then
    rm -rf "$build_dir"
fi

# Configure
cmake_configure_opts=(
    -S "$TIM_ROOT"
    -B "$build_dir"
    -DCMAKE_BUILD_TYPE="$build_type"
    -DCMAKE_INSTALL_PREFIX="$prefix"
)
[[ "$ninja" == true ]] && cmake_configure_opts+=(-G Ninja)
[[ "$run_tests" == true ]] && cmake_configure_opts+=(-DTIM_ADD_TEST_TARGETS=ON)
cmake "${cmake_configure_opts[@]}"

# Build
cmake_build_opts=("$build_dir")
[[ -n "$parallel_jobs" ]] && cmake_build_opts+=(--parallel "$parallel_jobs")
cmake --build "${cmake_build_opts[@]}"

# Install into $prefix so the build produces a consumable tree
# (TIMConfig.cmake, libs, modules) like the other turbo deps.
cmake --install "$build_dir"

# Run the C++ unit tests if requested.
if [[ "$run_tests" == true ]]; then
    ctest_opts=(--test-dir "$build_dir" --output-on-failure)
    [[ -n "$parallel_jobs" ]] && ctest_opts+=(-j "$parallel_jobs")
    ctest "${ctest_opts[@]}"
fi
