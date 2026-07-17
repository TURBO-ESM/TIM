# ParallelIO (PIO2) dependency. Prefer a prebuilt install (e.g. the
# Derecho parallelio module, which ships PIOConfig.cmake on CMAKE_PREFIX_PATH
# and defines the imported target PIO::PIO_C). The PIO env variable is also
# honored as a hint (set by the parallelio module and the CI container image).
# When none is found, fall back to building libpioc from source via FetchContent.

find_package(PIO QUIET COMPONENTS C HINTS $ENV{PIO})

if(NOT PIO_FOUND)
  if(TIM_FETCH_DEPS)
    message(STATUS "PIO not found; building libpioc from source via FetchContent (pio2_6_8)")
    include(FetchContent)
    FetchContent_Declare(
      parallelio
      URL https://github.com/NCAR/ParallelIO/archive/refs/tags/pio2_6_8.tar.gz
      DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    # C library only: TIM consumes just the PIOc_* C API.
    set(PIO_ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)
    set(PIO_ENABLE_TIMING  OFF CACHE BOOL "" FORCE)
    set(PIO_ENABLE_DOC     OFF CACHE BOOL "" FORCE)
    set(PIO_ENABLE_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(PIO_ENABLE_TESTS   OFF CACHE BOOL "" FORCE)
    set(WITH_PNETCDF       OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(parallelio)
    # Two embedded-build fixups (PIO's CMake assumes it is the top-level
    # project, built with MPI compiler wrappers), applied here instead of
    # patching PIO:
    #  - config.h is generated into PIO's PROJECT_BINARY_DIR but the include
    #    path carries the top-level CMAKE_BINARY_DIR; point pioc at the former.
    #  - pio.h includes mpi.h but pioc never links MPI::MPI_C (standalone
    #    builds get MPI from the mpicc wrapper); attach it explicitly. PUBLIC
    #    so pio.h consumers inherit the MPI include path too.
    target_include_directories(pioc PRIVATE ${parallelio_BINARY_DIR})
    cmake_policy(SET CMP0079 NEW) # allow linking a target defined in the PIO subdirectory
    target_link_libraries(pioc PUBLIC MPI::MPI_C)
    add_library(PIO::PIO_C ALIAS pioc)
  else()
    message(FATAL_ERROR
      "PIO not found. Load the parallelio module / add a PIO install to "
      "CMAKE_PREFIX_PATH, or configure with -DTIM_FETCH_DEPS=ON to build "
      "libpioc from source automatically.")
  endif()
endif()
