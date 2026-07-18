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
    find_package(MPI REQUIRED COMPONENTS C)
    FetchContent_MakeAvailable(parallelio)
    # Two embedded-build fixups, scoped to pioc:
    #  - pio.h includes mpi.h, but PIO's CMake assumes an MPI compiler wrapper and
    #    puts no MPI usage requirement on pioc; give it MPI's include path so it
    #    also compiles under a plain compiler.
    #  - PIO generates config.h into its own PROJECT_BINARY_DIR, but its include
    #    path carries the top-level CMAKE_BINARY_DIR (identical only when PIO is
    #    the top-level project); point pioc at the directory that holds it.
    target_include_directories(pioc SYSTEM PRIVATE ${MPI_C_INCLUDE_DIRS})
    target_include_directories(pioc PRIVATE ${parallelio_BINARY_DIR})
    add_library(PIO::PIO_C ALIAS pioc)
  else()
    message(FATAL_ERROR
      "PIO not found. Load the parallelio module / add a PIO install to "
      "CMAKE_PREFIX_PATH, or configure with -DTIM_FETCH_DEPS=ON to build "
      "libpioc from source automatically.")
  endif()
endif()
