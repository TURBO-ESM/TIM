# MOM C++ kernels library.
#
# Builds the AMReX-based MOM kernels as a static library

if(NOT TARGET AMReX::amrex_3d AND NOT TARGET amrex_3d)
    find_package(AMReX REQUIRED)
endif()

add_library(mom_continuity STATIC
    mom_continuity_ppm.cpp
)

if(AMReX_GPU_BACKEND STREQUAL "CUDA")
    set_source_files_properties(mom_continuity_ppm.cpp PROPERTIES LANGUAGE CUDA)
endif()

target_include_directories(mom_continuity PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(mom_continuity PUBLIC AMReX::amrex_3d)

target_compile_features(mom_continuity PUBLIC cxx_std_17)

# To ensure bit-identical results with the Fortran reference, we need to 
# disable FMA and reassociation.
target_compile_options(mom_continuity PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-ffp-contract=off>
    $<$<CXX_COMPILER_ID:IntelLLVM,Intel>:-fp-model=precise>
)
