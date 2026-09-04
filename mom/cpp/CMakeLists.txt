# MOM C++ kernels, and the Fortran<->C++ bridge MOM6 calls under -D_TIM.
#
#   mom_continuity  the AMReX kernels alone, for test_mom. EXCLUDE_FROM_ALL.
#   mom_bridge      kernels + bridge, exported as TIM::mom_bridge. Built only
#                   when the top level sets TIM_ENABLE_MOM_BRIDGE.

if(NOT TARGET AMReX::amrex_3d AND NOT TARGET amrex_3d)
    find_package(AMReX REQUIRED)
endif()

set(mom_turbotmp_dir ${CMAKE_CURRENT_SOURCE_DIR}/../../turbotmp)

# Settings shared by every target built from these sources. Factored into one
# place because of the FP-precision options: a target that misses them still
# compiles and links, it just stops being bit-identical to Fortran, which is the
# one property this directory exists to guarantee.
#
# Include directories are deliberately not set here -- mom_bridge is exported and
# so needs BUILD_INTERFACE/INSTALL_INTERFACE forms, which install(EXPORT) rejects
# a bare source path for.
function(mom_cpp_configure tgt)
    target_link_libraries(${tgt} PUBLIC AMReX::amrex_3d)

    target_compile_features(${tgt} PUBLIC cxx_std_17)

    # To ensure bit-identical results with the Fortran reference, we need to
    # disable FMA and reassociation.
    target_compile_options(${tgt} PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang>:-ffp-contract=off>
        $<$<CXX_COMPILER_ID:IntelLLVM,Intel>:-fp-model=precise>
    )
endfunction()

# The kernels alone, for test_mom. EXCLUDE_FROM_ALL because a full TIM build has
# no use for it -- test_mom brings this directory in through its own
# add_subdirectory(), and linking the target there is enough to get it built.
add_library(mom_continuity STATIC EXCLUDE_FROM_ALL
    mom_continuity_ppm.cpp
)
mom_cpp_configure(mom_continuity)
target_include_directories(mom_continuity PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# The extern "C" bridge MOM6 calls from its `#ifdef _TIM` branches in
# src/core/MOM_continuity_PPM.F90, exported as TIM::mom_bridge.
#
# Defined only when the top level asks for it, so that test_mom -- which
# configures this directory standalone and never sets the option -- does not
# build or install it.
#
# Carries mom_continuity_ppm.cpp rather than linking mom_continuity: a consumer
# resolves the bridge from one archive, and mom_continuity stays a test-only
# target that never has to join TIMExports.
if(TIM_ENABLE_MOM_BRIDGE)
    add_library(mom_bridge STATIC
        mom_continuity_ppm.cpp
        turbotmp_mom_continuity_ppm_bridge.cpp
        ${mom_turbotmp_dir}/turbotmp_helper.cpp
    )
    add_library(TIM::mom_bridge ALIAS mom_bridge)
    mom_cpp_configure(mom_bridge)

    target_include_directories(mom_bridge PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${mom_turbotmp_dir}>
        $<INSTALL_INTERFACE:include/tim/mom>)

    install(TARGETS mom_bridge
        EXPORT TIMExports
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
    # The bridge's public surface, and only that.  MOM6 declares the entry points
    # with Fortran `interface` blocks and includes nothing, but a C or C++ consumer
    # needs these two -- and needs no more, because
    # turbotmp_mom_continuity_ppm_bridge.h includes only turbotmp_bridge_c_types.h,
    # which in turn includes only <stdint.h> and <stddef.h>.
    #
    # turbotmp_helper.hpp is deliberately NOT installed: it is an implementation
    # header (AMReX Array4/MultiFab helpers) reached only from the .cpp files, so
    # installing it would put AMReX includes on the public path for no consumer.
    # It stays visible to the build through the BUILD_INTERFACE entry above.
    install(FILES
        turbotmp_mom_continuity_ppm_bridge.h
        turbotmp_bridge_c_types.h
        DESTINATION include/tim/mom)
endif()

# Every source here includes AMReX headers. With the CUDA backend those headers
# use __host__/__device__ qualifiers and CUDA kernel-launch syntax, so nvcc must
# compile them; CMake will not infer that from a .cpp extension. Source-file
# properties are directory-scoped, so this covers both targets above.
if(AMReX_GPU_BACKEND STREQUAL "CUDA")
    set_source_files_properties(
        mom_continuity_ppm.cpp
        turbotmp_mom_continuity_ppm_bridge.cpp
        ${mom_turbotmp_dir}/turbotmp_helper.cpp
        PROPERTIES LANGUAGE CUDA)
endif()
