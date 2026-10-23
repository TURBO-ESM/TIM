// GoogleTest-style assertions for AMReX types, shared across the MOM C++
// unit tests. Add new helpers here as more kernel tests land.
// SKILLS: 0.3.3
#pragma once

#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>

namespace test_mom {

// Element-wise EXPECT_DOUBLE_EQ between `got` and `expected` over the
// expected array's box. Asserts the boxes match before iterating.
//
// Both FABs must be host-resident (e.g. from CapturedFile::fab_host, or
// produced by to_host_fab()).
void expect_arrays_equal(const amrex::FArrayBox& expected,
                         const amrex::FArrayBox& got,
                         const char* label);

// Allocate a pinned-host FArrayBox with the same shape as `src` and copy
// `src` into it. Works whether `src` lives on host or device; on CPU builds
// the underlying Gpu::copy is a std::copy. Use this to bring a kernel-output
// device FAB back to host before calling expect_arrays_equal.
amrex::FArrayBox to_host_fab(const amrex::FArrayBox& src);

// Element-wise EXPECT_EQ between `got` and `expected` integer arrays (e.g. the
// integer-encoded logicals from CapturedFile::int_fab_host) over the expected
// array's box. Asserts the boxes match before iterating. Both IArrayBoxes must
// be host-resident.
void expect_int_arrays_equal(const amrex::IArrayBox& expected,
                             const amrex::IArrayBox& got,
                             const char* label);

// IArrayBox counterpart of to_host_fab().
amrex::IArrayBox to_host_int_fab(const amrex::IArrayBox& src);

} // namespace test_mom
