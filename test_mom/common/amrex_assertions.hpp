// GoogleTest-style assertions for AMReX types, shared across the MOM C++
// unit tests. Add new helpers here as more kernel tests land.
#pragma once

#include <AMReX_FArrayBox.H>

namespace test_mom {

// Element-wise EXPECT_DOUBLE_EQ between `got` and `expected` over the
// expected array's box. Asserts the boxes match before iterating. Host-side
// -- EXPECT_DOUBLE_EQ is not safe to call from a GPU lambda.
void expect_arrays_equal(const amrex::FArrayBox& expected,
                         const amrex::FArrayBox& got,
                         const char* label);

} // namespace test_mom
