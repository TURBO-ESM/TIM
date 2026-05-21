#include "amrex_assertions.hpp"

#include <gtest/gtest.h>

#include <AMReX_Arena.H>
#include <AMReX_Gpu.H>

#include <cstddef>

namespace test_mom {

void expect_arrays_equal(const amrex::FArrayBox& expected,
                         const amrex::FArrayBox& got,
                         const char* label) {

    static_assert(AMREX_SPACEDIM == 3, "expect_arrays_equal assumes 3D arrays");

    ASSERT_TRUE(got.box() == expected.box()) << label << ": box mismatch";
    const auto e  = expected.const_array();
    const auto g  = got.const_array();
    const auto lo = expected.box().smallEnd();
    const auto hi = expected.box().bigEnd();
    for (int k = lo[2]; k <= hi[2]; ++k)
        for (int j = lo[1]; j <= hi[1]; ++j)
            for (int i = lo[0]; i <= hi[0]; ++i)
                EXPECT_EQ(g(i, j, k), e(i, j, k))
                    << label << " mismatch at (" << i << "," << j << "," << k << ")";
}

amrex::FArrayBox to_host_fab(const amrex::FArrayBox& src) {
    amrex::FArrayBox dst(src.box(), src.nComp(), amrex::The_Pinned_Arena());
    const std::size_t n = static_cast<std::size_t>(src.box().numPts()) *
                          static_cast<std::size_t>(src.nComp());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     src.dataPtr(), src.dataPtr() + n,
                     dst.dataPtr());
    amrex::Gpu::streamSynchronize();
    return dst;
}

} // namespace test_mom
