// Unit tests for TIM::checksum: the Array4 entry points exercise the
// checksum core directly, and the MultiFab entry points exercise it once per
// local box, so together these cover the shared core through both surfaces.
// A live runtime (for Domain) is provided by the shared entry point
// (test_tim_main.cpp).

#include <cstring>

#include <gtest/gtest.h>

#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_MultiFab.H>

#include "core/tim_domain.hpp"
#include "tim_coms_infra.hpp"

namespace {

// Bitwise reinterpretation of a double as its 64-bit representation,
// mirroring the cast TIM::checksum itself performs; computed here
// independently to build the expected values below.
amrex::Long bits_of(double value) {
    amrex::Long bits;
    std::memcpy(&bits, &value, sizeof(bits));
    return bits;
}

// A field on a single-box domain, so tests can reach into its one Array4
// and their expected values stay independent of how many ranks ctest
// happens to run under (every other rank simply owns no boxes).
amrex::MultiFab make_single_box_field(int ni, int nj, int halo = 0) {
    const TIM::Domain domain({.ni_global = ni, .nj_global = nj,
                              .ni_halo = halo, .nj_halo = halo,
                              .periodic_x = true, .periodic_y = true,
                              .n_boxes = 1});
    return domain.make_field({.stagger = TIM::Stagger::Cell, .nk = 1});
}

TEST(Checksum, Array4MatchesManualBitSum) {
    amrex::MultiFab field = make_single_box_field(4, 3);
    field.setVal(2.5);

    amrex::MFIter mfi(field);
    ASSERT_TRUE(mfi.isValid());
    const amrex::Box bx = mfi.validbox();
    const amrex::Long expected = bits_of(2.5) * bx.numPts();
    EXPECT_EQ(TIM::checksum(bx, field.const_array(mfi)), expected);
}

TEST(Checksum, Array4MaskExcludesMatchingElements) {
    constexpr double mask_value = -999.0;
    constexpr double data_value = 2.5;
    constexpr int n_unmasked = 3;

    amrex::MultiFab field = make_single_box_field(4, 3);
    field.setVal(mask_value);

    amrex::MFIter mfi(field);
    ASSERT_TRUE(mfi.isValid());
    const amrex::Box bx = mfi.validbox();
    amrex::Array4<amrex::Real> const arr = field.array(mfi);
    amrex::ParallelFor(n_unmasked, [=] AMREX_GPU_DEVICE (int n) {
        arr(n, 0, 0) = data_value;
    });

    const amrex::Long expected = bits_of(data_value) * n_unmasked;
    EXPECT_EQ(TIM::checksum(bx, field.const_array(mfi), mask_value), expected);
}

// Rank count-independent: every valid cell in the global field contributes
// exactly once, regardless of how the domain's default decomposition splits
// it across boxes/ranks.
TEST(Checksum, MultiFabMatchesTotalCellCount) {
    constexpr double data_value = 1.25;
    const TIM::Domain domain({.ni_global = 6, .nj_global = 5});
    amrex::MultiFab field = domain.make_field({.stagger = TIM::Stagger::Cell, .nk = 2});
    field.setVal(data_value);

    const amrex::Long expected = bits_of(data_value) * 6 * 5 * 2;
    EXPECT_EQ(TIM::checksum(field), expected);
}

TEST(Checksum, MultiFabDefaultExcludesGhostCells) {
    constexpr double sentinel = 1.0e30;
    constexpr double data_value = 2.5;

    amrex::MultiFab field = make_single_box_field(4, 3, /*halo=*/1);
    field.setVal(sentinel);
    field.setVal(data_value, 0, 1, /*nghost=*/0);  // interior only

    const amrex::Long expected = bits_of(data_value) * 4 * 3;
    EXPECT_EQ(TIM::checksum(field), expected);
}

TEST(Checksum, MultiFabNghostIncludesGhostCells) {
    constexpr double data_value = 2.5;
    constexpr int ni = 4;
    constexpr int nj = 3;
    constexpr int halo = 1;

    const TIM::Domain domain({.ni_global = ni, .nj_global = nj,
                              .ni_halo = halo, .nj_halo = halo,
                              .periodic_x = true, .periodic_y = true,
                              .n_boxes = 1});
    amrex::MultiFab field = domain.make_field({.stagger = TIM::Stagger::Cell, .nk = 1});
    field.setVal(data_value);
    // Periodic on a single box: the halo exchange fills every ghost cell
    // from the (uniform) valid region, so it carries data_value too.
    field.FillBoundary(domain.periodicity());

    const amrex::Long expected =
        bits_of(data_value) * (ni + 2 * halo) * (nj + 2 * halo);
    EXPECT_EQ(TIM::checksum(field, domain.nghost()), expected);
}

TEST(Checksum, MultiFabMaskExcludesMatchingElements) {
    constexpr double mask_value = -999.0;
    constexpr double data_value = 2.5;
    constexpr int n_unmasked = 3;

    amrex::MultiFab field = make_single_box_field(4, 3);
    field.setVal(mask_value);

    amrex::MFIter mfi(field);
    ASSERT_TRUE(mfi.isValid());
    amrex::Array4<amrex::Real> const arr = field.array(mfi);
    amrex::ParallelFor(n_unmasked, [=] AMREX_GPU_DEVICE (int n) {
        arr(n, 0, 0) = data_value;
    });

    const amrex::Long expected = bits_of(data_value) * n_unmasked;
    EXPECT_EQ(TIM::checksum(field, mask_value), expected);
}

}  // namespace
