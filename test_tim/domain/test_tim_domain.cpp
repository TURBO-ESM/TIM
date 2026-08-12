// Unit tests for TIM::Domain, the horizontal computational domain and its
// decomposition. A Domain requires a live runtime, hence the main() below,
// which brings one up in owner mode.

#include <gtest/gtest.h>

#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include "core/tim_domain.hpp"
#include "core/tim_runtime.hpp"

namespace {

// The default decomposition is one box per rank (fewer only if the domain is
// too small to split, which this one is not).
TEST(Domain, DefaultDecompositionIsOneBoxPerRank) {
    const TIM::Domain domain({.ni_global = 10, .nj_global = 8,
                              .ni_halo = 2, .nj_halo = 3,
                              .periodic_x = true});
    EXPECT_EQ(domain.n_boxes(), amrex::ParallelDescriptor::NProcs());
}

// A requested box count above one-per-rank is honored.
TEST(Domain, RequestedBoxCountIsHonored) {
    const int n_boxes = 6;
    const TIM::Domain domain({.ni_global = 12, .nj_global = 12, .n_boxes = n_boxes});
    EXPECT_EQ(domain.n_boxes(), n_boxes);
}

// The recorded logical tile grid is consistent with the boxes: the layout
// shape multiplies out to the box count, every box gets in-range tile
// coordinates, and a box touches a global edge exactly when its tile sits on
// the layout boundary (the property the future I/O decomposition's
// tile-based edge ownership relies on).
TEST(Domain, LayoutRecordsTheTileGrid) {
    const int n_boxes = 6;
    const TIM::Domain domain({.ni_global = 12, .nj_global = 12, .n_boxes = n_boxes});
    const amrex::BoxArray boxes = domain.boxArray(1);
    ASSERT_EQ(domain.layout_nx() * domain.layout_ny(), domain.n_boxes());
    for (int b = 0; b < domain.n_boxes(); ++b) {
        const auto [tile_i, tile_j] = domain.tileOf(b);
        ASSERT_GE(tile_i, 0);
        ASSERT_LT(tile_i, domain.layout_nx());
        ASSERT_GE(tile_j, 0);
        ASSERT_LT(tile_j, domain.layout_ny());
        EXPECT_EQ(tile_i == 0, boxes[b].smallEnd(0) == 0);
        EXPECT_EQ(tile_i == domain.layout_nx() - 1,
                  boxes[b].bigEnd(0) == domain.ni_global() - 1);
        EXPECT_EQ(tile_j == 0, boxes[b].smallEnd(1) == 0);
        EXPECT_EQ(tile_j == domain.layout_ny() - 1,
                  boxes[b].bigEnd(1) == domain.nj_global() - 1);
    }
}

// A domain too small for the requested box count yields fewer boxes (the
// empty chunks are dropped), and the recorded tile grid matches what remains.
TEST(Domain, TooSmallDomainClampsBoxCount) {
    const int n_boxes = 7;
    const TIM::Domain domain({.ni_global = 3, .nj_global = 2, .n_boxes = n_boxes});
    const int actual_boxes = domain.n_boxes();
    EXPECT_GT(actual_boxes, 0);
    EXPECT_LT(actual_boxes, n_boxes);
    EXPECT_LE(actual_boxes, 3 * 2);
    EXPECT_EQ(domain.layout_nx() * domain.layout_ny(), actual_boxes);
}

// The periodicity product mirrors the connectivity flags.
TEST(Domain, PeriodicityFollowsTheConnectivityFlags) {
    const TIM::Domain domain({.ni_global = 10, .nj_global = 8,
                              .ni_halo = 2, .nj_halo = 3,
                              .periodic_x = true});
    const amrex::Periodicity periodicity = domain.periodicity();
    EXPECT_TRUE(periodicity.isPeriodic(0));
    EXPECT_FALSE(periodicity.isPeriodic(1));
    EXPECT_FALSE(periodicity.isPeriodic(2));
}

// The Domain's products create working distributed fields: a MultiFab built
// from boxArray()/distribution_mapping() supports a collective reduction
// and a periodic halo exchange driven by periodicity().
TEST(Domain, ProductsCreateWorkingFields) {
    const int ni = 8;
    const int nj = 8;
    const int n_levels = 2;
    const TIM::Domain domain({.ni_global = ni, .nj_global = nj,
                              .ni_halo = 1, .nj_halo = 1,
                              .periodic_x = true, .periodic_y = true});

    // The domain's halo metadata, as the ghost-cell vector (horizontal-only).
    const amrex::IntVect halo = domain.nghost();
    amrex::MultiFab field(domain.boxArray(n_levels),
                          domain.distribution_mapping(), 1, halo);
    constexpr double sentinel = 1.0e30;
    field.setVal(sentinel);
    field.setVal(2.5, 0, 1, /*nghost=*/0);    // interior only
    EXPECT_DOUBLE_EQ(field.sum(), 2.5 * ni * nj * n_levels);

    // Halo exchange: with both directions periodic, every halo cell has a
    // source cell, so no cell in the halo-grown region retains the sentinel.
    field.FillBoundary(domain.periodicity());
    EXPECT_DOUBLE_EQ(field.norminf(0, 1, halo), 2.5);
}

}  // namespace

// Test-binary entry point: initialize GTest, bring up the infrastructure
// layer in owner mode (Runtime calls MPI_Init itself), and run all tests.
// The Runtime is destroyed when main returns, before static teardown.
int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    const TIM::Runtime runtime(argc, argv);
    return RUN_ALL_TESTS();
}
