// Unit tests for TIM::Domain, the horizontal computational domain and its
// decomposition. A Domain requires a live runtime, hence the main() below,
// which brings one up in owner mode.

#include <gtest/gtest.h>

#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_Loop.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include "core/tim_domain.hpp"
#include "core/tim_runtime.hpp"

namespace {

// A valid decomposition at n_levels: the boxes tile the global index space
// exactly once, no box is split in k, and every box is assigned to
// a valid rank.
void expectValidDecomposition(const TIM::Domain& domain, const int n_levels) {
    const amrex::BoxArray box_array = domain.boxArray(n_levels);
    const amrex::Box whole_domain(
        amrex::IntVect(0, 0, 0),
        amrex::IntVect(domain.ni_global() - 1, domain.nj_global() - 1, n_levels - 1));
    EXPECT_EQ(box_array.numPts(), whole_domain.numPts());
    EXPECT_TRUE(box_array.contains(whole_domain));

    const amrex::DistributionMapping& distribution_mapping = domain.distribution_mapping();
    ASSERT_EQ(static_cast<int>(distribution_mapping.size()),
              static_cast<int>(box_array.size()));
    for (int b = 0; b < static_cast<int>(box_array.size()); ++b) {
        EXPECT_EQ(box_array[b].smallEnd(2), 0);
        EXPECT_EQ(box_array[b].bigEnd(2), n_levels - 1);
        EXPECT_GE(distribution_mapping[b], 0);
        EXPECT_LT(distribution_mapping[b], amrex::ParallelDescriptor::NProcs());
    }
}

// The default decomposition is one box per rank (fewer only if the domain is
// too small to split, which this one is not), tiles the domain, and is never
// split in k.
TEST(Domain, DefaultDecompositionIsOneBoxPerRank) {
    const TIM::Domain domain(10, 8, 2, 3, /*periodic_x=*/true,
                             /*periodic_y=*/false, /*tripolar_n=*/false);
    EXPECT_EQ(static_cast<int>(domain.boxArray(3).size()),
              amrex::ParallelDescriptor::NProcs());
    expectValidDecomposition(domain, 3);
}

// A requested box count above one-per-rank is honored, and the round-robin
// distribution mapping stays valid.
TEST(Domain, RequestedBoxCountIsHonored) {
    const int n_boxes = 6;
    const TIM::Domain domain(12, 12, 0, 0, false, false, false, n_boxes);
    EXPECT_EQ(static_cast<int>(domain.boxArray(1).size()), n_boxes);
    EXPECT_EQ(domain.n_boxes(), n_boxes);
    expectValidDecomposition(domain, 1);
}

// The recorded logical tile grid is consistent with the boxes: the layout
// shape multiplies out to the box count, every box gets in-range tile
// coordinates, and a box touches a global edge exactly when its tile sits on
// the layout boundary (the property the future I/O decomposition's
// tile-based edge ownership relies on).
TEST(Domain, LayoutRecordsTheTileGrid) {
    const int n_boxes = 6;
    const TIM::Domain domain(12, 12, 0, 0, false, false, false, n_boxes);
    const amrex::BoxArray boxes = domain.boxArray(1);
    ASSERT_EQ(domain.layout_nx() * domain.layout_ny(),
              static_cast<int>(boxes.size()));
    for (int b = 0; b < static_cast<int>(boxes.size()); ++b) {
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
// empty chunks are dropped), and what remains still tiles the domain as a
// regular tile grid.
TEST(Domain, TooSmallDomainClampsBoxCount) {
    const int n_boxes = 7;
    const TIM::Domain domain(3, 2, 0, 0, false, false, false, n_boxes);
    const int actual_boxes = static_cast<int>(domain.boxArray(1).size());
    EXPECT_GT(actual_boxes, 0);
    EXPECT_LT(actual_boxes, n_boxes);
    EXPECT_LE(actual_boxes, 3 * 2);
    EXPECT_EQ(domain.layout_nx() * domain.layout_ny(), actual_boxes);
    expectValidDecomposition(domain, 1);
}

// The periodicity product mirrors the connectivity flags.
TEST(Domain, PeriodicityFollowsTheConnectivityFlags) {
    const TIM::Domain domain(10, 8, 2, 3, true, false, false);
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
    const TIM::Domain domain(ni, nj, 1, 1, true, true, false);

    // Halos are horizontal-only, matching the domain's halo metadata.
    const amrex::IntVect halo(1, 1, 0);
    amrex::MultiFab field(domain.boxArray(n_levels),
                          domain.distribution_mapping(), 1, halo);
    field.setVal(0.0);                        // everywhere, incl. the halos
    field.setVal(2.5, 0, 1, /*nghost=*/0);    // interior only
    EXPECT_DOUBLE_EQ(field.sum(), 2.5 * ni * nj * n_levels);

    // Halo exchange: with both directions periodic, every halo cell has a
    // source cell, so the grown region becomes uniformly the interior value.
    field.FillBoundary(domain.periodicity());
    amrex::Gpu::streamSynchronize();
    for (amrex::MFIter mfi(field); mfi.isValid(); ++mfi) {
        const amrex::Box grown = amrex::grow(mfi.validbox(), halo);
        const auto a = field.const_array(mfi);
        amrex::LoopOnCpu(grown, [&](int i, int j, int k) {
            EXPECT_DOUBLE_EQ(a(i, j, k), 2.5) << "at (" << i << "," << j << "," << k << ")";
        });
    }
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
