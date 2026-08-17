// Unit tests for TIM::Domain, the horizontal computational domain and its
// decomposition. A Domain requires a live runtime, provided by the shared
// entry point (test_tim_main.cpp), which brings one up in owner mode.

#include <gtest/gtest.h>

#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include "core/tim_domain.hpp"

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

// make_field creates fields on the domain's decomposition with the requested
// staggering expressed as AMReX index-type nodality (one extra plane of
// points in each nodal direction), the requested vertical extent and
// component count, and the domain's halo widths as ghost cells.
TEST(Domain, MakeFieldStaggersAndSizesFields) {
    const int ni = 8;
    const int nj = 6;
    const TIM::Domain domain({.ni_global = ni, .nj_global = nj,
                              .ni_halo = 2, .nj_halo = 1});

    const amrex::MultiFab cell = domain.make_field(TIM::Stagger::Cell, 1, 1);
    const amrex::MultiFab x_face = domain.make_field(TIM::Stagger::XFace, 1, 1);
    const amrex::MultiFab y_face = domain.make_field(TIM::Stagger::YFace, 1, 1);
    const amrex::MultiFab node = domain.make_field(TIM::Stagger::Node, 3, 2);

    // The fields sit on this domain's decomposition: the same boxes on the
    // same ranks, staggering aside (which changes the boxes but not how many
    // there are, nor who owns them).
    EXPECT_EQ(cell.boxArray(), domain.boxArray(1));
    EXPECT_EQ(cell.DistributionMap(), domain.distribution_mapping());
    EXPECT_EQ(node.DistributionMap(), domain.distribution_mapping());
    EXPECT_EQ(node.boxArray().size(), cell.boxArray().size());

    EXPECT_TRUE(cell.ixType().cellCentered());
    EXPECT_EQ(x_face.ixType(), amrex::IndexType(amrex::IntVect(1, 0, 0)));
    EXPECT_EQ(y_face.ixType(), amrex::IndexType(amrex::IntVect(0, 1, 0)));
    EXPECT_EQ(node.ixType(), amrex::IndexType(amrex::IntVect(1, 1, 0)));

    EXPECT_EQ(cell.boxArray().minimalBox().length(0), ni);
    EXPECT_EQ(x_face.boxArray().minimalBox().length(0), ni + 1);
    EXPECT_EQ(y_face.boxArray().minimalBox().length(1), nj + 1);
    EXPECT_EQ(node.boxArray().minimalBox().length(0), ni + 1);
    EXPECT_EQ(node.boxArray().minimalBox().length(2), 3);  // n_levels

    EXPECT_EQ(cell.nComp(), 1);
    EXPECT_EQ(node.nComp(), 2);
    EXPECT_EQ(cell.nGrowVect(), domain.nghost());
    EXPECT_EQ(node.nGrowVect(), domain.nghost());

    // An explicit ghost-cell width overrides the domain's halo widths.
    const amrex::IntVect wide(4, 3, 0);
    const amrex::MultiFab no_halo = domain.make_field(TIM::Stagger::Cell, 1, 1,
                                                      amrex::IntVect(0, 0, 0));
    const amrex::MultiFab wide_halo = domain.make_field(TIM::Stagger::Cell, 1, 1, wide);
    EXPECT_EQ(no_halo.nGrowVect(), amrex::IntVect(0, 0, 0));
    EXPECT_EQ(wide_halo.nGrowVect(), wide);
}

// A staggered field does not tile its index space: the extra plane per nodal
// direction means the boxes on either side of an internal boundary each hold
// a valid copy of the plane between them. The duplicates are stored, and
// whole-field reductions see them.
TEST(Domain, MakeFieldStaggeredBoxesShareBoundaryPlanes) {
    const int ni = 8;
    const int nj = 6;
    const TIM::Domain domain({.ni_global = ni, .nj_global = nj, .n_boxes = 4});
    ASSERT_EQ(domain.n_boxes(), 4);

    const amrex::MultiFab cell = domain.make_field(TIM::Stagger::Cell, 1, 1);
    const amrex::MultiFab node = domain.make_field(TIM::Stagger::Node, 1, 1);

    // A cell-centered decomposition tiles exactly: the boxes' point counts add
    // up to the distinct points they cover.
    EXPECT_EQ(cell.boxArray().numPts(), cell.boxArray().minimalBox().numPts());
    EXPECT_EQ(cell.boxArray().minimalBox().numPts(), ni * nj);

    // A node field covers (ni+1) x (nj+1) distinct points, but its boxes hold
    // more than that between them: the shared planes are stored twice.
    EXPECT_EQ(node.boxArray().minimalBox().numPts(), (ni + 1) * (nj + 1));
    EXPECT_GT(node.boxArray().numPts(), node.boxArray().minimalBox().numPts());

    // So a reduction counts every stored point, shared planes included, rather
    // than every distinct point.
    amrex::MultiFab ones = domain.make_field(TIM::Stagger::Node, 1, 1);
    ones.setVal(1.0);
    EXPECT_DOUBLE_EQ(ones.sum(), static_cast<double>(node.boxArray().numPts()));
}

}  // namespace
