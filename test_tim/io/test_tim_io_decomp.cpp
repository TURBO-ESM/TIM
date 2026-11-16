/**
 * @file test_tim_io_decomp.cpp
 * @brief Unit tests for TIM::IoDecomp: window sizes per staggering and the
 * exactly-once write partition, including random masks
 * and compute ranks that hold several tiles.
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "io/tim_io_decomp.hpp"

using TIM::IndexBox;
using TIM::IoDecomp;
using TIM::IoDecompSpec;
using TIM::Stagger;

namespace {

using Tile = std::array<int, 2>;

constexpr std::array<Stagger, 4> kStaggers{Stagger::Cell, Stagger::XFace, Stagger::YFace,
                                           Stagger::Node};

// A box from its inclusive i range and j range.
IndexBox ij_box(int ilo, int ihi, int jlo, int jhi) { return IndexBox{{ilo, jlo}, {ihi, jhi}}; }

std::string name(Stagger s) {
    switch (s) {
        case Stagger::Cell: return "Cell";
        case Stagger::XFace: return "XFace";
        case Stagger::YFace: return "YFace";
        case Stagger::Node: return "Node";
    }
    return "?";
}

testing::AssertionResult writes_exactly(const IoDecomp& decomp, int b, Stagger s,
                                        const std::vector<IndexBox>& boxes) {
    const IndexBox window = decomp.window(b, s);
    const std::vector<bool> mask = decomp.write_mask(b, s);
    if (mask.size() != static_cast<std::size_t>(window.n_points())) {
        return testing::AssertionFailure() << "the mask has " << mask.size() << " entries";
    }
    std::size_t k = 0;
    for (int j = window.lo[1]; j <= window.hi[1]; ++j) {
        for (int i = window.lo[0]; i <= window.hi[0]; ++i, ++k) {
            const bool want = std::any_of(boxes.begin(), boxes.end(),
                                          [i, j](const IndexBox& box) { return box.contains(i, j); });
            if (mask[k] != want) {
                return testing::AssertionFailure() << name(s) << " tile " << b << " index (" << i << ", "
                                                   << j << "): mask is " << mask[k];
            }
        }
    }
    return testing::AssertionSuccess();
}

// A 2 x 2 layout: tile columns of 3 and 4 cells, tile rows of 2 and 3 cells.
//   tile (0, 0): i 0..2, j 0..1     tile (1, 0): i 3..6, j 0..1
//   tile (0, 1): i 0..2, j 2..4     tile (1, 1): i 3..6, j 2..4
IoDecompSpec two_by_two(bool symmetric, std::vector<Tile> local,
                        std::vector<Tile> eliminated = {}) {
    return {.tile_ni = {3, 4},
            .tile_nj = {2, 3},
            .symmetric = symmetric,
            .eliminated_tiles = std::move(eliminated),
            .local_tiles = std::move(local)};
}

// A random layout and the test's own model of the ownership rule. Up to 5 x 5
// tiles of 1 to 4 cells each, about a third of them eliminated, with the live
// tiles dealt out to 1 to 4 compute ranks in random order, so that compute
// ranks hold any number of tiles (including none).
struct Layout {
    std::vector<int> tile_ni;
    std::vector<int> tile_nj;
    std::vector<Tile> eliminated;
    std::vector<std::vector<Tile>> compute_ranks;  // the local tiles of each compute rank
    std::vector<int> column_of;                    // the tile column of each cell column
    std::vector<int> row_of;                       // the tile row of each cell row

    explicit Layout(std::mt19937& rng) {
        std::uniform_int_distribution<int> n_tiles(1, 5);
        std::uniform_int_distribution<int> width(1, 4);
        std::uniform_int_distribution<int> n_compute_ranks(1, 4);
        std::bernoulli_distribution eliminate(0.35);
        tile_ni.resize(static_cast<std::size_t>(n_tiles(rng)));
        tile_nj.resize(static_cast<std::size_t>(n_tiles(rng)));
        for (std::vector<int>* widths : {&tile_ni, &tile_nj}) {
            for (int& w : *widths) {
                w = width(rng);
            }
        }
        for (int t = 0; t < static_cast<int>(tile_ni.size()); ++t) {
            column_of.insert(column_of.end(), static_cast<std::size_t>(tile_ni[static_cast<std::size_t>(t)]), t);
        }
        for (int t = 0; t < static_cast<int>(tile_nj.size()); ++t) {
            row_of.insert(row_of.end(), static_cast<std::size_t>(tile_nj[static_cast<std::size_t>(t)]), t);
        }
        std::vector<Tile> live_tiles;
        for (int tj = 0; tj < static_cast<int>(tile_nj.size()); ++tj) {
            for (int ti = 0; ti < static_cast<int>(tile_ni.size()); ++ti) {
                (eliminate(rng) ? eliminated : live_tiles).push_back({ti, tj});
            }
        }
        std::shuffle(live_tiles.begin(), live_tiles.end(), rng);
        compute_ranks.resize(static_cast<std::size_t>(n_compute_ranks(rng)));
        std::uniform_int_distribution<std::size_t> pick(0, compute_ranks.size() - 1);
        for (const Tile& t : live_tiles) {
            compute_ranks[pick(rng)].push_back(t);
        }
    }

    bool live(const Tile& t) const {
        return std::find(eliminated.begin(), eliminated.end(), t) == eliminated.end();
    }

    // The tile that writes file index (p, q), or nullopt if no tile does: the
    // first live cell touching the index, in the order (p, q), (p - 1, q),
    // (p, q - 1), (p - 1, q - 1).
    std::optional<Tile> owner(bool symmetric, Stagger s, int p, int q) const {
        const bool two_i = symmetric && (s == Stagger::XFace || s == Stagger::Node);
        const bool two_j = symmetric && (s == Stagger::YFace || s == Stagger::Node);
        const int ni = static_cast<int>(column_of.size());
        const int nj = static_cast<int>(row_of.size());
        for (const auto& [di, dj] : {Tile{0, 0}, Tile{1, 0}, Tile{0, 1}, Tile{1, 1}}) {
            const int i = p - di;
            const int j = q - dj;
            if ((di == 1 && !two_i) || (dj == 1 && !two_j) || i < 0 || i >= ni || j < 0 || j >= nj) {
                continue;
            }
            const Tile t{column_of[static_cast<std::size_t>(i)], row_of[static_cast<std::size_t>(j)]};
            if (live(t)) {
                return t;
            }
        }
        return std::nullopt;
    }
};

}  // namespace

// ---- IndexBox ----

TEST(IndexBox, CountsPoints) {
    // The default box is empty.
    EXPECT_TRUE(IndexBox{}.empty());
    EXPECT_EQ(IndexBox{}.n_points(), 0);
    const IndexBox box = ij_box(2, 4, 3, 7);
    EXPECT_FALSE(box.empty());
    EXPECT_EQ(box.length_i(), 3);
    EXPECT_EQ(box.length_j(), 5);
    EXPECT_EQ(box.n_points(), 15);
    // Empty in one direction or both: each length is 0, so two negative raw
    // lengths never multiply into a positive count.
    EXPECT_EQ(ij_box(0, -1, 5, 9).n_points(), 0);
    const IndexBox both = ij_box(5, 2, 5, 2);
    EXPECT_TRUE(both.empty());
    EXPECT_EQ(both.length_i(), 0);
    EXPECT_EQ(both.length_j(), 0);
    EXPECT_EQ(both.n_points(), 0);
    // Negative indices are ordinary bounds, and the count is 64-bit.
    EXPECT_EQ(ij_box(-5, -2, -3, -1).n_points(), 12);
    EXPECT_EQ(ij_box(0, 99999, 0, 99999).n_points(), std::int64_t{100000} * 100000);
}

TEST(IndexBox, Contains) {
    EXPECT_FALSE(IndexBox{}.contains(0, 0));
    const IndexBox box = ij_box(2, 4, 3, 7);
    EXPECT_TRUE(box.contains(2, 3));
    EXPECT_TRUE(box.contains(4, 7));
    EXPECT_FALSE(box.contains(1, 3));
    EXPECT_FALSE(box.contains(4, 8));
    EXPECT_TRUE(ij_box(-5, -2, -3, -1).contains(-5, -3));
    EXPECT_FALSE(ij_box(-5, -2, -3, -1).contains(0, 0));
}

// ---- sizes and windows ----

TEST(IoDecomp, GlobalSizesPerStagger) {
    const IoDecomp plain(two_by_two(false, {{0, 0}}));
    const IoDecomp symmetric(two_by_two(true, {{0, 0}}));
    EXPECT_EQ(plain.layout_nx(), 2);
    EXPECT_EQ(plain.layout_ny(), 2);
    EXPECT_FALSE(plain.symmetric());
    EXPECT_TRUE(symmetric.symmetric());
    // A non-symmetric file is 7 x 5 for every staggering. A symmetric file
    // adds one index along each staggered direction.
    const std::array<std::tuple<Stagger, int, int>, 4> symmetric_sizes{
        {{Stagger::Cell, 7, 5}, {Stagger::XFace, 8, 5}, {Stagger::YFace, 7, 6}, {Stagger::Node, 8, 6}}};
    for (const auto& [s, ni, nj] : symmetric_sizes) {
        EXPECT_EQ(plain.ni_global(s), 7) << name(s);
        EXPECT_EQ(plain.nj_global(s), 5) << name(s);
        EXPECT_EQ(symmetric.ni_global(s), ni) << name(s);
        EXPECT_EQ(symmetric.nj_global(s), nj) << name(s);
    }
}

TEST(IoDecomp, WindowsFollowLocalTileOrder) {
    const IoDecomp symmetric(two_by_two(true, {{1, 1}, {0, 0}}));
    ASSERT_EQ(symmetric.n_local_tiles(), 2);
    EXPECT_EQ(symmetric.tile(0), (Tile{1, 1}));
    EXPECT_EQ(symmetric.tile(1), (Tile{0, 0}));
    EXPECT_EQ(symmetric.window(0), ij_box(3, 6, 2, 4));
    EXPECT_EQ(symmetric.window(1), ij_box(0, 2, 0, 1));
    EXPECT_EQ(symmetric.window(0, Stagger::XFace), ij_box(3, 7, 2, 4));
    EXPECT_EQ(symmetric.window(0, Stagger::YFace), ij_box(3, 6, 2, 5));
    EXPECT_EQ(symmetric.window(0, Stagger::Node), ij_box(3, 7, 2, 5));
    EXPECT_EQ(symmetric.window(1, Stagger::Node), ij_box(0, 3, 0, 2));

    // In a non-symmetric file the window of every staggering is the tile's cells.
    const IoDecomp plain(two_by_two(false, {{1, 1}}));
    for (const Stagger s : kStaggers) {
        EXPECT_EQ(plain.window(0, s), ij_box(3, 6, 2, 4)) << name(s);
    }

    // A compute rank may hold no tiles.
    EXPECT_EQ(IoDecomp(two_by_two(true, {})).n_local_tiles(), 0);
}

TEST(IoDecomp, EqualityComparesValues) {
    const IoDecomp a(two_by_two(true, {{0, 0}}, {{1, 0}, {0, 1}}));
    const IoDecomp same_elimination_other_order(two_by_two(true, {{0, 0}}, {{0, 1}, {1, 0}}));
    const IoDecomp other_local(two_by_two(true, {{1, 1}}, {{1, 0}, {0, 1}}));
    const IoDecomp other_symmetry(two_by_two(false, {{0, 0}}, {{1, 0}, {0, 1}}));
    EXPECT_EQ(a, same_elimination_other_order);
    EXPECT_NE(a, other_local);
    EXPECT_NE(a, other_symmetry);
}

// ---- write ownership ----

TEST(IoDecomp, UnmaskedSymmetricWrites) {
    // One compute rank holding all four tiles. The windows of tiles 0 and 1
    // overlap on the face line i = 3 between them ...
    const IoDecomp decomp(two_by_two(true, {{0, 0}, {1, 0}, {0, 1}, {1, 1}}));
    EXPECT_EQ(decomp.window(0, Stagger::XFace), ij_box(0, 3, 0, 1));
    EXPECT_EQ(decomp.window(1, Stagger::XFace), ij_box(3, 7, 0, 1));
    // ... but every face is written once. Each tile writes its cells, plus the
    // face lines on the global east and north edges.
    EXPECT_TRUE(writes_exactly(decomp, 0, Stagger::XFace, {ij_box(0, 2, 0, 1)}));
    EXPECT_TRUE(writes_exactly(decomp, 1, Stagger::XFace, {ij_box(3, 7, 0, 1)}));
    EXPECT_TRUE(writes_exactly(decomp, 0, Stagger::Node, {ij_box(0, 2, 0, 1)}));
    EXPECT_TRUE(writes_exactly(decomp, 1, Stagger::Node, {ij_box(3, 7, 0, 1)}));
    EXPECT_TRUE(writes_exactly(decomp, 2, Stagger::Node, {ij_box(0, 2, 2, 5)}));
    EXPECT_TRUE(writes_exactly(decomp, 3, Stagger::Node, {ij_box(3, 7, 2, 5)}));
    EXPECT_TRUE(writes_exactly(decomp, 0, Stagger::Cell, {ij_box(0, 2, 0, 1)}));
}

TEST(IoDecomp, NonSymmetricWritesAreTheCells) {
    // Even next to an eliminated tile: the face it shares belongs to its own cell.
    const IoDecomp decomp(two_by_two(false, {{0, 0}}, {{1, 0}, {0, 1}, {1, 1}}));
    for (const Stagger s : kStaggers) {
        EXPECT_TRUE(writes_exactly(decomp, 0, s, {ij_box(0, 2, 0, 1)}));
    }
}

TEST(IoDecomp, EliminatedNeighborsLeaveTheirEdgesToLiveTiles) {
    // A 3 x 2 layout of 2 x 2 tiles with tile (1, 1) eliminated.
    const IoDecompSpec spec{.tile_ni = {2, 2, 2},
                            .tile_nj = {2, 2},
                            .symmetric = true,
                            .eliminated_tiles = {{1, 1}},
                            .local_tiles = {{1, 0}, {0, 1}}};
    const IoDecomp decomp(spec);
    // Tile (1, 0) writes the north row its eliminated north neighbor would
    // have written, except the west end, which the live northwest tile (0, 1)
    // writes on its east face line. What tile (1, 0) writes is not a rectangle.
    EXPECT_TRUE(writes_exactly(decomp, 0, Stagger::Node, {ij_box(2, 3, 0, 1), ij_box(3, 3, 2, 2)}));
    EXPECT_TRUE(writes_exactly(decomp, 0, Stagger::YFace, {ij_box(2, 3, 0, 2)}));
    // Tile (0, 1) has an eliminated east neighbor and sits on the north edge,
    // so it writes its whole window.
    EXPECT_TRUE(writes_exactly(decomp, 1, Stagger::Node, {ij_box(0, 2, 2, 4)}));
    EXPECT_TRUE(writes_exactly(decomp, 1, Stagger::XFace, {ij_box(0, 2, 2, 3)}));
}

TEST(IoDecomp, WriteMaskListsTheWindowRowByRow) {
    // Tile (1, 0) of the 2 x 2 layout, with every tile live. Its Node window
    // is i 3..7 by j 0..2. It writes rows 0 and 1 in full (index 7 is on the
    // global east edge) but not row 2, which belongs to the north tile.
    const IoDecomp decomp(two_by_two(true, {{1, 0}}));
    const std::vector<bool> expected{true,  true,  true,  true,  true,
                                     true,  true,  true,  true,  true,
                                     false, false, false, false, false};
    EXPECT_EQ(decomp.write_mask(0, Stagger::Node), expected);
}

TEST(IoDecomp, WritesCoverEachPointOnceUnderRandomMasks) {
    std::mt19937 rng(20261114);
    for (int trial = 0; trial < 400; ++trial) {
        const Layout layout(rng);
        for (const bool symmetric : {false, true}) {
            std::vector<IoDecomp> decomps;
            for (const std::vector<Tile>& local : layout.compute_ranks) {
                decomps.emplace_back(IoDecompSpec{.tile_ni = layout.tile_ni,
                                                  .tile_nj = layout.tile_nj,
                                                  .symmetric = symmetric,
                                                  .eliminated_tiles = layout.eliminated,
                                                  .local_tiles = local});
            }
            for (const Stagger s : kStaggers) {
                const int ni = decomps.front().ni_global(s);
                const int nj = decomps.front().nj_global(s);
                std::vector<int> times_written(static_cast<std::size_t>(ni) * nj, 0);
                std::vector<Tile> writer(times_written.size(), Tile{-1, -1});
                for (const IoDecomp& decomp : decomps) {
                    for (int b = 0; b < decomp.n_local_tiles(); ++b) {
                        const IndexBox window = decomp.window(b, s);
                        const std::vector<bool> mask = decomp.write_mask(b, s);
                        ASSERT_EQ(mask.size(), static_cast<std::size_t>(window.n_points()));
                        std::size_t k = 0;
                        for (int j = window.lo[1]; j <= window.hi[1]; ++j) {
                            for (int i = window.lo[0]; i <= window.hi[0]; ++i, ++k) {
                                if (mask[k]) {
                                    const std::size_t at = static_cast<std::size_t>(j) * ni + i;
                                    ++times_written[at];
                                    writer[at] = decomp.tile(b);
                                }
                            }
                        }
                    }
                }
                for (int q = 0; q < nj; ++q) {
                    for (int p = 0; p < ni; ++p) {
                        const std::size_t at = static_cast<std::size_t>(q) * ni + p;
                        const std::optional<Tile> expected = layout.owner(symmetric, s, p, q);
                        const bool ok = times_written[at] == (expected ? 1 : 0) &&
                                        (!expected || writer[at] == *expected);
                        ASSERT_TRUE(ok) << "trial " << trial << " " << name(s)
                                        << (symmetric ? " symmetric" : "") << " point (" << p << ", " << q
                                        << ") written " << times_written[at] << " times";
                    }
                }
            }
        }
    }
}

// ---- invalid input ----

TEST(IoDecompDeathTest, InvalidSpecsAbort) {
    EXPECT_DEATH(IoDecomp({}), "tile_ni is empty");
    EXPECT_DEATH(IoDecomp({.tile_ni = {3}}), "tile_nj is empty");
    EXPECT_DEATH(IoDecomp({.tile_ni = {3, 0}, .tile_nj = {2}}), "tile_ni has a width of 0");
    EXPECT_DEATH(IoDecomp({.tile_ni = {3}, .tile_nj = {-1}}), "tile_nj has a width of -1");
    EXPECT_DEATH(IoDecomp({.tile_ni = {1 << 30, 1 << 30}, .tile_nj = {2}}),
                 "add up to more than an int");
    EXPECT_DEATH(IoDecomp(two_by_two(false, {}, {{2, 0}})), "eliminated tile .* is outside the 2 x 2");
    EXPECT_DEATH(IoDecomp(two_by_two(false, {}, {{1, 0}, {1, 0}})), "eliminated tile .* is listed twice");
    EXPECT_DEATH(IoDecomp(two_by_two(false, {{0, -1}})), "local tile .* is outside the 2 x 2");
    EXPECT_DEATH(IoDecomp(two_by_two(false, {{1, 0}}, {{1, 0}})), "local tile .* is also eliminated");
    EXPECT_DEATH(IoDecomp(two_by_two(false, {{0, 0}, {0, 0}})), "local tile .* is listed twice");
}

TEST(IoDecompDeathTest, OutOfRangeTileIndexAborts) {
    const IoDecomp decomp(two_by_two(false, {{0, 0}, {1, 1}}));
    EXPECT_DEATH(decomp.tile(2), "local tile index 2 is out of range");
    EXPECT_DEATH(decomp.window(-1), "local tile index -1 is out of range");
    EXPECT_DEATH(decomp.write_mask(2, Stagger::Node), "local tile index 2 is out of range");
}
