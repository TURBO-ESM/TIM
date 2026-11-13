#pragma once
/**
 * @file tim_io_decomp.hpp
 * @brief The I/O decomposition: which part of a horizontal field each
 * compute rank holds, in the index space of a file.
 */

#include <array>
#include <cstdint>
#include <vector>

#include "core/tim_stagger.hpp"

namespace TIM {

/// @brief A rectangle of 2-D indices with inclusive bounds: (i, j) is inside
/// when lo[0] <= i <= hi[0] and lo[1] <= j <= hi[1]. The box is empty when
/// hi is below lo in either direction. The default box is empty. Each
/// length, hi - lo + 1, must fit in an int.
///
/// An IndexBox carries no staggering: its indices are in the file index
/// space of one variable, and the IoDecomp notes say what an index names for
/// each staggering. Unlike amrex::Box, two boxes with the same bounds are
/// equal whatever their staggering.
struct IndexBox {
    std::array<int, 2> lo{0, 0};    ///< The lowest (i, j) in the box.
    std::array<int, 2> hi{-1, -1};  ///< The highest (i, j) in the box.

    /// @brief The number of points along i.
    /// @return hi[0] - lo[0] + 1, or 0 if that is negative.
    int length_i() const { return hi[0] < lo[0] ? 0 : hi[0] - lo[0] + 1; }

    /// @brief The number of points along j.
    /// @return hi[1] - lo[1] + 1, or 0 if that is negative.
    int length_j() const { return hi[1] < lo[1] ? 0 : hi[1] - lo[1] + 1; }

    /// @brief Whether the box holds no points.
    /// @return True when hi is below lo in either direction.
    bool empty() const { return hi[0] < lo[0] || hi[1] < lo[1]; }

    /// @brief The number of points in the box.
    /// @return length_i() * length_j(), as a 64-bit count.
    std::int64_t n_points() const { return std::int64_t{length_i()} * length_j(); }

    /// @brief Whether a point lies in the box.
    /// @param i The point's i index.
    /// @param j The point's j index.
    /// @return True when (i, j) is within the bounds.
    bool contains(const int i, const int j) const {
        return lo[0] <= i && i <= hi[0] && lo[1] <= j && j <= hi[1];
    }

    /// Member-wise equality (two empty boxes with different bounds differ).
    friend bool operator==(const IndexBox&, const IndexBox&) = default;
};

/// @brief The construction specification of an IoDecomp.
///
/// The global cell grid is cut into a layout of tiles: columns of cells, west
/// to east, by rows of cells, south to north. Tile (ti, tj) is in tile column
/// ti and tile row tj, both 0-based. Each tile that is not eliminated is held
/// by exactly one compute rank.
struct IoDecompSpec {
    /// @brief The width in cells of each tile column, west to east. Must not
    /// be empty, and every width must be positive. The widths add up to the
    /// global i-extent.
    std::vector<int> tile_ni;
    /// @brief The height in cells of each tile row, south to north. Must not
    /// be empty, and every height must be positive. The heights add up to the
    /// global j-extent.
    std::vector<int> tile_nj;
    /// @brief Whether files also hold the faces on the global west and south
    /// edges. See the IoDecomp notes.
    bool symmetric = false;
    /// @brief The tiles (ti, tj) that no compute rank holds, such as the
    /// land-only tiles dropped by a mask table. Each one listed once.
    std::vector<std::array<int, 2>> eliminated_tiles;
    /// @brief The tiles (ti, tj) this compute rank holds, each listed once and
    /// none eliminated. Entry b becomes local tile b. May be empty.
    std::vector<std::array<int, 2>> local_tiles;
};

/// @brief Where each compute rank's part of a horizontal field sits in the
/// index space of a file, for every staggering. (Analogue of what FMS I/O takes
/// from a domain2D: the global extents, the layout, each PE's compute domain
/// and the mask table.)
///
/// Terms used here:
/// - **cell**: (i, j) of the global grid of ni x nj cells.
/// - **tile**: one rectangle of cells in the layout, at (ti, tj), live or
///   eliminated. A **local tile** is one this compute rank holds, numbered b
///   in the order of IoDecompSpec::local_tiles.
/// - **layout**: the grid of tiles, layout_nx() columns by layout_ny() rows.
/// - **file index**: the global position (i, j) of a value in the file,
///   counting from 0. For cells, and for every staggering in a symmetric
///   file, it equals Domain's global index. In a non-symmetric file, a face
///   or corner has a file index one less than its Domain index along each
///   staggered direction, because the file leaves out the faces on the
///   global west and south edges.
/// - **window**: the rectangle of file indices that a local tile holds
///   values at, for one staggering. A compute rank has one window for each
///   local tile and staggering. window(b, s) returns it.
///
/// Windows are IndexBox values.
///
/// An IoDecomp describes compute ranks, the ranks that hold a field's
/// values. I/O ranks, which PIO picks to read and write the file, play no
/// part in it. An IoDecomp makes no MPI or PIO calls, so its methods never
/// wait on other compute ranks. Each compute rank builds its own value, with
/// the same layout and eliminated tiles and its own local tiles. The local
/// tiles of all compute ranks together must cover each live tile exactly
/// once. One compute rank cannot check that, so it is a precondition.
///
/// One IoDecomp serves every variable in a file: the staggering is an
/// argument of each query.
///
/// **File indices.** File index (i, j) of a Cell variable is cell (i, j). The
/// axis lengths of a variable in the file are:
///
/// | Stagger | i-length                      | j-length                      |
/// |---------|-------------------------------|-------------------------------|
/// | Cell    | ni                            | nj                            |
/// | XFace   | ni + 1 if symmetric, else ni  | nj                            |
/// | YFace   | ni                            | nj + 1 if symmetric, else nj  |
/// | Node    | ni + 1 if symmetric, else ni  | nj + 1 if symmetric, else nj  |
///
/// Along a staggered direction, what a file index names depends on
/// symmetric():
/// - In a symmetric file, index p is the face on the west (or south) side of
///   cell p. Index 0 is the global west (south) edge, and index ni (nj) the
///   east (north) edge. This is Stagger's convention and AMReX's nodal index.
/// - In a non-symmetric file, index p is the face on the east (or north) side
///   of cell p, which is AMReX nodal index p + 1. The global west (south)
///   edge is not in the file. This is MOM6's non-symmetric memory, where u(i)
///   sits east of h(i).
///
/// **Windows.** In a symmetric file a window reaches one index further at the
/// high end along each staggered direction, so it holds the faces on both
/// sides of the tile's cells, and two neighboring windows overlap on the face
/// line between them. In a non-symmetric file a window is the indices of the
/// tile's cells. For example, along i with ni = 7 and two local tiles holding
/// cells 0..2 and 3..6, the windows of an XFace variable are:
/// @verbatim
///   cells:            |  0  |  1  |  2  |  3  |  4  |  5  |  6  |
///   symmetric:        0     1     2     3     4     5     6     7
///     tile 0 window:  0 ............... 3
///     tile 1 window:                    3 ..................... 7
///   non-symmetric:          0     1     2     3     4     5     6
///     tile 0 window:        0 ......... 2
///     tile 1 window:                          3 ............... 6
/// @endverbatim
/// Reading and writing use windows differently. In a symmetric file, a face
/// or corner can lie in up to 4 windows. A reader lists each tile's whole
/// window: PIO's read-only decompositions accept an index that several
/// compute ranks list. A writer must list each index only once, so it keeps
/// only the indices of a window that write_mask() marks.
///
/// **Write ownership.** Along a staggered direction of a symmetric file,
/// index p touches cells p - 1 and p. Along any other direction it touches
/// only cell p. A cell is live when its tile is not eliminated. Each point is
/// written by the tile of the first live cell touching it, checked in the
/// order (i, j), (i - 1, j), (i, j - 1), (i - 1, j - 1). As a result:
/// - each point that touches a live cell is written exactly once. This holds
///   for the faces next to an eliminated tile, which a live neighbor writes,
///   and for a face between two local tiles of one compute rank;
/// - a point that touches only eliminated cells is written by nobody and
///   keeps the file's fill value;
/// - without eliminated tiles, each tile writes its own cells' indices, and
///   the tiles on the global east and north edges also write the face line
///   there.
///
/// Equality compares this compute rank's value only. Two compute ranks can
/// disagree on whether two decompositions are equal (when these differ only
/// in local tiles), so a choice that leads to a collective call must not rest
/// on operator== alone.
class IoDecomp {
public:
    /// @brief Build the decomposition. Aborts on an invalid spec: an empty
    /// tile list or a width that is not positive, a global extent too large
    /// for an int, a listed tile outside the layout, a tile listed twice, or
    /// a local tile that is also eliminated.
    /// @param spec The specification; see IoDecompSpec.
    explicit IoDecomp(IoDecompSpec spec);

    /// @brief The global i-length of a variable in the file.
    /// @param s The variable's staggering.
    /// @return ni, plus 1 for XFace and Node in a symmetric file.
    int ni_global(Stagger s = Stagger::Cell) const {
        return tile_lo_i_.back() + (grows_i(s) ? 1 : 0);
    }

    /// @brief The global j-length of a variable in the file.
    /// @param s The variable's staggering.
    /// @return nj, plus 1 for YFace and Node in a symmetric file.
    int nj_global(Stagger s = Stagger::Cell) const {
        return tile_lo_j_.back() + (grows_j(s) ? 1 : 0);
    }

    /// @brief Whether files hold the faces on the global west and south edges.
    /// @return The symmetric flag of the spec.
    bool symmetric() const { return symmetric_; }

    /// @brief The number of tile columns.
    /// @return The size of the spec's tile_ni.
    int layout_nx() const { return static_cast<int>(tile_lo_i_.size()) - 1; }

    /// @brief The number of tile rows.
    /// @return The size of the spec's tile_nj.
    int layout_ny() const { return static_cast<int>(tile_lo_j_.size()) - 1; }

    /// @brief The number of tiles this compute rank holds.
    /// @return The size of the spec's local_tiles.
    int n_local_tiles() const { return static_cast<int>(local_tiles_.size()); }

    /// @brief The layout position of a local tile.
    /// @param b The local tile index. Aborts if out of range.
    /// @return {ti, tj}.
    std::array<int, 2> tile(int b) const;

    /// @brief A local tile's window for a staggering; see the notes above.
    /// @param b The local tile index. Aborts if out of range.
    /// @param s The staggering.
    /// @return The tile's cells, grown by one at the high end along each
    ///         staggered direction in a symmetric file.
    IndexBox window(int b, Stagger s = Stagger::Cell) const;

    /// @brief Which indices of a local tile's window this tile writes, under
    /// the write-ownership rule above.
    ///
    /// The mask lists the window's indices row by row, with i running
    /// fastest, and is true where this tile writes the index. A tile always
    /// writes its own cells' indices. A writer builds its PIO map from the
    /// window and puts 0 wherever the mask is false.
    /// @param b The local tile index. Aborts if out of range.
    /// @param s The staggering.
    /// @return One entry per index of window(b, s).
    std::vector<bool> write_mask(int b, Stagger s) const;

    /// Member-wise equality of this compute rank's value; see the class notes.
    friend bool operator==(const IoDecomp&, const IoDecomp&) = default;

private:
    /// @brief Whether a tile is held by some compute rank: false for an
    /// eliminated tile or one outside the layout.
    bool tile_live(int ti, int tj) const;

    /// @brief Whether a staggering adds a point along i in this file.
    bool grows_i(Stagger s) const {
        return symmetric_ && (s == Stagger::XFace || s == Stagger::Node);
    }
    /// @brief Whether a staggering adds a point along j in this file.
    bool grows_j(Stagger s) const {
        return symmetric_ && (s == Stagger::YFace || s == Stagger::Node);
    }

    /// @brief The first cell of each tile column, then the global i-extent.
    std::vector<int> tile_lo_i_;
    /// @brief The first cell of each tile row, then the global j-extent.
    std::vector<int> tile_lo_j_;
    bool symmetric_ = false;
    /// @brief Whether each tile is live, indexed tj * layout_nx() + ti.
    std::vector<bool> live_;
    /// @brief The tiles this compute rank holds, in local tile order.
    std::vector<std::array<int, 2>> local_tiles_;
};

}  // namespace TIM
