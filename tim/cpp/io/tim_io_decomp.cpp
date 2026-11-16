/**
 * @file tim_io_decomp.cpp
 * @brief Implementation of TIM::IoDecomp: spec checks, windows and the
 * write-ownership test.
 */

#include "io/tim_io_decomp.hpp"

#include <climits>
#include <cstddef>
#include <string>
#include <utility>

#include "core/tim_abort.hpp"

namespace TIM {

namespace {

std::string tile_name(const std::array<int, 2>& t) {
    return "(" + std::to_string(t[0]) + ", " + std::to_string(t[1]) + ")";
}

// The first cell of each tile column (or row), then the global extent.
// Aborts on an empty list, a width that is not positive, or an extent that
// leaves no room for the extra point of a staggered axis.
std::vector<int> tile_origins(const std::vector<int>& widths, const char* name) {
    if (widths.empty()) {
        TIM::abort(std::string("TIM::IoDecomp: ") + name + " is empty.");
    }
    std::vector<int> lo{0};
    std::int64_t extent = 0;
    for (const int w : widths) {
        if (w <= 0) {
            TIM::abort(std::string("TIM::IoDecomp: ") + name + " has a width of " +
                       std::to_string(w) + "; widths must be positive.");
        }
        extent += w;
        if (extent > INT_MAX - 1) {
            TIM::abort(std::string("TIM::IoDecomp: the widths in ") + name +
                       " add up to more than an int can index.");
        }
        lo.push_back(static_cast<int>(extent));
    }
    return lo;
}

}  // namespace

IoDecomp::IoDecomp(IoDecompSpec spec)
    : tile_lo_i_(tile_origins(spec.tile_ni, "tile_ni")),
      tile_lo_j_(tile_origins(spec.tile_nj, "tile_nj")),
      symmetric_(spec.symmetric),
      local_tiles_(std::move(spec.local_tiles)) {
    const int nx = layout_nx();
    const int ny = layout_ny();
    const auto in_layout = [nx, ny](const std::array<int, 2>& t) {
        return 0 <= t[0] && t[0] < nx && 0 <= t[1] && t[1] < ny;
    };
    const auto index = [nx](const std::array<int, 2>& t) {
        return static_cast<std::size_t>(t[1]) * nx + t[0];
    };
    const std::string layout = std::to_string(nx) + " x " + std::to_string(ny) + " layout";

    live_.assign(static_cast<std::size_t>(nx) * ny, true);
    for (const auto& t : spec.eliminated_tiles) {
        if (!in_layout(t)) {
            TIM::abort("TIM::IoDecomp: eliminated tile " + tile_name(t) + " is outside the " +
                       layout + ".");
        }
        if (!live_[index(t)]) {
            TIM::abort("TIM::IoDecomp: eliminated tile " + tile_name(t) + " is listed twice.");
        }
        live_[index(t)] = false;
    }

    std::vector<bool> held(live_.size(), false);
    for (const auto& t : local_tiles_) {
        if (!in_layout(t)) {
            TIM::abort("TIM::IoDecomp: local tile " + tile_name(t) + " is outside the " + layout +
                       ".");
        }
        if (!live_[index(t)]) {
            TIM::abort("TIM::IoDecomp: local tile " + tile_name(t) + " is also eliminated.");
        }
        if (held[index(t)]) {
            TIM::abort("TIM::IoDecomp: local tile " + tile_name(t) + " is listed twice.");
        }
        held[index(t)] = true;
    }
}

bool IoDecomp::tile_live(const int ti, const int tj) const {
    if (ti < 0 || ti >= layout_nx() || tj < 0 || tj >= layout_ny()) {
        return false;
    }
    return live_[static_cast<std::size_t>(tj) * layout_nx() + ti];
}

std::array<int, 2> IoDecomp::tile(const int b) const {
    if (b < 0 || b >= n_local_tiles()) {
        TIM::abort("TIM::IoDecomp: local tile index " + std::to_string(b) +
                   " is out of range; this compute rank holds " + std::to_string(n_local_tiles()) +
                   " tiles.");
    }
    return local_tiles_[static_cast<std::size_t>(b)];
}

IndexBox IoDecomp::window(const int b, const Stagger s) const {
    const auto [ti, tj] = tile(b);
    IndexBox result{{tile_lo_i_[ti], tile_lo_j_[tj]},
                    {tile_lo_i_[ti + 1] - 1, tile_lo_j_[tj + 1] - 1}};
    if (grows_i(s)) {
        result.hi[0] += 1;
    }
    if (grows_j(s)) {
        result.hi[1] += 1;
    }
    return result;
}

// The ownership rule in the header, worked out once per tile. Tiles form a
// layout, so every index on a face line along a tile's side has the same
// neighbors, and one liveness test decides the whole line.
std::vector<bool> IoDecomp::write_mask(const int b, const Stagger s) const {
    const IndexBox w = window(b, s);
    const IndexBox c = window(b, Stagger::Cell);
    const auto [ti, tj] = tile(b);

    // The tile beyond a face line comes before this tile in the order, so
    // this tile writes the line only if that tile is eliminated.
    const bool east_line = !tile_live(ti + 1, tj);
    const bool north_line = !tile_live(ti, tj + 1);
    // The northeast corner comes after the northeast, north and east tiles.
    const bool corner = east_line && north_line && !tile_live(ti + 1, tj + 1);
    // For Node, the north line's west end also touches the northwest tile,
    // which comes before this tile too.
    const bool northwest_takes_west_end = grows_i(s) && tile_live(ti - 1, tj + 1);

    std::vector<bool> mask(static_cast<std::size_t>(w.n_points()));
    std::size_t k = 0;
    for (int j = w.lo[1]; j <= w.hi[1]; ++j) {
        const bool on_north_line = j > c.hi[1];
        for (int i = w.lo[0]; i <= w.hi[0]; ++i, ++k) {
            const bool on_east_line = i > c.hi[0];
            if (on_east_line && on_north_line) {
                mask[k] = corner;
            } else if (on_east_line) {
                mask[k] = east_line;
            } else if (on_north_line) {
                mask[k] = north_line && !(northwest_takes_west_end && i == c.lo[0]);
            } else {
                mask[k] = true;  // the tile's own cells
            }
        }
    }
    return mask;
}

}  // namespace TIM
