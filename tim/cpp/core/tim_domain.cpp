/**
 * @file tim_domain.cpp
 * @brief Implementation of TIM::Domain, the horizontal computational domain
 * and its decomposition.
 */

#include <algorithm>
#include <set>
#include <utility>

#include <AMReX_Box.H>
#include <AMReX_IntVect.H>
#include <AMReX_ParallelDescriptor.H>

#include "tim_domain.hpp"
#include "tim_runtime.hpp"

namespace TIM {

Domain::Domain(const int ni_global, const int nj_global,
               const int ni_halo, const int nj_halo,
               const bool periodic_x, const bool periodic_y,
               const bool tripolar_n, const int n_boxes_in)
    : ni_global_(ni_global), nj_global_(nj_global),
      ni_halo_(ni_halo), nj_halo_(nj_halo),
      periodic_x_(periodic_x), periodic_y_(periodic_y),
      tripolar_n_(tripolar_n) {

    if (!Runtime::active()) {
        TIM::abort("TIM::Domain: the infrastructure runtime is not up; "
                   "construct a TIM::Runtime before constructing a Domain.");
    }
    if (ni_global <= 0 || nj_global <= 0) {
        TIM::abort("TIM::Domain: global extents must be positive.");
    }
    if (ni_halo < 0 || nj_halo < 0) {
        TIM::abort("TIM::Domain: halo widths must be non-negative.");
    }
    if (periodic_y && tripolar_n) {
        TIM::abort("TIM::Domain: periodic_y and tripolar_n cannot both be true.");
    }
    if (tripolar_n) {
        TIM::abort("TIM::Domain: tripolar connectivity (tripolar_n) is not "
                   "implemented yet.");
    }
    if (n_boxes_in < 0) {
        TIM::abort("TIM::Domain: n_boxes_in must be positive, or 0 for one box per rank.");
    }

    // The global cell-centered index space. The decomposition is
    // horizontal-only, hence the single-level (k = 0) box.
    const amrex::Box domain_2d(amrex::IntVect(0, 0, 0),
                               amrex::IntVect(ni_global - 1, nj_global - 1, 0));

    // Split into near-square horizontal boxes, one per rank by default (the
    // requested count is an upper bound: a domain too small to split that
    // finely yields fewer boxes).
    const int nranks = amrex::ParallelDescriptor::NProcs();
    const int target_boxes = (n_boxes_in == 0) ? nranks : n_boxes_in;
    box_array_2d_ = amrex::decompose(domain_2d, target_boxes,
                                     {true, true, false});
    const int actual_boxes = n_boxes();

    // Deterministic round-robin assignment of boxes to ranks (sort=false:
    // sorting orders ranks by current load via a collective, making the map
    // nondeterministic). With the default one-box-per-rank decomposition
    // this is the identity map.
    // todo: explore other decomposition strategies available in AMReX (e.g.,
    //       DistributionMapping::makeSFC, makeKnapSack, etc.) particularly
    //       for many-boxes-per-rank scenarios (GPU tiling, load balancing, etc.)
    distribution_mapping_.RoundRobinProcessorMap(
        actual_boxes, nranks, /*sort=*/false);

    // Record the logical tile grid: the sorted distinct box corners are the
    // tile-column/row origins. The I/O decomposition model (TIM::IoDecomp) is
    // tile-based and needs the layout shape and tile coordinates. The
    // invariant it needs is that every box maps to a distinct (tile_i, tile_j).
    std::set<int> lo_i, lo_j;
    std::set<std::pair<int, int>> lo_pairs;
    for (int b = 0; b < actual_boxes; ++b) {
        const amrex::Box box = box_array_2d_[b];
        lo_i.insert(box.smallEnd(0));
        lo_j.insert(box.smallEnd(1));
        lo_pairs.emplace(box.smallEnd(0), box.smallEnd(1));
    }
    tile_lo_i_.assign(lo_i.begin(), lo_i.end());
    tile_lo_j_.assign(lo_j.begin(), lo_j.end());
    // A box's tile coordinates are determined by its lower corner, so
    // corner-pair uniqueness IS tile uniqueness.
    if (static_cast<int>(lo_pairs.size()) != actual_boxes) {
        TIM::abort("TIM::Domain: two boxes map to the same tile of the "
                   "logical tile grid; the I/O decomposition model requires "
                   "one box per tile.");
    }
}

std::array<int, 2> Domain::tileOf(const int box_index) const {
    if (box_index < 0 || box_index >= n_boxes()) {
        TIM::abort("TIM::Domain::tileOf: box index out of range.");
    }
    const amrex::Box box = box_array_2d_[box_index];
    const auto tile_i = std::lower_bound(tile_lo_i_.begin(), tile_lo_i_.end(),
                                         box.smallEnd(0)) - tile_lo_i_.begin();
    const auto tile_j = std::lower_bound(tile_lo_j_.begin(), tile_lo_j_.end(),
                                         box.smallEnd(1)) - tile_lo_j_.begin();
    return {static_cast<int>(tile_i), static_cast<int>(tile_j)};
}

amrex::BoxArray Domain::boxArray(const int n_levels) const {
    if (n_levels <= 0) {
        TIM::abort("TIM::Domain::boxArray: n_levels must be positive.");
    }
    amrex::BoxArray box_array = box_array_2d_;
    box_array.growHi(2, n_levels - 1);  // k-range [0,0] -> [0,n_levels-1]
    return box_array;
}

amrex::Periodicity Domain::periodicity() const {
    return amrex::Periodicity(amrex::IntVect(periodic_x_ ? ni_global_ : 0,
                                             periodic_y_ ? nj_global_ : 0,
                                             0));
}

}  // namespace TIM
