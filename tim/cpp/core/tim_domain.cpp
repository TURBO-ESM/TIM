/**
 * @file tim_domain.cpp
 * @brief Implementation of TIM::Domain, the horizontal computational domain
 * and its decomposition.
 */

#include <algorithm>
#include <set>
#include <string>
#include <utility>

#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_IntVect.H>
#include <AMReX_ParallelDescriptor.H>

#include "tim_domain.hpp"
#include "tim_runtime.hpp"

namespace TIM {

Domain::Domain(const DomainSpec& spec)
    : ni_halo_(spec.ni_halo), nj_halo_(spec.nj_halo),
      tripolar_n_(spec.tripolar_n) {

    if (!Runtime::active()) {
        TIM::abort("TIM::Domain: the infrastructure runtime is not up; "
                   "construct a TIM::Runtime before constructing a Domain.");
    }
    if (spec.ni_global <= 0 || spec.nj_global <= 0) {
        TIM::abort("TIM::Domain: global extents must be positive.");
    }
    if (spec.ni_halo < 0 || spec.nj_halo < 0) {
        TIM::abort("TIM::Domain: halo widths must be non-negative.");
    }
    if (spec.periodic_y && spec.tripolar_n) {
        TIM::abort("TIM::Domain: periodic_y and tripolar_n cannot both be true.");
    }
    if (spec.tripolar_n) {
        TIM::abort("TIM::Domain: tripolar connectivity (tripolar_n) is not "
                   "implemented yet.");
    }
    if (spec.n_boxes.has_value() && *spec.n_boxes <= 0) {
        TIM::abort("TIM::Domain: n_boxes must be positive when specified.");
    }

    // The global cell-centered index space. The decomposition is
    // horizontal-only (k=0).
    const amrex::Box domain_2d(amrex::IntVect(0, 0, 0),
                               amrex::IntVect(spec.ni_global - 1, spec.nj_global - 1, 0));

    // The AMReX geometry of the index space, which owns the connectivity/
    // periodicity detail. Cartesian, always: this Geometry is index-space
    // bookkeeping. Physical metrics arrive as discrete fields with the
    // model-specific horizontal grid, as in MOM6, where the grid
    // (not the domain) carries dx/dy/area arrays.
    geometry_2d_.define(domain_2d,
                        amrex::RealBox({AMREX_D_DECL(0., 0., 0.)},
                                       {AMREX_D_DECL(1., 1., 1.)}),
                        amrex::CoordSys::cartesian,
                        {spec.periodic_x ? 1 : 0, spec.periodic_y ? 1 : 0, 0});

    // Split into near-square horizontal boxes, one per rank by default (the
    // requested count is an upper bound: a domain too small to split that
    // finely yields fewer boxes).
    const int nranks = amrex::ParallelDescriptor::NProcs();
    const int target_boxes = spec.n_boxes.value_or(nranks);
    box_array_2d_ = amrex::decompose(domain_2d, target_boxes,
                                     {true, true, false});
    const int actual_boxes = n_boxes();

    // The default one-box-per-rank decomposition can fall short of the rank
    // count (e.g., for a small domain on many ranks), leaving some ranks with
    // no boxes. This is valid but likely a misconfiguration, so warn the user
    // about it. (An explicit n_boxes request is an upper bound by contract,
    // so no warning in that case.)
    if (!spec.n_boxes.has_value() && actual_boxes < nranks) {
        amrex::Warning("TIM::Domain: the default one-box-per-rank decomposition "
                       "yielded " + std::to_string(actual_boxes) + " boxes for " +
                       std::to_string(nranks) + " ranks; the ranks beyond the box "
                       "count own no cells.");
    }

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
    if (static_cast<int>(tile_lo_i_.size() * tile_lo_j_.size()) != actual_boxes) {
        TIM::abort("TIM::Domain: the decomposition does not form a complete "
                   "rectangular tile grid (layout_nx * layout_ny != n_boxes).");
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

amrex::BoxArray Domain::boxArray(const int nk) const {
    if (nk <= 0) {
        TIM::abort("TIM::Domain::boxArray: nk must be positive.");
    }
    amrex::BoxArray box_array = box_array_2d_;
    box_array.growHi(2, nk - 1);  // k-range [0,0] -> [0,nk-1]
    return box_array;
}

amrex::MultiFab Domain::make_field(const Stagger stagger, const int nk,
                                   const FieldOpts opts) const {
    if (nk <= 0) {
        TIM::abort("TIM::Domain::make_field: nk must be positive.");
    }
    if (opts.ncomp <= 0) {
        TIM::abort("TIM::Domain::make_field: ncomp must be positive.");
    }
    if (opts.nghost && *opts.nghost < 0) {
        TIM::abort("TIM::Domain::make_field: the ghost-cell width must be non-negative.");
    }
    const amrex::IntVect ng =
        opts.nghost ? amrex::IntVect(*opts.nghost, *opts.nghost, 0) : this->nghost();
    return amrex::MultiFab(amrex::convert(boxArray(nk), nodality(stagger)),
                           distribution_mapping_, opts.ncomp, ng);
}

amrex::Periodicity Domain::periodicity() const {
    return geometry_2d_.periodicity();
}

}  // namespace TIM
