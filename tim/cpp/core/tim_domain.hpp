#pragma once
/**
 * @file tim_domain.hpp
 * @brief The horizontal computational domain: global index space, connectivity,
 * and the decomposition into rank-assigned boxes.
 */

#include <array>
#include <optional>
#include <vector>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Periodicity.H>

#include "tim_stagger.hpp"

namespace TIM {

/// @brief The construction specification of a Domain.
struct DomainSpec {
    int ni_global = 0;   ///< Global number of cells in the i-direction (x).
    int nj_global = 0;   ///< Global number of cells in the j-direction (y).
    int ni_halo = 0;     ///< Halo width in the i-direction (metadata; see Domain docs).
    int nj_halo = 0;     ///< Halo width in the j-direction (metadata; see Domain docs).
    bool periodic_x = false;  ///< Whether the i-direction is periodic (cyclic).
    bool periodic_y = false;  ///< Whether the j-direction is periodic (cyclic).
    bool tripolar_n = false;  ///< Whether the domain has tripolar connectivity (a
                              ///< fold) at the northern edge. Not implemented yet.
    /// @brief Number of boxes to decompose the domain into; nullopt (the
    /// default) means one box per rank. The actual box count can be smaller
    /// if the domain is too small to split that finely.
    std::optional<int> n_boxes = std::nullopt;
};

/// @brief The optional arguments of Domain::make_field.
struct FieldOpts {
    int ncomp = 1;  ///< Number of field components. Must be positive.
    /// @brief Ghost-cell width of the field, applied in both horizontal
    /// directions. The default (nullopt) is the domain's halo widths.
    std::optional<int> nghost = std::nullopt;
};

/// @brief The horizontal computational domain and its decomposition
/// (Analogue of FMS mpp_domains, rebuilt on AMReX). A Domain owns the global
/// cell-centered index space, the connectivity flags, and the decomposition of
/// that index space into boxes assigned to MPI ranks (an amrex::BoxArray and an
/// amrex::DistributionMapping). Clients create distributed fields from
/// these products.
///
/// Design notes:
/// - This Domain class is model-agnostic.
/// - The decomposition is horizontal-only: boxes span the full k-range and
///   are never split in the vertical. By default the global domain is split
///   into one box per rank. A different box count can be requested for testing
///   or finer-grained load balancing.
/// - Halo widths are carried as metadata only. In AMReX, halos are a
///   per-field property (the ghost cells of a MultiFab), so the halo widths
///   stored here are consumed at field-creation sites rather than baked into
///   the domain's index space.
/// - A staggered field is not a clean tiling of its index space. Converting
///   the cell-centered decomposition to a face or node index type grows every
///   box by one plane per nodal direction, so the boxes on either side of an
///   internal boundary each hold a copy of the plane between them, and on a
///   periodic axis the plane at the far edge is the image of the one at the
///   near edge. This mirrors FMS symmetric memory, where neighboring PEs
///   likewise both hold a shared u/v/q edge. The contract for these shared
///   planes:
///   - Ownership follows AMReX's convention: the lowest-index box containing
///     a point owns it. That is the rule FabArray::OwnerMask, sum_unique,
///     and OverrideSync all apply, so a reduction that must count each point
///     once uses sum_unique (or masks with OwnerMask); a plain reduction
///     counts the stored copies.
///   - The copies are kept in agreement by computing shared planes
///     redundantly, from a rule both boxes evaluate identically (e.g.
///     analytically from global indices). Where redundant computation is not
///     possible, OverrideSync(periodicity()) reconciles them afterwards by
///     the same ownership rule.
///   FillBoundary fills ghost cells from valid ones and leaves disagreeing
///   valid copies alone, so agreement is the producer's obligation; halo
///   exchanges do not restore it.
/// - todo: this class is the intended producer of TIM::IoDecomp values
///   ("Decomp2D produced from AMReX distribution maps"); an io_decomp()
///   method will be added here once TIM's parallel IO layer merges. Its
///   initial precondition is the default one-box-per-rank decomposition
///   (IoDecomp describes one window per rank); supporting many boxes per
///   rank (GPU tiling, makeSFC/makeKnapSack load balancing) requires
///   generalizing the I/O DofMap generation to a list of windows per rank.
/// - todo: land-block elimination (the FMS mask_table analogue, e.g.
///   cesm_t232's masked 25x40 layout) is a planned Domain feature: drop
///   all-land boxes from the BoxArray before distribution. The I/O layer's
///   write-edge ownership is already mask-aware. This class would supply the
///   tile-liveness mask. At that point the tile layout is recorded at
///   decomposition time rather than derived from surviving box corners (the
///   constructor's derived-origin consistency checks assume an unmasked
///   grid).
class Domain {
public:
    /// @brief Build the domain and its horizontal decomposition.
    /// @param spec The domain specification; see DomainSpec. Aborts on an
    ///        invalid or inconsistent specification.
    /// @pre The runtime is up (a TIM::Runtime exists); aborts otherwise.
    explicit Domain(const DomainSpec& spec);

    /// @brief Global number of cells in the i-direction (x).
    /// @return The global i-extent.
    int ni_global() const { return geometry_2d_.Domain().length(0); }

    /// @brief Global number of cells in the j-direction (y).
    /// @return The global j-extent.
    int nj_global() const { return geometry_2d_.Domain().length(1); }

    /// @brief Halo width in the i-direction (metadata for field creation).
    /// @return The i-direction halo width.
    int ni_halo() const { return ni_halo_; }

    /// @brief Halo width in the j-direction (metadata for field creation).
    /// @return The j-direction halo width.
    int nj_halo() const { return nj_halo_; }

    /// @brief The halo widths as a MultiFab ghost-cell vector. Horizontal-only:
    /// the k-component is always 0.
    /// @return {ni_halo, nj_halo, 0}.
    amrex::IntVect nghost() const { return amrex::IntVect(ni_halo_, nj_halo_, 0); }

    /// @brief True if the i-direction is periodic (cyclic).
    /// @return The i-direction periodicity flag.
    bool periodic_x() const { return geometry_2d_.isPeriodic(0); }

    /// @brief True if the j-direction is periodic (cyclic).
    /// @return The j-direction periodicity flag.
    bool periodic_y() const { return geometry_2d_.isPeriodic(1); }

    /// @brief True if the domain has tripolar connectivity at the northern edge.
    /// @return The tripolar connectivity flag.
    bool tripolar_n() const { return tripolar_n_; }

    /// @brief The cell-centered decomposition of the global domain, extended
    /// to the requested vertical extent. Every box spans the full
    /// k-range [0, nk): the decomposition is horizontal-only.
    /// @param nk Number of vertical points of the field to be created:
    ///        1 for 2-D fields, NK for 3-D layer fields, NK+1 for interface
    ///        fields, etc. Aborts if not positive.
    /// @return The cell-centered BoxArray with the requested k-extent.
    amrex::BoxArray boxArray(int nk) const;

    /// @brief The assignment of the horizontal boxes to MPI ranks.
    /// @return The distribution mapping of the horizontal decomposition.
    const amrex::DistributionMapping& distribution_mapping() const {
        return distribution_mapping_;
    }

    /// @brief Create a distributed field on this domain's decomposition: a
    /// MultiFab with the requested staggering, vertical extent, and component
    /// count, carrying the domain's halo widths as ghost cells unless overridden.
    /// @param stagger Where the field's values sit within a grid cell.
    /// @param nk Number of vertical points of the field. Must be positive.
    /// @param opts The optional arguments (component count, ghost-cell
    ///        width); see FieldOpts.
    /// @return The newly created field, with uninitialized contents.
    /// @note A non-Cell staggering duplicates the planes that neighboring
    ///       boxes share; see the Domain class notes before reducing over
    ///       such a field.
    amrex::MultiFab make_field(Stagger stagger, int nk,
                               FieldOpts opts = {}) const;

    /// @brief The domain's periodicity in index space, for halo exchanges
    /// (e.g. MultiFab::FillBoundary).
    /// @return The global extent in each periodic direction, 0 in each
    ///         non-periodic one.
    amrex::Periodicity periodicity() const;

    /// @brief Number of boxes in the horizontal decomposition.
    /// @return The box count.
    int n_boxes() const { return static_cast<int>(box_array_2d_.size()); }

    /// @brief Number of tiles in the i-direction of the logical tile grid.
    /// @return Tiles in i.
    int layout_nx() const { return static_cast<int>(tile_lo_i_.size()); }

    /// @brief Number of tiles in the j-direction of the logical tile grid.
    /// @return Tiles in j.
    int layout_ny() const { return static_cast<int>(tile_lo_j_.size()); }

    /// @brief The tile coordinates of a box in the logical
    /// layout_nx() x layout_ny() tile grid.
    /// @param box_index Index of the box in boxArray() /
    ///        distribution_mapping(). Aborts if out of range.
    /// @return {tile_i, tile_j}, 0-based.
    std::array<int, 2> tileOf(int box_index) const;

private:
    int ni_halo_;
    int nj_halo_;
    bool tripolar_n_;

    /// @brief The AMReX geometry of the 2D (horizontal) global index space:
    /// owns the connectivity/periodicity detail. Index-space only: its
    /// RealBox stays a placeholder unit box, since physical metrics live
    /// with the horizontal grid.
    amrex::Geometry geometry_2d_;
    /// @brief The 2D (horizontal) decomposition.
    amrex::BoxArray box_array_2d_;
    /// @brief The assignment of the horizontal boxes to MPI ranks.
    amrex::DistributionMapping distribution_mapping_;
    /// @brief Sorted distinct box-corner coordinates: the tile-column and
    /// tile-row origins of the logical tile grid (sizes = layout shape).
    std::vector<int> tile_lo_i_;
    std::vector<int> tile_lo_j_;
};

}  // namespace TIM
