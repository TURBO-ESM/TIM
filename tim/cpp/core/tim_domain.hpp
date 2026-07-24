#pragma once
/**
 * @file tim_domain.hpp
 * @brief The horizontal computational domain: global index space, connectivity,
 * and the decomposition into rank-assigned boxes.
 */

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Periodicity.H>

namespace TIM {

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
/// - todo: this class is the intended producer of TIM::IoDecomp values
///   ("Decomp2D produced from AMReX distribution maps"); an io_decomp()
///   method will be added here once TIM's parallel IO layer merges.
class Domain {
public:
    /// @brief Build the domain and its horizontal decomposition.
    /// @param ni_global Global number of cells in the i-direction (x).
    /// @param nj_global Global number of cells in the j-direction (y).
    /// @param ni_halo   Halo width in the i-direction (metadata; see class docs).
    /// @param nj_halo   Halo width in the j-direction (metadata; see class docs).
    /// @param periodic_x Whether the i-direction is periodic (cyclic).
    /// @param periodic_y Whether the j-direction is periodic (cyclic).
    /// @param tripolar_n Whether the domain has tripolar connectivity (a
    ///        fold) at the northern edge. Not implemented yet; aborts if true.
    /// @param n_boxes Number of boxes to decompose the domain into, or 0
    ///        (the default) for one box per rank. The actual box count can be
    ///        smaller if the domain is too small to split that finely.
    /// @pre The runtime is up (a TIM::Runtime exists); aborts otherwise.
    Domain(int ni_global, int nj_global,
           int ni_halo, int nj_halo,
           bool periodic_x, bool periodic_y,
           bool tripolar_n, int n_boxes = 0);

    /// @brief Global number of cells in the i-direction (x).
    /// @return The global i-extent.
    int ni_global() const { return ni_global_; }

    /// @brief Global number of cells in the j-direction (y).
    /// @return The global j-extent.
    int nj_global() const { return nj_global_; }

    /// @brief Halo width in the i-direction (metadata for field creation).
    /// @return The i-direction halo width.
    int ni_halo() const { return ni_halo_; }

    /// @brief Halo width in the j-direction (metadata for field creation).
    /// @return The j-direction halo width.
    int nj_halo() const { return nj_halo_; }

    /// @brief True if the i-direction is periodic (cyclic).
    /// @return The i-direction periodicity flag.
    bool periodic_x() const { return periodic_x_; }

    /// @brief True if the j-direction is periodic (cyclic).
    /// @return The j-direction periodicity flag.
    bool periodic_y() const { return periodic_y_; }

    /// @brief True if the domain has tripolar connectivity at the northern edge.
    /// @return The tripolar connectivity flag.
    bool tripolar_n() const { return tripolar_n_; }

    /// @brief The cell-centered decomposition of the global domain, extended
    /// to the requested number of vertical levels. Every box spans the full
    /// k-range [0, n_levels): the decomposition is horizontal-only.
    /// @param n_levels Number of vertical levels of the field to be created:
    ///        1 for 2-D fields, NK for 3-D layer fields, NK+1 for interface
    ///        fields, etc. Aborts if not positive.
    /// @return The cell-centered BoxArray with the requested k-extent.
    amrex::BoxArray boxArray(int n_levels) const;

    /// @brief The assignment of the horizontal boxes to MPI ranks.
    /// @return The distribution mapping of the horizontal decomposition.
    const amrex::DistributionMapping& distribution_mapping() const {
        return distribution_mapping_;
    }

    /// @brief The domain's periodicity in index space, for halo exchanges
    /// (e.g. MultiFab::FillBoundary).
    /// @return The global extent in each periodic direction, 0 in each
    ///         non-periodic one.
    amrex::Periodicity periodicity() const;

private:
    int ni_global_;
    int nj_global_;
    int ni_halo_;
    int nj_halo_;
    bool periodic_x_;
    bool periodic_y_;
    bool tripolar_n_;

    /// @brief The horizontal (single-level) decomposition.
    amrex::BoxArray box_array_2d_;
    /// @brief The assignment of the horizontal boxes to MPI ranks.
    amrex::DistributionMapping distribution_mapping_;
};

}  // namespace TIM
