// mom_continuity_ppm.hpp
#pragma once
/**
 * @file mom_continuity_ppm.hpp
 * @brief Box-level AMReX kernel declarations for MOM6 PPM continuity
 *        (piecewise parabolic reconstruction and edge-thickness routines).
 */

#include "mom_continuity_ppm_kernel.hpp"

struct OceanOBC;    // Undefined at the moment

/// @brief Options controlling the transport adjustment and barotropic-consistency
/// iteration used by the continuity solver. Field-for-field mirror of the Fortran
/// `bind(C)` type `transport_adjust_CS_C` -- order and types must not change.
struct transport_adjust_CS_C {
    double tol_eta;            ///< Tolerance for free-surface height discrepancies.
    double tol_vel;            ///< Tolerance for barotropic velocity discrepancies.
    double CFL_limit_adjust;   ///< Maximum CFL of the adjusted velocities.
    bool   aggress_adjust;     ///< If true, allow a larger relative CFL change.
    bool   vol_CFL;            ///< If true, use the ratio of open face lengths to
                                ///< tracer cell areas when estimating CFL numbers.
    bool   better_iter;        ///< If true, use a velocity-based iteration criterion.
    bool   use_visc_rem_max;   ///< If true, use limiting bounds for viscous columns.
    bool   marginal_faces;     ///< If true, use marginal face areas as barotropic weights.
};

/// @brief AMReX ports of MOM6 numerical kernels.
namespace MOM {
using amrex::Box;
using amrex::Array4;
/**
 * @brief Piecewise parabolic limiter
 */
void ppm_limit_pos(const Box &,
                   Array4<const Real> const&,
                   Array4<Real> const&,
                   Array4<Real> const&,
                   const Real);

/**
 * @brief Piecewise parabolic limiter of Colella and Woodward, 1984
 */
void ppm_limit_cw84(const Box&,
                    Array4<const Real> const&,
                    Array4<Real> const&,
                    Array4<Real> const&);

/**
 * @brief Piecewise reconstruction in the y dimension
 */
void PPM_reconstruction_y(
    const Box&,
    Array4<const Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<const Real> const&,
    Real,
    bool,
    bool,
    OceanOBC*);

/**
 * @brief Piecewise reconstruction in the x dimension
 */
void PPM_reconstruction_x(
    const Box&,
    Array4<const Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<const Real> const&,
    Real,
    bool,
    bool,
    OceanOBC*);

/**
 * @brief Zonal edge thickness — upwind copy or x-direction PPM reconstruction
 */
void zonal_edge_thickness(
    const Box&,
    Array4<const Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<const Real> const&,
    Real,
    bool,
    bool,
    bool,
    OceanOBC*);

/**
 * @brief Meridional edge thickness — upwind copy or y-direction PPM reconstruction
 */
void meridional_edge_thickness(
    const Box&,
    Array4<const Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<const Real> const&,
    Real,
    bool,
    bool,
    bool,
    OceanOBC*);
}
