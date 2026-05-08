// mom_continuity_ppm.hpp
#pragma once

#include "mom_continuity_ppm_kernel.hpp"

struct OceanOBC;    // Undefined at the moment

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
