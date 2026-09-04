#pragma once
/**
 * @file mom_continuity_ppm_kernel.hpp
 * @brief Per-cell device primitives for the MOM6 PPM continuity kernels.
 */

#include <AMReX_Box.H>
#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>
#include <AMReX_Math.H>

namespace MOM {
using amrex::Real;
using namespace amrex::literals;
/**
 * @brief Piecewise parabolic limiter
 *
 *  This function limits the left/right edge values of the PPM reconstruction
 *  to give a reconstruction that is positive-definite.  Here this is
 *  reinterpreted as giving a constant thickness if the mean thickness is less
 *  than @p h_min, with a minimum of @p h_min otherwise.
 *
 *  @pre h_in >= 0
 *  @pre h_min > 0
 *
 *  @param h_in Layer thickness [H ~> m or kg m-2].
 *  @param h_L  Left thickness in the reconstruction [H ~> m or kg m-2].
 *  @param h_R Right thickness in the reconstruction [H ~> m or kg m-2].
 *  @param h_min The minimum thickness that can be obtained by a concave
 *              parabolic fit [H ~> m or kg m-2]
 *  @note On return, @p h_L and @p h_R hold the modified thickness values.
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void ppm_limit_pos_point(Real& h_L,
                         Real& h_R,
                         Real const h_in,
                         Real const h_min 
                         ) noexcept
{
    /// This limiter prevents undershooting minima within the domain with
    /// values less than h_min.
    /// The grid-normalized curvature of the three thicknesses  [H ~> m or kg m-2]
    Real const curv = 3.0_rt * ((h_L + h_R) - 2.0_rt * h_in);

    if (curv > 0.0_rt) { /// Only minima are limited.
        Real const dh = h_R - h_L; ///< The difference between the edge thicknesses [H ~> m or kg m-2]

        if (amrex::Math::abs(dh) < curv) { /// The parabola's minimum is within the cell.
            if (h_in <= h_min) {
                h_L = h_in;
                h_R = h_in;
            }
            else if (12.0_rt * curv * (h_in - h_min) < ((curv * curv) + (3.0_rt * (dh * dh)))) {
                /// The minimum value is h_in - (curv^2 + 3*dh^2)/(12*curv), and must
                /// be limited in this case.  0 < scale < 1.
		/// A scaling factor to reduce the curvature of the fit     [nondim]
                Real const scale = 12.0_rt * curv * (h_in - h_min) / ((curv * curv) + (3.0_rt * (dh * dh)));

                h_L = h_in + scale * (h_L - h_in);
                h_R = h_in + scale * (h_R - h_in);
            }
        }
    }
}

/**
 * @brief Peacewise parabolic limiter of Colella and Woodward, 1984
 *
 *  This subroutine limits the left/right edge values of the PPM reconstruction
 *  according to the monotonic prescription of Colella and Woodward, 1984.
 *
 *  @pre h_in >= 0
 *
 *  @param h_in Layer thickness [H ~> m or kg m-2].
 *  @param h_L  Left thickness in the reconstruction [H ~> m or kg m-2].
 *  @param h_R Right thickness in the reconstruction [H ~> m or kg m-2].
 *
 *  @note On return, @p h_L and @p h_R hold the modified thickness values.
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void ppm_limit_cw84_point(Real& h_L,
                          Real& h_R,
                          Real const h_in
			  ) noexcept
{
    /// This limiter monotonizes the parabola following
    /// Colella and Woodward, 1984, Eq. 1.10
    Real h_i = h_in;  ///< A copy of the cell-average layer thickness

    if ( ( h_R - h_i ) * ( h_i - h_L ) <= 0.0_rt ) {
        h_L = h_i;
        h_R = h_i;
    } else {
        Real const RLdiff  = h_R - h_L;           /// The difference between the input edge values
        Real const RLmean  = 0.5_rt * ( h_R + h_L );  /// The average of the input edge thicknesses
        Real const FunFac  = 6.0_rt * RLdiff * ( h_i - RLmean ); /// A curious product of the thickness slope and curvature
        Real const RLdiff2 = RLdiff * RLdiff;  //// The squared difference between the input edge values

        if ( FunFac >  RLdiff2 ) h_L = 3.0_rt * h_i - 2.0_rt * h_R;
        if ( FunFac < -RLdiff2 ) h_R = 3.0_rt * h_i - 2.0_rt * h_L;
    }
}

/**
 * @brief Upwind (1st-order) edge thickness: copy the cell thickness to both edges.
 *
 *  @param h_L  Left edge thickness [H ~> m or kg m-2].
 *  @param h_R  Right edge thickness [H ~> m or kg m-2].
 *  @param h_in Cell-average layer thickness [H ~> m or kg m-2].
 *
 *  @note On return, @p h_L and @p h_R are both set to @p h_in.
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void edge_thickness_upwind_point(Real& h_L, Real& h_R, Real const h_in) noexcept
{
    h_L = h_in;
    h_R = h_in;
}

/**
 * @brief Advances a layer thickness by the convergence of a thickness flux.
 *
 *  Shared by continuity_zonal_convergence (flux_out/flux_in = uh(i,j,k)/uh(i-1,j,k))
 *  and continuity_meridional_convergence (flux_out/flux_in = vh(i,j,k)/vh(i,j-1,k)),
 *  and by both branches of each (h_prev = hin(i,j,k) when present, else h(i,j,k)
 *  itself) -- the stencil read happens at the call site, this primitive only
 *  combines already-indexed scalars.
 *
 *  @param h_prev   Thickness before this update -- hin(i,j,k), or h(i,j,k) itself
 *                  when hin is absent [H ~> m or kg m-2].
 *  @param flux_out Thickness flux out of the cell's high-index face [H L2 T-1 ~> m3 s-1].
 *  @param flux_in  Thickness flux out of the cell's low-index face [H L2 T-1 ~> m3 s-1].
 *  @param dt       Time increment [T ~> s].
 *  @param IareaT   1/areaT for this column [L-2 ~> m-2].
 *  @param h_min    The minimum layer thickness [H ~> m or kg m-2].
 *  @return         The updated layer thickness, floored at @p h_min [H ~> m or kg m-2].
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
Real continuity_convergence_point(Real const h_prev, Real const flux_out, Real const flux_in,
                                  Real const dt, Real const IareaT, Real const h_min) noexcept
{
    return amrex::max(h_prev - dt * IareaT * (flux_out - flux_in), h_min);
}

/**
 * @brief PPM-reconstructed volume/mass transport and its velocity derivative
 *        across one zonal or meridional face, for one candidate face velocity.
 *
 *  @param u             Zonal or meridional velocity [L T-1 ~> m s-1].
 *  @param h             Layer thickness [H ~> m or kg m-2].
 *  @param h_p1          Layer thickness, offset by 1 [H ~> m or kg m-2].
 *  @param h_L           West/South edge thickness [H ~> m or kg m-2].
 *  @param h_L_p1        West/South edge thickness, offset by 1 [H ~> m or kg m-2].
 *  @param h_R           East/North edge thickness [H ~> m or kg m-2].
 *  @param h_R_p1        East/North edge thickness, offset by 1 [H ~> m or kg m-2].
 *  @param uh            Zonal or meridional mass/volume transport [H L2 T-1 ~> m3 s-1 or kg s-1].
 *  @param duhdu         Partial derivative of uh with u [H L ~> m2 or kg m-1].
 *  @param visc_rem      Fraction of momentum/barotropic acceleration remaining after viscosity [nondim].
 *  @param G_dy_Cu       Unblocked u/v-face length of the h-cell [L ~> m].
 *  @param G_IareaT      1/areaT [L-2 ~> m-2].
 *  @param G_IareaT_p1   1/areaT, offset by 1 [L-2 ~> m-2].
 *  @param G_IdxT        1/dxT [L-1 ~> m-1].
 *  @param G_IdxT_p1     1/dxT, offset by 1 [L-1 ~> m-1].
 *  @param dt            Time increment [T ~> s].
 *  @param vol_CFL       If true, rescale the ratio of face areas to cell areas when estimating CFL.
 *  @param por_face_area Fractional open area of the U/V-face [nondim].
 *
 *  @note On return, @p uh and @p duhdu hold the computed transport and derivative.
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void flux_elem_point(Real const u,
                     Real const h,
                     Real const h_p1,
                     Real const h_L,
                     Real const h_L_p1,
                     Real const h_R,
                     Real const h_R_p1,
                     Real& uh,
                     Real& duhdu,
                     Real const visc_rem,
                     Real const G_dy_Cu,
                     Real const G_IareaT,
                     Real const G_IareaT_p1,
                     Real const G_IdxT,
                     Real const G_IdxT_p1,
                     Real const dt,
                     bool const vol_CFL,
                     Real const por_face_area) noexcept
{
    Real const tmp = G_dy_Cu * por_face_area;
    Real CFL, curv_3, h_marg, dh;

    if (u > 0.0_rt) {
        if (vol_CFL) {
            CFL = (u * dt) * (G_dy_Cu * G_IareaT);
        } else {
            CFL = u * dt * G_IdxT;
        }
        curv_3 = (h_L + h_R) - 2.0_rt*h;
        dh = h_L - h_R;
        uh = tmp * u * (h_R + CFL * (0.5_rt*dh + curv_3*(CFL - 1.5_rt)));
        h_marg = h_R + CFL * (dh + 3.0_rt*curv_3*(CFL - 1.0_rt));
    } else if (u < 0.0_rt) {
        if (vol_CFL) {
            CFL = (-u * dt) * (G_dy_Cu * G_IareaT_p1);
        } else {
            CFL = -u * dt * G_IdxT_p1;
        }
        curv_3 = (h_L_p1 + h_R_p1) - 2.0_rt*h_p1;
        dh = h_R_p1 - h_L_p1;
        uh = tmp * u * (h_L_p1 + CFL * (0.5_rt*dh + curv_3*(CFL - 1.5_rt)));
        h_marg = h_L_p1 + CFL * (dh + 3.0_rt*curv_3*(CFL - 1.0_rt));
    } else {
        uh = 0.0_rt;
        h_marg = 0.5_rt * (h_L_p1 + h_R);
    }
    duhdu = tmp * h_marg * visc_rem;
}
}
