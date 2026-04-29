#pragma once

#include <AMReX_Box.H>
#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>
#include <AMReX_Math.H>

namespace MOM {
using amrex::Real;
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
 *  @return Modified thickness values @p HL and @p HR.
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
    Real const curv = 3.0 * ((h_L + h_R) - 2.0 * h_in);

    if (curv > 0.0) { /// Only minima are limited.
        Real const dh = h_R - h_L; ///< The difference between the edge thicknesses [H ~> m or kg m-2]

        if (amrex::Math::abs(dh) < curv) { /// The parabola's minimum is within the cell.
            if (h_in <= h_min) {
                h_L = h_in;
                h_R = h_in;
            }
            else if (12.0 * curv * (h_in - h_min) < (curv * curv + 3.0 * dh * dh)) {
                /// The minimum value is h_in - (curv^2 + 3*dh^2)/(12*curv), and must
                /// be limited in this case.  0 < scale < 1.
		/// A scaling factor to reduce the curvature of the fit     [nondim]
                Real const scale = 12.0 * curv * (h_in - h_min) / (curv * curv + 3.0 * dh * dh);

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
 *  @return  Modified thickness values for @p h_L and @p h_R.
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

    if ( ( h_R - h_i ) * ( h_i - h_L ) <= 0.0 ) {
        h_L = h_i;
        h_R = h_i;
    } else {
        Real const RLdiff  = h_R - h_L;           /// The difference between the input edge values
        Real const RLmean  = 0.5 * ( h_R + h_L );  /// The average of the input edge thicknesses
        Real const FunFac  = 6.0 * RLdiff * ( h_i - RLmean ); /// A curious product of the thickness slope and curvature
        Real const RLdiff2 = RLdiff * RLdiff;  //// The squared difference between the input edge values

        if ( FunFac >  RLdiff2 ) h_L = 3.0 * h_i - 2.0 * h_R;
        if ( FunFac < -RLdiff2 ) h_R = 3.0 * h_i - 2.0 * h_L;
    }
}
}
