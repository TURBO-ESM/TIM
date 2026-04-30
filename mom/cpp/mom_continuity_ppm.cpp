// mom_continuity_ppm.cpp
#define AMREX_ABORT_LOC(msg) \
	amrex::Abort(std::string(msg) + " [" + __FILE__ + ":" + std::to_string(__LINE__) + "]")
#include <AMReX.H>
#include <AMReX_FArrayBox.H>

#include "mom_continuity_ppm.hpp"


namespace MOM {
using amrex::FArrayBox;
/**
 * @brief Piecewise parabolic limiter
 */
void ppm_limit_pos(const Box & bx,
		  Array4<const Real> const& h_in,
		  Array4<Real> const& h_L,
		  Array4<Real> const& h_R,
                  const Real h_min)
{
    BL_PROFILE("ppm_limit_pos");

    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // This limiter prevents undershooting minima within the domain with
        //  values less than h_min.
        ppm_limit_pos_point(h_L(i,j,k), h_R(i,j,k), h_in(i,j,k), h_min);
    });
}

/**
 * @brief Peacewise parabolic limiter of Colella and Woodward, 1984
 */
void ppm_limit_cw84(const Box & bx,
		   Array4<const Real> const& h_in,
		   Array4<Real> const& h_L,
		   Array4<Real> const& h_R)
{
    BL_PROFILE("ppm_limit_cw84");

    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // This limiter monotonizes the parabola following
        // Colella and Woodward, 1984, Eq. 1.10
        ppm_limit_cw84_point(h_L(i,j,k), h_R(i,j,k), h_in(i,j,k));
    });
}

//> Calculates left/right edge values for PPM reconstruction.
void PPM_reconstruction_y(
    const Box& bxH,                 //!< H-grid iteration Box
    const Array4<const Real>& h_in,   //!< Layer thickness
    const Array4<Real>& h_S,          //!< South edge thickness
    const Array4<Real>& h_N,          //!< North edge thickness
    const Array4<const Real>& mask2dT,//!< 0 for land, 1 for ocean
    Real h_min,                     //!< Minimum thickness
    bool monotonic,                       //!< Use CW84 limiter if true
    bool simple_2nd,                      //!< Use simple 2nd order if true
    OceanOBC* OBC                         //!< Open boundary control structure
)
{
    BL_PROFILE("PPM_reconstruction_y");

    // Local variables
    const Real oneSixth = 1.0 / 6.0;

    // NOTE: OBC support temporarily disabled.
    // OceanOBC is forward-declared only.
    // All boundary-condition logic removed for initial port validation.
    if (OBC != nullptr) {
       AMREX_ABORT_LOC("OBC pointer provided but not yet implemented");
    }
    /*
    bool local_open_BC = false;
    if (OBC != nullptr) {
        local_open_BC = OBC->open_v_BCs_exist_globally;
    }
    */

    // Local iteration box extends the h-grid by one element
    Box bx  = grow(bxH, 1, 1);  // grow in y-direction (dim=1)

    // Extended iteration box extends the h-grid by two elements
    Box bxE = grow(bxH, 1, 2); // grow in y-dimension (dim=1)

    // Temporary slope array
    FArrayBox slp_fab(bxE, 1);
    Array4<Real> slp = slp_fab.array(); 

    if (simple_2nd) {

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
	    Real h_jm1 = mask2dT(i,j-1,0) * h_in(i,j-1,k)
                  + (1.0 - mask2dT(i,j-1,0)) * h_in(i,j,k);

	    Real h_jp1 = mask2dT(i,j+1,0) * h_in(i,j+1,k)
                  + (1.0 - mask2dT(i,j+1,0)) * h_in(i,j,k);

            h_S(i,j,k) = 0.5 * (h_jm1 + h_in(i,j,k));
            h_N(i,j,k) = 0.5 * (h_jp1 + h_in(i,j,k));
        });

    } else {

        // Compute slopes on expanded box
        ParallelFor(bxE, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if ((mask2dT(i,j-1,0) * mask2dT(i,j,0) * mask2dT(i,j+1,0)) == 0.0) {
                slp(i,j,k) = 0.0;
            } else {
                // Simple 2nd order slope
	        Real slope = 0.5 * (h_in(i,j+1,k) - h_in(i,j-1,k));

                // Monotonic constraint (Lin 1994, Eq. B2)
		Real dMx = amrex::max(amrex::max(h_in(i,j+1,k), h_in(i,j-1,k)), h_in(i,j,k)) - h_in(i,j,k);
		Real dMn = h_in(i,j,k) - amrex::min(amrex::min(h_in(i,j+1,k), h_in(i,j-1,k)), h_in(i,j,k));

                slp(i,j,k) = amrex::Math::copysign(
                    amrex::min(amrex::Math::abs(slope), 2.0 * amrex::min(dMx, dMn)),
                    slope
                );
            }
        });

	/*
        // Apply open boundary condition to slopes
        if (local_open_BC) {
            for (int n = 0; n < OBC->number_of_segments; ++n) {
                auto& segment = OBC->segment[n];
                if (!segment.on_pe) continue;

                if (segment.direction == OBC_DIRECTION_S ||
                    segment.direction == OBC_DIRECTION_N) {

                    int j = segment.HI.JsdB;

                    ParallelFor(Box(IntVect(segment.HI.isd, j, bx.smallEnd(2)),
                                    IntVect(segment.HI.ied, j, bx.bigEnd(2))),
                    [=] AMREX_GPU_DEVICE (int i, int jj, int k)
                    {
                        slp(i,j+1,k) = 0.0;
                        slp(i,j,k)   = 0.0;
                    });
                }
            }
        }
	*/

        // Compute edge values
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
	    Real h_jm1 = mask2dT(i,j-1,0) * h_in(i,j-1,k)
                  + (1.0 - mask2dT(i,j-1,0)) * h_in(i,j,k);

	    Real h_jp1 = mask2dT(i,j+1,0) * h_in(i,j+1,k)
                  + (1.0 - mask2dT(i,j+1,0)) * h_in(i,j,k);

            // Left/right values (Lin 1994 Eq. B2)
            h_S(i,j,k) = 0.5*(h_jm1 + h_in(i,j,k))
                       + oneSixth*(slp(i,j-1,k) - slp(i,j,k));

            h_N(i,j,k) = 0.5*(h_jp1 + h_in(i,j,k))
                       + oneSixth*(slp(i,j,k) - slp(i,j+1,k));
        });
    }

    /*
    // Apply open boundary condition to final values
    if (local_open_BC) {
        for (int n = 0; n < OBC->number_of_segments; ++n) {
            auto& segment = OBC->segment[n];
            if (!segment.on_pe) continue;

            int j = segment.HI.JsdB;

            if (segment.direction == OBC_DIRECTION_N) {

                ParallelFor(Box(IntVect(segment.HI.isd, j, bx.smallEnd(2)),
                                IntVect(segment.HI.ied, j, bx.bigEnd(2))),
                [=] AMREX_GPU_DEVICE (int i, int jj, int k)
                {
                    h_S(i,j+1,k) = h_in(i,j,k);
                    h_N(i,j+1,k) = h_in(i,j,k);
                    h_S(i,j,k)   = h_in(i,j,k);
                    h_N(i,j,k)   = h_in(i,j,k);
                });

            } else if (segment.direction == OBC_DIRECTION_S) {

                ParallelFor(Box(IntVect(segment.HI.isd, j, bx.smallEnd(2)),
                                IntVect(segment.HI.ied, j, bx.bigEnd(2))),
                [=] AMREX_GPU_DEVICE (int i, int jj, int k)
                {
                    h_S(i,j,k)   = h_in(i,j+1,k);
                    h_N(i,j,k)   = h_in(i,j+1,k);
                    h_S(i,j+1,k) = h_in(i,j+1,k);
                    h_N(i,j+1,k) = h_in(i,j+1,k);
                });
            }
        }
	
    }
    */

    // Apply limiters
    if (monotonic) {
        ppm_limit_cw84(bx, h_in, h_S, h_N);
    } else {
        ppm_limit_pos(bx, h_in, h_S, h_N, h_min);
    }
}
}
