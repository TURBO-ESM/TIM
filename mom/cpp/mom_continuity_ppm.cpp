// mom_continuity_ppm.cpp
/**
 * @file mom_continuity_ppm.cpp
 * @brief Box-level AMReX kernel implementations for MOM6 PPM continuity.
 */
/// @brief Abort with @p msg annotated by the source file and line number.
// SKILLS: 0.3.1
#define AMREX_ABORT_LOC(msg) \
	amrex::Abort(std::string(msg) + " [" + __FILE__ + ":" + std::to_string(__LINE__) + "]")
#include <AMReX.H>
#include <AMReX_FArrayBox.H>

#include "mom_continuity_ppm.hpp"


namespace MOM {
using amrex::FArrayBox;
using namespace amrex::literals;
/**
 * @brief Piecewise parabolic limiter (positive-definite) over a Box.
 *
 * @param bx    Iteration Box.
 * @param h_in  Layer thickness [H ~> m or kg m-2].
 * @param h_L   Left edge thickness of the reconstruction [H ~> m or kg m-2].
 * @param h_R   Right edge thickness of the reconstruction [H ~> m or kg m-2].
 * @param h_min Minimum thickness allowed by the parabolic fit [H ~> m or kg m-2].
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
 * @brief Piecewise parabolic limiter of Colella and Woodward, 1984, over a Box.
 *
 * @param bx   Iteration Box.
 * @param h_in Layer thickness [H ~> m or kg m-2].
 * @param h_L  Left edge thickness of the reconstruction [H ~> m or kg m-2].
 * @param h_R  Right edge thickness of the reconstruction [H ~> m or kg m-2].
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
    const Real oneSixth = 1.0_rt / 6.0_rt;

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
                  + (1.0_rt - mask2dT(i,j-1,0)) * h_in(i,j,k);

	    Real h_jp1 = mask2dT(i,j+1,0) * h_in(i,j+1,k)
                  + (1.0_rt - mask2dT(i,j+1,0)) * h_in(i,j,k);

            h_S(i,j,k) = 0.5_rt * (h_jm1 + h_in(i,j,k));
            h_N(i,j,k) = 0.5_rt * (h_jp1 + h_in(i,j,k));
        });

    } else {

        // Compute slopes on expanded box
        ParallelFor(bxE, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if ((mask2dT(i,j-1,0) * mask2dT(i,j,0) * mask2dT(i,j+1,0)) == 0.0_rt) {
                slp(i,j,k) = 0.0_rt;
            } else {
                // Simple 2nd order slope
	        Real slope = 0.5_rt * (h_in(i,j+1,k) - h_in(i,j-1,k));

                // Monotonic constraint (Lin 1994, Eq. B2)
		Real dMx = amrex::max(amrex::max(h_in(i,j+1,k), h_in(i,j-1,k)), h_in(i,j,k)) - h_in(i,j,k);
		Real dMn = h_in(i,j,k) - amrex::min(amrex::min(h_in(i,j+1,k), h_in(i,j-1,k)), h_in(i,j,k));

                slp(i,j,k) = amrex::Math::copysign(
                    amrex::min(amrex::Math::abs(slope), 2.0_rt * amrex::min(dMx, dMn)),
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
                        slp(i,j+1,k) = 0.0_rt;
                        slp(i,j,k)   = 0.0_rt;
                    });
                }
            }
        }
	*/

        // Compute edge values
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
	    Real h_jm1 = mask2dT(i,j-1,0) * h_in(i,j-1,k)
                  + (1.0_rt - mask2dT(i,j-1,0)) * h_in(i,j,k);

	    Real h_jp1 = mask2dT(i,j+1,0) * h_in(i,j+1,k)
                  + (1.0_rt - mask2dT(i,j+1,0)) * h_in(i,j,k);

            // Left/right values (Lin 1994 Eq. B2)
            h_S(i,j,k) = 0.5_rt*(h_jm1 + h_in(i,j,k))
                       + oneSixth*(slp(i,j-1,k) - slp(i,j,k));

            h_N(i,j,k) = 0.5_rt*(h_jp1 + h_in(i,j,k))
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

//> Calculates west/east edge values for PPM reconstruction.
void PPM_reconstruction_x(
    const Box& bxH,                  //!< H-grid iteration Box
    const Array4<const Real>& h_in,  //!< Layer thickness
    const Array4<Real>& h_W,         //!< West edge thickness
    const Array4<Real>& h_E,         //!< East edge thickness
    const Array4<const Real>& mask2dT,//!< 0 for land, 1 for ocean
    Real h_min,                      //!< Minimum thickness
    bool monotonic,                  //!< Use CW84 limiter if true
    bool simple_2nd,                 //!< Use simple 2nd order if true
    OceanOBC* OBC                    //!< Open boundary control structure
)
{
    BL_PROFILE("PPM_reconstruction_x");
    const Real oneSixth = 1.0_rt / 6.0_rt;

    // NOTE: OBC support temporarily disabled.
    if (OBC != nullptr) {
       AMREX_ABORT_LOC("OBC pointer provided but not yet implemented");
    }

    // Local iteration box extends the h-grid by one element in x
    Box bx  = grow(bxH, 0, 1);  // grow in x-direction (dim=0)

    // Extended iteration box extends the h-grid by two elements in x
    Box bxE = grow(bxH, 0, 2);  // grow in x-direction (dim=0)

    // Temporary slope array
    FArrayBox slp_fab(bxE, 1);
    Array4<Real> slp = slp_fab.array();

    if (simple_2nd) {

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real h_im1 = mask2dT(i-1,j,0) * h_in(i-1,j,k)
                      + (1.0_rt - mask2dT(i-1,j,0)) * h_in(i,j,k);
            Real h_ip1 = mask2dT(i+1,j,0) * h_in(i+1,j,k)
                      + (1.0_rt - mask2dT(i+1,j,0)) * h_in(i,j,k);
            h_W(i,j,k) = 0.5_rt * (h_im1 + h_in(i,j,k));
            h_E(i,j,k) = 0.5_rt * (h_ip1 + h_in(i,j,k));
        });

    } else {

        // Compute slopes on expanded box
        ParallelFor(bxE, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if ((mask2dT(i-1,j,0) * mask2dT(i,j,0) * mask2dT(i+1,j,0)) == 0.0_rt) {
                slp(i,j,k) = 0.0_rt;
            } else {
                Real slope = 0.5_rt * (h_in(i+1,j,k) - h_in(i-1,j,k));
                Real dMx = amrex::max(amrex::max(h_in(i+1,j,k), h_in(i-1,j,k)), h_in(i,j,k)) - h_in(i,j,k);
                Real dMn = h_in(i,j,k) - amrex::min(amrex::min(h_in(i+1,j,k), h_in(i-1,j,k)), h_in(i,j,k));
                slp(i,j,k) = amrex::Math::copysign(
                    amrex::min(amrex::Math::abs(slope), 2.0_rt * amrex::min(dMx, dMn)),
                    slope
                );
            }
        });

        // Compute edge values (Lin 1994 Eq. B2)
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real h_im1 = mask2dT(i-1,j,0) * h_in(i-1,j,k)
                      + (1.0_rt - mask2dT(i-1,j,0)) * h_in(i,j,k);
            Real h_ip1 = mask2dT(i+1,j,0) * h_in(i+1,j,k)
                      + (1.0_rt - mask2dT(i+1,j,0)) * h_in(i,j,k);
            h_W(i,j,k) = 0.5_rt * (h_im1 + h_in(i,j,k))
                       + oneSixth * (slp(i-1,j,k) - slp(i,j,k));
            h_E(i,j,k) = 0.5_rt * (h_ip1 + h_in(i,j,k))
                       + oneSixth * (slp(i,j,k) - slp(i+1,j,k));
        });
    }

    // Apply limiters
    if (monotonic) {
        ppm_limit_cw84(bx, h_in, h_W, h_E);
    } else {
        ppm_limit_pos(bx, h_in, h_W, h_E, h_min);
    }
}
//> Zonal edge thickness: upwind copy or x-direction PPM reconstruction.
void zonal_edge_thickness(
    const Box& bxC,                   //!< Continuity iteration box
    const Array4<const Real>& h_in,   //!< Layer thickness
    const Array4<Real>& h_W,          //!< West edge thickness
    const Array4<Real>& h_E,          //!< East edge thickness
    const Array4<const Real>& mask2dT,//!< 0 for land, 1 for ocean
    Real h_min,                       //!< Minimum thickness
    bool upwind_1st,                  //!< If true, use 1st-order upwind
    bool monotonic,                   //!< Use CW84 limiter if true
    bool simple_2nd,                  //!< Use simple 2nd order if true
    OceanOBC* obc                     //!< Open boundary control structure
)
{
    BL_PROFILE("zonal_edge_thickness");
    if (upwind_1st) {
        // 1st-order upwind: set both edges to cell-centre value over box grown by 1 in x
        Box bx = grow(bxC, 0, 1);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            edge_thickness_upwind_point(h_W(i,j,k), h_E(i,j,k), h_in(i,j,k));
        });
    } else {
        PPM_reconstruction_x(bxC, h_in, h_W, h_E, mask2dT, h_min, monotonic, simple_2nd, obc);
    }
}

//> Meridional edge thickness: upwind copy or y-direction PPM reconstruction.
void meridional_edge_thickness(
    const Box& bxC,                   //!< Continuity iteration box
    const Array4<const Real>& h_in,   //!< Layer thickness
    const Array4<Real>& h_S,          //!< South edge thickness
    const Array4<Real>& h_N,          //!< North edge thickness
    const Array4<const Real>& mask2dT,//!< 0 for land, 1 for ocean
    Real h_min,                       //!< Minimum thickness
    bool upwind_1st,                  //!< If true, use 1st-order upwind
    bool monotonic,                   //!< Use CW84 limiter if true
    bool simple_2nd,                  //!< Use simple 2nd order if true
    OceanOBC* obc                     //!< Open boundary control structure
)
{
    BL_PROFILE("meridional_edge_thickness");
    if (upwind_1st) {
        // 1st-order upwind: set both edges to cell-centre value over box grown by 1 in y
        Box bx = grow(bxC, 1, 1);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            edge_thickness_upwind_point(h_S(i,j,k), h_N(i,j,k), h_in(i,j,k));
        });
    } else {
        PPM_reconstruction_y(bxC, h_in, h_S, h_N, mask2dT, h_min, monotonic, simple_2nd, obc);
    }
}

//> Zonal continuity update: advances layer thickness by the convergence
//  of the zonal thickness flux.
void continuity_zonal_convergence(
    const Box& bxC,                   //!< Iteration box for continuity solver
    Array4<Real> const& h,            //!< Final layer thickness
    Array4<const Real> const& uh,     //!< Zonal thickness flux, u*h*dy
    Real dt,                          //!< Time increment
    Array4<const Real> const& IareaT, //!< 1/areaT (2D, addressed at k=0)
    Array4<const Real> const& hin,    //!< Initial layer thickness; may be absent (.p == nullptr),
                                       //!< in which case the final thickness is also the initial
                                       //!< thickness
    Real h_min)                       //!< The minimum layer thickness
{
    BL_PROFILE("continuity_zonal_convergence");

    if (hin.p != nullptr) {
        ParallelFor(bxC, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            h(i,j,k) = continuity_convergence_point(hin(i,j,k), uh(i,j,k), uh(i-1,j,k), dt, IareaT(i,j,0), h_min);
        });
    } else {
        // untested
        ParallelFor(bxC, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            h(i,j,k) = continuity_convergence_point(h(i,j,k), uh(i,j,k), uh(i-1,j,k), dt, IareaT(i,j,0), h_min);
        });
    }
}

//> Meridional continuity update: advances layer thickness by the
//  convergence of the meridional thickness flux.
void continuity_meridional_convergence(
    const Box& bxC,                   //!< Iteration box for continuity solver
    Array4<Real> const& h,            //!< Final layer thickness
    Array4<const Real> const& vh,     //!< Meridional thickness flux, v*h*dx
    Real dt,                          //!< Time increment
    Array4<const Real> const& IareaT, //!< 1/areaT (2D, addressed at k=0)
    Array4<const Real> const& hin,    //!< Initial layer thickness; may be absent (.p == nullptr),
                                       //!< in which case the final thickness is also the initial
                                       //!< thickness
    Real h_min)                       //!< The minimum layer thickness
{
    BL_PROFILE("continuity_meridional_convergence");

    if (hin.p != nullptr) {
        // untested
        ParallelFor(bxC, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            h(i,j,k) = continuity_convergence_point(hin(i,j,k), vh(i,j,k), vh(i,j-1,k), dt, IareaT(i,j,0), h_min);
        });
    } else {
        ParallelFor(bxC, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            h(i,j,k) = continuity_convergence_point(h(i,j,k), vh(i,j,k), vh(i,j-1,k), dt, IareaT(i,j,0), h_min);
        });
    }
}

//> Zonal volume/thickness flux: PPM-reconstructed edge thickness advected
//  by the zonal velocity, scaled by viscosity remnant and open-face area.
void zonal_flux_thickness(
    const Box& bxC,                        //!< Continuity iteration box
    Array4<const Real> const& u,           //!< Zonal velocity
    Array4<const Real> const& h,           //!< Layer thickness used to calculate fluxes
    Array4<const Real> const& h_W,         //!< West edge thickness in the reconstruction
    Array4<const Real> const& h_E,         //!< East edge thickness in the reconstruction
    Array4<Real> const& h_u,               //!< Effective thickness at zonal faces
    Real dt,                               //!< Time increment
    Array4<const Real> const& dy_Cu,       //!< Unblocked u-face length (2D, addressed at k=0)
    Array4<const Real> const& IareaT,      //!< 1/areaT (2D, addressed at k=0)
    Array4<const Real> const& IdxT,        //!< 1/dxT (2D, addressed at k=0)
    bool vol_CFL,                          //!< If true, rescale face/cell area ratio for CFL
    bool marginal,                         //!< If true, report marginal (not transport-averaged) thickness
    OceanOBC* OBC,                         //!< Open boundary control structure
    Array4<const Real> const& por_face_areaU, //!< Fractional open area of U-faces
    Array4<const Real> const& visc_rem_u)  //!< Both the fraction of the momentum originally in a layer
                                            //!< that remains after a time-step of viscosity, and the
                                            //!< fraction of a time-step's worth of a barotropic
                                            //!< acceleration that a layer experiences after viscosity is
                                            //!< applied [nondim]. Between 0 (at the bottom) and 1 (far
                                            //!< above the bottom).
{
    BL_PROFILE("zonal_flux_thickness");

    // NOTE: OBC support temporarily disabled.
    // OceanOBC is forward-declared only.
    // All boundary-condition logic removed for initial port validation.
    if (OBC != nullptr) {
       AMREX_ABORT_LOC("OBC pointer provided but not yet implemented");
    }
    /*
    bool local_open_BC = false;
    if (OBC != nullptr) local_open_BC = OBC->open_u_BCs_exist_globally;
    */

    const bool has_visc_rem_u = (visc_rem_u.p != nullptr);

    // Increase the lower extent of the x-dimension (U-grid)
    Box bxU = growLo(bxC, 0, 1);

    ParallelFor(bxU, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real CFL, curv_3, dh;
        if (u(i,j,k) > 0.0_rt) {
            if (vol_CFL) {
                CFL = (u(i,j,k) * dt) * (dy_Cu(i,j,0) * IareaT(i,j,0));
            } else {
                CFL = u(i,j,k) * dt * IdxT(i,j,0);
            }
            curv_3 = (h_W(i,j,k) + h_E(i,j,k)) - 2.0_rt*h(i,j,k);
            dh = h_W(i,j,k) - h_E(i,j,k);
            if (marginal) {
                h_u(i,j,k) = h_E(i,j,k) + CFL * (dh + 3.0_rt*curv_3*(CFL - 1.0_rt));
            } else {
                h_u(i,j,k) = h_E(i,j,k) + CFL * (0.5_rt*dh + curv_3*(CFL - 1.5_rt));
            }
        } else if (u(i,j,k) < 0.0_rt) {
            if (vol_CFL) {
                CFL = (-u(i,j,k) * dt) * (dy_Cu(i,j,0) * IareaT(i+1,j,0));
            } else {
                CFL = -u(i,j,k) * dt * IdxT(i+1,j,0);
            }
            curv_3 = (h_W(i+1,j,k) + h_E(i+1,j,k)) - 2.0_rt*h(i+1,j,k);
            dh = h_E(i+1,j,k) - h_W(i+1,j,k);
            if (marginal) {
                h_u(i,j,k) = h_W(i+1,j,k) + CFL * (dh + 3.0_rt*curv_3*(CFL - 1.0_rt));
            } else {
                h_u(i,j,k) = h_W(i+1,j,k) + CFL * (0.5_rt*dh + curv_3*(CFL - 1.5_rt));
            }
        } else {
            // The choice to use the arithmetic mean here is somewhat arbitrarily, but
            // it should be noted that h_W(i+1,j,k) and h_E(i,j,k) are usually the same.
            h_u(i,j,k) = 0.5_rt * (h_W(i+1,j,k) + h_E(i,j,k));
            // real h_marg = (2.0 * h_W(i+1,j,k) * h_E(i,j,k)) /
            //               (h_W(i+1,j,k) + h_E(i,j,k) + GV%H_subroundoff)
        }

        if (has_visc_rem_u) {
            // Scale back the thickness to account for the effects of viscosity and the
            // fractional open thickness to give an appropriate non-normalized weight for
            // each layer in determining the barotropic acceleration.
            h_u(i,j,k) = h_u(i,j,k) * (visc_rem_u(i,j,k) * por_face_areaU(i,j,k));
        } else {
            h_u(i,j,k) = h_u(i,j,k) * por_face_areaU(i,j,k);
        }
    });

    /*
    // untested
    if (local_open_BC) {
        for (int n = 0; n < OBC->number_of_segments; ++n) {
            auto& segment = OBC->segment[n];
            if (segment.open && segment.is_E_or_W) {
                int i = segment.HI.IsdB;
                Box bxSeg(IntVect(i, segment.HI.jsd, bxU.smallEnd(2)),
                          IntVect(i, segment.HI.jed, bxU.bigEnd(2)));
                if (segment.direction == OBC_DIRECTION_E) {
                    if (has_visc_rem_u) {
                        ParallelFor(bxSeg, [=] AMREX_GPU_DEVICE (int ii, int j, int k)
                        {
                            h_u(i,j,k) = h(i,j,k) * (visc_rem_u(i,j,k) * por_face_areaU(i,j,k));
                        });
                    } else {
                        ParallelFor(bxSeg, [=] AMREX_GPU_DEVICE (int ii, int j, int k)
                        {
                            h_u(i,j,k) = h(i,j,k) * por_face_areaU(i,j,k);
                        });
                    }
                } else {
                    if (has_visc_rem_u) {
                        ParallelFor(bxSeg, [=] AMREX_GPU_DEVICE (int ii, int j, int k)
                        {
                            h_u(i,j,k) = h(i+1,j,k) * (visc_rem_u(i,j,k) * por_face_areaU(i,j,k));
                        });
                    } else {
                        ParallelFor(bxSeg, [=] AMREX_GPU_DEVICE (int ii, int j, int k)
                        {
                            h_u(i,j,k) = h(i+1,j,k) * por_face_areaU(i,j,k);
                        });
                    }
                }
            }
        }
    }
    */
}

//> Meridional volume/thickness flux: PPM-reconstructed edge thickness
//  advected by the meridional velocity, scaled by viscosity remnant and
//  open-face area.
void meridional_flux_thickness(
    const Box& bxC,                        //!< Continuity iteration box
    Array4<const Real> const& v,           //!< Meridional velocity
    Array4<const Real> const& h,           //!< Layer thickness used to calculate fluxes
    Array4<const Real> const& h_S,         //!< South edge thickness in the reconstruction
    Array4<const Real> const& h_N,         //!< North edge thickness in the reconstruction
    Array4<Real> const& h_v,               //!< Effective thickness at meridional faces [out]
    Real dt,                               //!< Time increment
    Array4<const Real> const& dx_Cv,       //!< Unblocked v-face length (2D, addressed at k=0)
    Array4<const Real> const& IareaT,      //!< 1/areaT (2D, addressed at k=0)
    Array4<const Real> const& IdyT,        //!< 1/dyT (2D, addressed at k=0)
    bool vol_CFL,                          //!< If true, rescale face/cell area ratio for CFL
    bool marginal,                         //!< If true, report marginal (not transport-averaged) thickness
    OceanOBC* OBC,                         //!< Open boundary control structure
    Array4<const Real> const& por_face_areaV, //!< Fractional open area of V-faces
    Array4<const Real> const& visc_rem_v)  //!< Both the fraction of the momentum originally in a layer
                                            //!< that remains after a time-step of viscosity, and the
                                            //!< fraction of a time-step's worth of a barotropic
                                            //!< acceleration that a layer experiences after viscosity is
                                            //!< applied [nondim]. Between 0 (at the bottom) and 1 (far
                                            //!< above the bottom).
{
    BL_PROFILE("meridional_flux_thickness");

    // NOTE: OBC support temporarily disabled.
    // OceanOBC is forward-declared only.
    // All boundary-condition logic removed for initial port validation.
    if (OBC != nullptr) {
       AMREX_ABORT_LOC("OBC pointer provided but not yet implemented");
    }
    /*
    bool local_open_BC = false;
    if (OBC != nullptr) local_open_BC = OBC->open_v_BCs_exist_globally;
    */

    const bool has_visc_rem_v = (visc_rem_v.p != nullptr);

    // Increase the lower extent of the y-dimension (V-grid)
    Box bxV = growLo(bxC, 1, 1);

    ParallelFor(bxV, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real CFL, curv_3, dh;
        if (v(i,j,k) > 0.0_rt) {
            if (vol_CFL) {
                CFL = (v(i,j,k) * dt) * (dx_Cv(i,j,0) * IareaT(i,j,0));
            } else {
                CFL = v(i,j,k) * dt * IdyT(i,j,0);
            }
            curv_3 = (h_S(i,j,k) + h_N(i,j,k)) - 2.0_rt*h(i,j,k);
            dh = h_S(i,j,k) - h_N(i,j,k);
            if (marginal) {
                h_v(i,j,k) = h_N(i,j,k) + CFL * (dh + 3.0_rt*curv_3*(CFL - 1.0_rt));
            } else {
                h_v(i,j,k) = h_N(i,j,k) + CFL * (0.5_rt*dh + curv_3*(CFL - 1.5_rt));
            }
        } else if (v(i,j,k) < 0.0_rt) {
            if (vol_CFL) {
                CFL = (-v(i,j,k) * dt) * (dx_Cv(i,j,0) * IareaT(i,j+1,0));
            } else {
                CFL = -v(i,j,k) * dt * IdyT(i,j+1,0);
            }
            curv_3 = (h_S(i,j+1,k) + h_N(i,j+1,k)) - 2.0_rt*h(i,j+1,k);
            dh = h_N(i,j+1,k) - h_S(i,j+1,k);
            if (marginal) {
                h_v(i,j,k) = h_S(i,j+1,k) + CFL * (dh + 3.0_rt*curv_3*(CFL - 1.0_rt));
            } else {
                h_v(i,j,k) = h_S(i,j+1,k) + CFL * (0.5_rt*dh + curv_3*(CFL - 1.5_rt));
            }
        } else {
            // The choice to use the arithmetic mean here is somewhat arbitrarily, but
            // it should be noted that h_S(i+1,j,k) and h_N(i,j,k) are usually the same.
            h_v(i,j,k) = 0.5_rt * (h_S(i,j+1,k) + h_N(i,j,k));
            // real h_marg = (2.0 * h_S(i,j+1,k) * h_N(i,j,k)) /
            //               (h_S(i,j+1,k) + h_N(i,j,k) + GV%H_subroundoff)
        }

        if (has_visc_rem_v) {
            // Scale back the thickness to account for the effects of viscosity and the
            // fractional open thickness to give an appropriate non-normalized weight for
            // each layer in determining the barotropic acceleration.
            h_v(i,j,k) = h_v(i,j,k) * (visc_rem_v(i,j,k) * por_face_areaV(i,j,k));
        } else {
            h_v(i,j,k) = h_v(i,j,k) * por_face_areaV(i,j,k);
        }
    });

    /*
    // untested - will need to be refactored to be performant on GPUs
    if (local_open_BC) {
        for (int n = 0; n < OBC->number_of_segments; ++n) {
            auto& segment = OBC->segment[n];
            if (segment.open && segment.is_N_or_S) {
                int j = segment.HI.JsdB;
                Box bxSeg(IntVect(segment.HI.isd, j, bxV.smallEnd(2)),
                          IntVect(segment.HI.ied, j, bxV.bigEnd(2)));
                if (segment.direction == OBC_DIRECTION_N) {
                    if (has_visc_rem_v) {
                        ParallelFor(bxSeg, [=] AMREX_GPU_DEVICE (int i, int jj, int k)
                        {
                            h_v(i,j,k) = h(i,j,k) * (visc_rem_v(i,j,k) * por_face_areaV(i,j,k));
                        });
                    } else {
                        ParallelFor(bxSeg, [=] AMREX_GPU_DEVICE (int i, int jj, int k)
                        {
                            h_v(i,j,k) = h(i,j,k) * por_face_areaV(i,j,k);
                        });
                    }
                } else {
                    if (has_visc_rem_v) {
                        ParallelFor(bxSeg, [=] AMREX_GPU_DEVICE (int i, int jj, int k)
                        {
                            h_v(i,j,k) = h(i,j+1,k) * (visc_rem_v(i,j,k) * por_face_areaV(i,j,k));
                        });
                    } else {
                        ParallelFor(bxSeg, [=] AMREX_GPU_DEVICE (int i, int jj, int k)
                        {
                            h_v(i,j,k) = h(i,j+1,k) * por_face_areaV(i,j,k);
                        });
                    }
                }
            }
        }
    }
    */
}

//> Sets the effective open face areas and barotropic velocity corrections at
//  zonal faces that reproduce the summed layer transports for three
//  test barotropic velocities (westerly, easterly, and zero-adjustment),
//  for use in the barotropic-consistency iteration.
void set_zonal_BT_cont(
    const Box& bxC,                          //!< Iteration box for continuity solver
    Array4<const Real> const& u,             //!< Zonal velocity
    Array4<const Real> const& h_in,          //!< Layer thickness used to calculate fluxes
    Array4<const Real> const& h_W,           //!< West edge thickness in the reconstruction
    Array4<const Real> const& h_E,           //!< East edge thickness in the reconstruction
    Array4<Real> const& FA_u_W0,             //!< Effective open face area, west, 0 transport
    Array4<Real> const& FA_u_E0,             //!< Effective open face area, east, 0 transport
    Array4<Real> const& FA_u_WW,             //!< Effective open face area, westerly test velocity
    Array4<Real> const& FA_u_EE,             //!< Effective open face area, easterly test velocity
    Array4<Real> const& uBT_WW,              //!< Westerly correction to the barotropic velocity
    Array4<Real> const& uBT_EE,              //!< Easterly correction to the barotropic velocity
    Array4<const Real> const& du0,           //!< Barotropic velocity increment that gives 0 transport
    Array4<const Real> const& uh_tot_0,      //!< Summed transport with 0 adjustment
    Array4<const Real> const& duhdu_tot_0,   //!< Partial derivative of du_err with du at 0 adjustment
    Array4<const Real> const& du_max_CFL,    //!< Maximum acceptable value of du
    Array4<const Real> const& du_min_CFL,    //!< Minimum acceptable value of du
    Real dt,                                 //!< Time increment
    Array4<const Real> const& dxCu,          //!< The grid cell's u-point x-extent
    Array4<const Real> const& dy_Cu,         //!< Unblocked u-face length (2D, addressed at k=0)
    Array4<const Real> const& IareaT,        //!< 1/areaT (2D, addressed at k=0)
    Array4<const Real> const& IdxT,          //!< 1/dxT (2D, addressed at k=0)
    const transport_adjust_CS_C& CS,         //!< Transport-adjustment and barotropic-consistency options
    Array4<const Real> const& visc_rem,      //!< Fraction of momentum/barotropic acceleration
                                              //!< remaining after viscosity
    Array4<const Real> const& visc_rem_max,  //!< Maximum allowable viscosity remnant
    Array4<const int> const& do_I,           //!< Logical flag (0/1) indicating which I values to work on
    Array4<const Real> const& por_face_areaU)//!< Fractional open area of U-faces
{
    BL_PROFILE("set_zonal_BT_cont");

    const Real Idt = 1.0_rt / dt;
    const Real min_visc_rem = 0.1_rt;
    const Real CFL_min = 1.0e-6_rt;

    const int kmin = bxC.smallEnd(2);
    const int kmax = bxC.bigEnd(2);

    // Iteration box for u-point (U-grid) fields: grown by 1 at the lower x-extent
    Box bxU = growLo(bxC, 0, 1);
    Box bx2d(IntVect(bxU.smallEnd(0), bxU.smallEnd(1), 0),
             IntVect(bxU.bigEnd(0),   bxU.bigEnd(1),   0));

    ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        const bool active = (do_I(i,j,0) != 0);

        // Determine the westerly- and easterly- fluxes. Choose a sufficiently
        // negative velocity correction for the easterly-flux, and a sufficiently
        // positive correction for the westerly-flux.
        const Real du_CFL = (CFL_min * Idt) * dxCu(i,j,0);
        Real duR = amrex::min(0.0_rt, du0(i,j,0) - du_CFL);
        Real duL = amrex::max(0.0_rt, du0(i,j,0) + du_CFL);
        Real FAmt_L = 0.0_rt, FAmt_R = 0.0_rt, FAmt_0 = 0.0_rt;
        Real uhtot_L = 0.0_rt, uhtot_R = 0.0_rt;

        if (active) {
            for (int k = kmin; k <= kmax; ++k) {
                Real const visc_rem_lim = amrex::max(visc_rem(i,j,k), min_visc_rem*visc_rem_max(i,j,0));
                if (visc_rem_lim > 0.0_rt) { // This is almost always true for ocean points.
                    if (u(i,j,k) + duR*visc_rem_lim > -du_CFL*visc_rem(i,j,k)) {
                        duR = -(u(i,j,k) + du_CFL*visc_rem(i,j,k)) / visc_rem_lim;
                    }
                    if (u(i,j,k) + duL*visc_rem_lim < du_CFL*visc_rem(i,j,k)) {
                        duL = -(u(i,j,k) - du_CFL*visc_rem(i,j,k)) / visc_rem_lim;
                    }
                }
            }

            for (int k = kmin; k <= kmax; ++k) {
                Real const u_L = u(i,j,k) + duL * visc_rem(i,j,k);
                Real const u_R = u(i,j,k) + duR * visc_rem(i,j,k);
                Real const u_0 = u(i,j,k) + du0(i,j,0) * visc_rem(i,j,k);
                Real uh_0, uh_L, uh_R, duhdu_0, duhdu_L, duhdu_R;
                flux_elem_point(u_0, h_in(i,j,k), h_in(i+1,j,k), h_W(i,j,k), h_W(i+1,j,k),
                                h_E(i,j,k), h_E(i+1,j,k), uh_0, duhdu_0, visc_rem(i,j,k),
                                dy_Cu(i,j,0), IareaT(i,j,0), IareaT(i+1,j,0), IdxT(i,j,0), IdxT(i+1,j,0),
                                dt, CS.vol_CFL, por_face_areaU(i,j,k));
                flux_elem_point(u_L, h_in(i,j,k), h_in(i+1,j,k), h_W(i,j,k), h_W(i+1,j,k),
                                h_E(i,j,k), h_E(i+1,j,k), uh_L, duhdu_L, visc_rem(i,j,k),
                                dy_Cu(i,j,0), IareaT(i,j,0), IareaT(i+1,j,0), IdxT(i,j,0), IdxT(i+1,j,0),
                                dt, CS.vol_CFL, por_face_areaU(i,j,k));
                flux_elem_point(u_R, h_in(i,j,k), h_in(i+1,j,k), h_W(i,j,k), h_W(i+1,j,k),
                                h_E(i,j,k), h_E(i+1,j,k), uh_R, duhdu_R, visc_rem(i,j,k),
                                dy_Cu(i,j,0), IareaT(i,j,0), IareaT(i+1,j,0), IdxT(i,j,0), IdxT(i+1,j,0),
                                dt, CS.vol_CFL, por_face_areaU(i,j,k));
                FAmt_0 += duhdu_0;
                FAmt_L += duhdu_L;
                FAmt_R += duhdu_R;
                uhtot_L += uh_L;
                uhtot_R += uh_R;
            }

            Real FA_0 = FAmt_0, FA_avg = FAmt_0;
            if ((duL - du0(i,j,0)) != 0.0_rt) {
                FA_avg = uhtot_L / (duL - du0(i,j,0));
            }
            if (FA_avg > amrex::max(FA_0, FAmt_L)) {
                FA_avg = amrex::max(FA_0, FAmt_L);
            } else if (FA_avg < amrex::min(FA_0, FAmt_L)) {
                FA_0 = FA_avg;
            }

            FA_u_W0(i,j,0) = FA_0; FA_u_WW(i,j,0) = FAmt_L;
            if (amrex::Math::abs(FA_0 - FAmt_L) <= 1.0e-12_rt*FA_0) {
                uBT_WW(i,j,0) = 0.0_rt;
            } else {
                uBT_WW(i,j,0) = (1.5_rt * (duL - du0(i,j,0))) * ((FAmt_L - FA_avg) / (FAmt_L - FA_0));
            }

            FA_0 = FAmt_0; FA_avg = FAmt_0;
            if ((duR - du0(i,j,0)) != 0.0_rt) {
                FA_avg = uhtot_R / (duR - du0(i,j,0));
            }
            if (FA_avg > amrex::max(FA_0, FAmt_R)) {
                FA_avg = amrex::max(FA_0, FAmt_R);
            } else if (FA_avg < amrex::min(FA_0, FAmt_R)) {
                FA_0 = FA_avg;
            }

            FA_u_E0(i,j,0) = FA_0; FA_u_EE(i,j,0) = FAmt_R;
            if (amrex::Math::abs(FAmt_R - FA_0) <= 1.0e-12_rt*FA_0) {
                uBT_EE(i,j,0) = 0.0_rt;
            } else {
                uBT_EE(i,j,0) = (1.5_rt * (duR - du0(i,j,0))) * ((FAmt_R - FA_avg) / (FAmt_R - FA_0));
            }
        } else {
            FA_u_W0(i,j,0) = 0.0_rt; FA_u_WW(i,j,0) = 0.0_rt;
            FA_u_E0(i,j,0) = 0.0_rt; FA_u_EE(i,j,0) = 0.0_rt;
            uBT_WW(i,j,0) = 0.0_rt; uBT_EE(i,j,0) = 0.0_rt;
        }
    });
}

//> Sets the effective open face areas and barotropic velocity corrections at
//  meridional faces that reproduce the summed layer transports for three
//  test barotropic velocities (southerly, northerly, and zero-adjustment),
//  for use in the barotropic-consistency iteration.
void set_merid_BT_cont(
    const Box& bxC,                          //!< Iteration box for continuity solver
    Array4<const Real> const& v,             //!< Meridional velocity
    Array4<const Real> const& h_in,          //!< Layer thickness used to calculate fluxes
    Array4<const Real> const& h_S,           //!< South edge thickness in the reconstruction
    Array4<const Real> const& h_N,           //!< North edge thickness in the reconstruction
    Array4<Real> const& FA_v_S0,             //!< Effective open face area, south, 0 transport
    Array4<Real> const& FA_v_N0,             //!< Effective open face area, north, 0 transport
    Array4<Real> const& FA_v_SS,             //!< Effective open face area, southerly test velocity
    Array4<Real> const& FA_v_NN,             //!< Effective open face area, northerly test velocity
    Array4<Real> const& vBT_SS,              //!< Southerly correction to the barotropic velocity
    Array4<Real> const& vBT_NN,              //!< Northerly correction to the barotropic velocity
    Array4<const Real> const& dv0,           //!< Barotropic velocity increment that gives 0 transport
    Array4<const Real> const& vh_tot_0,      //!< Summed transport with 0 adjustment
    Array4<const Real> const& dvhdv_tot_0,   //!< Partial derivative of du_err with dv at 0 adjustment
    Array4<const Real> const& dv_max_CFL,    //!< Maximum acceptable value of dv
    Array4<const Real> const& dv_min_CFL,    //!< Minimum acceptable value of dv
    Real dt,                                 //!< Time increment
    Array4<const Real> const& dyCv,          //!< The grid cell's v-point y-extent
    Array4<const Real> const& dx_Cv,         //!< Unblocked v-face length (2D, addressed at k=0)
    Array4<const Real> const& IareaT,        //!< 1/areaT (2D, addressed at k=0)
    Array4<const Real> const& IdyT,          //!< 1/dyT (2D, addressed at k=0)
    const transport_adjust_CS_C& CS,         //!< Transport-adjustment and barotropic-consistency options
    Array4<const Real> const& visc_rem,      //!< Fraction of momentum/barotropic acceleration
                                              //!< remaining after viscosity
    Array4<const Real> const& visc_rem_max,  //!< Maximum allowable viscosity remnant
    Array4<const int> const& do_I,           //!< Logical flag (0/1) indicating which I values to work on
    Array4<const Real> const& por_face_areaV)//!< Fractional open area of V-faces
{
    BL_PROFILE("set_merid_BT_cont");

    const Real Idt = 1.0_rt / dt;
    const Real min_visc_rem = 0.1_rt;
    const Real CFL_min = 1.0e-6_rt;

    const int kmin = bxC.smallEnd(2);
    const int kmax = bxC.bigEnd(2);

    // Iteration box for v-point (V-grid) fields: grown by 1 at the lower y-extent
    Box bxV = growLo(bxC, 1, 1);
    Box bx2d(IntVect(bxV.smallEnd(0), bxV.smallEnd(1), 0),
             IntVect(bxV.bigEnd(0),   bxV.bigEnd(1),   0));

    ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        const bool active = (do_I(i,j,0) != 0);

        // Determine the southerly- and northerly- fluxes. Choose a sufficiently
        // negative velocity correction for the northerly-flux, and a sufficiently
        // positive correction for the southerly-flux.
        const Real dv_CFL = (CFL_min * Idt) * dyCv(i,j,0);
        Real dvR = amrex::min(0.0_rt, dv0(i,j,0) - dv_CFL);
        Real dvL = amrex::max(0.0_rt, dv0(i,j,0) + dv_CFL);
        Real FAmt_L = 0.0_rt, FAmt_R = 0.0_rt, FAmt_0 = 0.0_rt;
        Real vhtot_L = 0.0_rt, vhtot_R = 0.0_rt;

        if (active) {
            for (int k = kmin; k <= kmax; ++k) {
                Real const visc_rem_lim = amrex::max(visc_rem(i,j,k), min_visc_rem*visc_rem_max(i,j,0));
                if (visc_rem_lim > 0.0_rt) { // This is almost always true for ocean points.
                    if (v(i,j,k) + dvR*visc_rem_lim > -dv_CFL*visc_rem(i,j,k)) {
                        dvR = -(v(i,j,k) + dv_CFL*visc_rem(i,j,k)) / visc_rem_lim;
                    }
                    if (v(i,j,k) + dvL*visc_rem_lim < dv_CFL*visc_rem(i,j,k)) {
                        dvL = -(v(i,j,k) - dv_CFL*visc_rem(i,j,k)) / visc_rem_lim;
                    }
                }
            }

            for (int k = kmin; k <= kmax; ++k) {
                Real const v_L = v(i,j,k) + dvL * visc_rem(i,j,k);
                Real const v_R = v(i,j,k) + dvR * visc_rem(i,j,k);
                Real const v_0 = v(i,j,k) + dv0(i,j,0) * visc_rem(i,j,k);
                Real vh_0, vh_L, vh_R, dvhdv_0, dvhdv_L, dvhdv_R;
                flux_elem_point(v_0, h_in(i,j,k), h_in(i,j+1,k), h_S(i,j,k), h_S(i,j+1,k),
                                h_N(i,j,k), h_N(i,j+1,k), vh_0, dvhdv_0, visc_rem(i,j,k),
                                dx_Cv(i,j,0), IareaT(i,j,0), IareaT(i,j+1,0), IdyT(i,j,0), IdyT(i,j+1,0),
                                dt, CS.vol_CFL, por_face_areaV(i,j,k));
                flux_elem_point(v_L, h_in(i,j,k), h_in(i,j+1,k), h_S(i,j,k), h_S(i,j+1,k),
                                h_N(i,j,k), h_N(i,j+1,k), vh_L, dvhdv_L, visc_rem(i,j,k),
                                dx_Cv(i,j,0), IareaT(i,j,0), IareaT(i,j+1,0), IdyT(i,j,0), IdyT(i,j+1,0),
                                dt, CS.vol_CFL, por_face_areaV(i,j,k));
                flux_elem_point(v_R, h_in(i,j,k), h_in(i,j+1,k), h_S(i,j,k), h_S(i,j+1,k),
                                h_N(i,j,k), h_N(i,j+1,k), vh_R, dvhdv_R, visc_rem(i,j,k),
                                dx_Cv(i,j,0), IareaT(i,j,0), IareaT(i,j+1,0), IdyT(i,j,0), IdyT(i,j+1,0),
                                dt, CS.vol_CFL, por_face_areaV(i,j,k));
                FAmt_0 += dvhdv_0;
                FAmt_L += dvhdv_L;
                FAmt_R += dvhdv_R;
                vhtot_L += vh_L;
                vhtot_R += vh_R;
            }

            Real FA_0 = FAmt_0, FA_avg = FAmt_0;
            if ((dvL - dv0(i,j,0)) != 0.0_rt) {
                FA_avg = vhtot_L / (dvL - dv0(i,j,0));
            }
            if (FA_avg > amrex::max(FA_0, FAmt_L)) {
                FA_avg = amrex::max(FA_0, FAmt_L);
            } else if (FA_avg < amrex::min(FA_0, FAmt_L)) {
                FA_0 = FA_avg;
            }

            FA_v_S0(i,j,0) = FA_0; FA_v_SS(i,j,0) = FAmt_L;
            if (amrex::Math::abs(FA_0 - FAmt_L) <= 1.0e-12_rt*FA_0) {
                vBT_SS(i,j,0) = 0.0_rt;
            } else {
                vBT_SS(i,j,0) = (1.5_rt * (dvL - dv0(i,j,0))) * ((FAmt_L - FA_avg) / (FAmt_L - FA_0));
            }

            FA_0 = FAmt_0; FA_avg = FAmt_0;
            if ((dvR - dv0(i,j,0)) != 0.0_rt) {
                FA_avg = vhtot_R / (dvR - dv0(i,j,0));
            }
            if (FA_avg > amrex::max(FA_0, FAmt_R)) {
                FA_avg = amrex::max(FA_0, FAmt_R);
            } else if (FA_avg < amrex::min(FA_0, FAmt_R)) {
                FA_0 = FA_avg;
            }

            FA_v_N0(i,j,0) = FA_0; FA_v_NN(i,j,0) = FAmt_R;
            if (amrex::Math::abs(FAmt_R - FA_0) <= 1.0e-12_rt*FA_0) {
                vBT_NN(i,j,0) = 0.0_rt;
            } else {
                vBT_NN(i,j,0) = (1.5_rt * (dvR - dv0(i,j,0))) * ((FAmt_R - FA_avg) / (FAmt_R - FA_0));
            }
        } else {
            FA_v_S0(i,j,0) = 0.0_rt; FA_v_SS(i,j,0) = 0.0_rt;
            FA_v_N0(i,j,0) = 0.0_rt; FA_v_NN(i,j,0) = 0.0_rt;
            vBT_SS(i,j,0) = 0.0_rt; vBT_NN(i,j,0) = 0.0_rt;
        }
    });
}
}
