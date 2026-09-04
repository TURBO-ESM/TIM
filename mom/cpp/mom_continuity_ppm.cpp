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

//> Newton-iterates a barotropic velocity correction per zonal face so that
//  the vertically-summed zonal mass/volume transport matches the target
//  barotropic transport, to within the transport-adjustment iteration's
//  tolerance. Always completes the fixed-count itt-loop rather than exiting
//  early once every column in a row has converged, matching the Fortran
//  source's own OpenMP-target-compiled path -- the alternative (a
//  data-dependent per-row early exit) is disabled there because it
//  serializes on GPU-style parallel execution; do_I masks further updates
//  to a column once it has converged.
void zonal_flux_adjust(
    const Box& bxC,                          //!< Iteration box for continuity solver
    Array4<const Real> const& u,             //!< Zonal velocity [L T-1 ~> m s-1]
    Array4<const Real> const& h_in,          //!< Layer thickness used to calculate fluxes [H ~> m or kg m-2]
    Array4<const Real> const& h_W,           //!< West edge thickness in the reconstruction [H ~> m or kg m-2]
    Array4<const Real> const& h_E,           //!< East edge thickness in the reconstruction [H ~> m or kg m-2]
    Array4<const Real> const& uh_tot_0,      //!< Summed transport with 0 adjustment [H L2 T-1 ~> m3 s-1 or kg s-1]
    Array4<const Real> const& duhdu_tot_0,   //!< Partial derivative of du_err with du at 0 adjustment
                                              //!< [H L ~> m2 or kg m-1]
    Array4<Real> const& du,                  //!< The barotropic velocity adjustment [L T-1 ~> m s-1]
    Array4<const Real> const& du_max_CFL,    //!< Maximum acceptable value of du [L T-1 ~> m s-1]
    Array4<const Real> const& du_min_CFL,    //!< Minimum acceptable value of du [L T-1 ~> m s-1]
    Real dt,                                 //!< Time increment [T ~> s]
    Array4<const Real> const& dy_Cu,         //!< The grid cell's unblocked lengths of the u-faces
                                              //!< of the h-cell [L ~> m]
    Array4<const Real> const& IareaT,        //!< The grid cell's 1/areaT [L-2 ~> m-2]
    Array4<const Real> const& IdxT,          //!< The grid cell's 1/dxT [L-1 ~> m-1]
    const transport_adjust_CS_C& CS,         //!< Options controlling the transport adjustment
                                              //!< and barotropic-consistency iteration
    Array4<const Real> const& visc_rem,      //!< Fraction of momentum/barotropic acceleration remaining
                                              //!< after viscosity [nondim]
    Array4<const int> const& do_I_in,        //!< Logical flag (0/1) indicating which I values to work on
    Array4<const Real> const& por_face_areaU,//!< Fractional open area of U-faces [nondim]
    Array4<const Real> const& uhbt,          //!< Summed volume flux through zonal faces
                                              //!< [H L2 T-1 ~> m3 s-1 or kg s-1]; may be absent (.p == nullptr)
    Array4<Real> const& uh_3d,               //!< Volume flux through zonal faces, u*h*dy
                                              //!< [H L2 T-1 ~> m3 s-1 or kg s-1]; may be absent (.p == nullptr)
    OceanOBC* obc)                           //!< Open boundary control structure
{
    BL_PROFILE("zonal_flux_adjust");

    // NOTE: OBC support temporarily disabled.
    // OceanOBC is forward-declared only.
    // All boundary-condition logic removed for initial port validation.
    if (obc != nullptr) {
       AMREX_ABORT_LOC("OBC pointer provided but not yet implemented");
    }
    /*
    bool local_OBC = false;
    if (obc != nullptr) local_OBC = obc->open_u_BCs_exist_globally;
    */

    const bool use_uhbt  = (uhbt.p != nullptr);
    const bool use_uh_3d = (uh_3d.p != nullptr);

    const Real tol_vel   = CS.tol_vel;
    const int  max_itts  = 20;

    const int kmin = bxC.smallEnd(2);
    const int kmax = bxC.bigEnd(2);

    // Iteration box for u-point (U-grid) fields: grown by 1 at the lower x-extent
    Box bxU = growLo(bxC, 0, 1);
    Box bx2d(IntVect(bxU.smallEnd(0), bxU.smallEnd(1), 0),
             IntVect(bxU.bigEnd(0),   bxU.bigEnd(1),   0));

    ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        bool do_I = (do_I_in(i,j,0) != 0);
        Real du_val = 0.0_rt;
        Real du_max = du_max_CFL(i,j,0);
        Real du_min = du_min_CFL(i,j,0);
        Real uh_err = uh_tot_0(i,j,0);
        if (use_uhbt) uh_err -= uhbt(i,j,0);
        Real duhdu_tot = duhdu_tot_0(i,j,0);
        Real uh_err_best = amrex::Math::abs(uh_err);

        for (int itt = 1; itt <= max_itts; ++itt) {
            Real tol_eta;
            if (itt <= 1) {
                tol_eta = 1.0e-6_rt * CS.tol_eta;
            } else if (itt == 2) {
                tol_eta = 1.0e-4_rt * CS.tol_eta;
            } else if (itt == 3) {
                tol_eta = 1.0e-2_rt * CS.tol_eta;
            } else {
                tol_eta = CS.tol_eta;
            }

            if (do_I) {
                if (uh_err > 0.0_rt) {
                    du_max = du_val;
                } else if (uh_err < 0.0_rt) {
                    du_min = du_val;
                } else {
                    do_I = false;
                }

                if ((dt * amrex::min(IareaT(i,j,0), IareaT(i+1,j,0)) * amrex::Math::abs(uh_err) > tol_eta) ||
                    (CS.better_iter &&
                     ((amrex::Math::abs(uh_err) > tol_vel * duhdu_tot) ||
                      (amrex::Math::abs(uh_err) > uh_err_best)))) {
                    // Use Newton's method, provided it stays bounded. Otherwise bisect
                    // the value with the appropriate bound.
                    Real const ddu = -uh_err / duhdu_tot;
                    Real const du_prev = du_val;
                    du_val = du_val + ddu;
                    if (amrex::Math::abs(ddu) < 1.0e-15_rt * amrex::Math::abs(du_val)) {
                        do_I = false; // ddu is small enough to quit.
                    } else if (ddu > 0.0_rt) {
                        if (du_val >= du_max) {
                            du_val = 0.5_rt * (du_prev + du_max);
                            if (du_max - du_prev < 1.0e-15_rt * amrex::Math::abs(du_val)) do_I = false;
                        }
                    } else { // ddu < 0.0
                        if (du_val <= du_min) {
                            du_val = 0.5_rt * (du_prev + du_min);
                            if (du_prev - du_min < 1.0e-15_rt * amrex::Math::abs(du_val)) do_I = false;
                        }
                    }
                } else {
                    do_I = false;
                }
            }

            if ((itt < max_itts) || use_uh_3d) {
                uh_err = 0.0_rt;
                duhdu_tot = 0.0_rt;
                if (use_uhbt) uh_err = -uhbt(i,j,0);
                if (do_I) {
                    for (int k = kmin; k <= kmax; ++k) {
                        Real const u_new = u(i,j,k) + du_val * visc_rem(i,j,k);
                        Real uh_val, duhdu_val;
                        flux_elem_point(u_new, h_in(i,j,k), h_in(i+1,j,k), h_W(i,j,k), h_W(i+1,j,k),
                                        h_E(i,j,k), h_E(i+1,j,k), uh_val, duhdu_val, visc_rem(i,j,k),
                                        dy_Cu(i,j,0), IareaT(i,j,0), IareaT(i+1,j,0), IdxT(i,j,0), IdxT(i+1,j,0),
                                        dt, CS.vol_CFL, por_face_areaU(i,j,k));
                        /*
                        if (local_OBC) {
                            int const l_seg = obc->segnum_u(i,j,k);
                            if (l_seg != 0 && obc->segment[amrex::Math::abs(l_seg)].open) {
                                if (l_seg > 0) {
                                    uh_val    = (dy_Cu(i,j,0) * por_face_areaU(i,j,k)) * u_new * h_in(i,j,k);
                                    duhdu_val = (dy_Cu(i,j,0) * por_face_areaU(i,j,k)) * h_in(i,j,k) * visc_rem(i,j,k);
                                } else {
                                    uh_val    = (dy_Cu(i,j,0) * por_face_areaU(i,j,k)) * u_new * h_in(i+1,j,k);
                                    duhdu_val = (dy_Cu(i,j,0) * por_face_areaU(i,j,k)) * h_in(i+1,j,k) * visc_rem(i,j,k);
                                }
                            }
                        }
                        */
                        if (use_uh_3d) uh_3d(i,j,k) = uh_val;
                        uh_err += uh_val;
                        duhdu_tot += duhdu_val;
                    }
                }
                uh_err_best = amrex::min(uh_err_best, amrex::Math::abs(uh_err));
            }
        }

        du(i,j,0) = du_val;
    });
}

//> Newton-iterates a barotropic velocity correction per meridional face so
//  that the vertically-summed meridional mass/volume transport matches the
//  target barotropic transport, to within the transport-adjustment
//  iteration's tolerance. Always completes the fixed-count itt-loop rather
//  than exiting early once every column in a row has converged, matching
//  the Fortran source's own OpenMP-target-compiled path -- the alternative
//  (a data-dependent per-row early exit) is disabled there because it
//  serializes on GPU-style parallel execution; do_I masks further updates
//  to a column once it has converged.
void meridional_flux_adjust(
    const Box& bxC,                          //!< Iteration box for continuity solver
    Array4<const Real> const& v,             //!< Meridional velocity [L T-1 ~> m s-1]
    Array4<const Real> const& h_in,          //!< Layer thickness used to calculate fluxes [H ~> m or kg m-2]
    Array4<const Real> const& h_S,           //!< South edge thickness in the reconstruction [H ~> m or kg m-2]
    Array4<const Real> const& h_N,           //!< North edge thickness in the reconstruction [H ~> m or kg m-2]
    Array4<const Real> const& vh_tot_0,      //!< Summed transport with 0 adjustment [H L2 T-1 ~> m3 s-1 or kg s-1]
    Array4<const Real> const& dvhdv_tot_0,   //!< Partial derivative of dv_err with dv at 0 adjustment
                                              //!< [H L ~> m2 or kg m-1]
    Array4<Real> const& dv,                  //!< The barotropic velocity adjustment [L T-1 ~> m s-1]
    Array4<const Real> const& dv_max_CFL,    //!< Maximum acceptable value of dv [L T-1 ~> m s-1]
    Array4<const Real> const& dv_min_CFL,    //!< Minimum acceptable value of dv [L T-1 ~> m s-1]
    Real dt,                                 //!< Time increment [T ~> s]
    Array4<const Real> const& dx_Cv,         //!< The grid cell's unblocked lengths of the v-faces
                                              //!< of the h-cell [L ~> m]
    Array4<const Real> const& IareaT,        //!< The grid cell's 1/areaT [L-2 ~> m-2]
    Array4<const Real> const& IdyT,          //!< The grid cell's 1/dyT [L-1 ~> m-1]
    const transport_adjust_CS_C& CS,         //!< Options controlling the transport adjustment
                                              //!< and barotropic-consistency iteration
    Array4<const Real> const& visc_rem,      //!< Fraction of momentum/barotropic acceleration remaining
                                              //!< after viscosity [nondim]
    Array4<const int> const& do_I_in,        //!< Logical flag (0/1) indicating which I values to work on
    Array4<const Real> const& por_face_areaV,//!< Fractional open area of V-faces [nondim]
    Array4<const Real> const& vhbt,          //!< Summed volume flux through meridional faces
                                              //!< [H L2 T-1 ~> m3 s-1 or kg s-1]; may be absent (.p == nullptr)
    Array4<Real> const& vh_3d,               //!< Volume flux through meridional faces, v*h*dx
                                              //!< [H L2 T-1 ~> m3 s-1 or kg s-1]; may be absent (.p == nullptr)
    OceanOBC* obc)                           //!< Open boundary control structure
{
    BL_PROFILE("meridional_flux_adjust");

    // NOTE: OBC support temporarily disabled.
    // OceanOBC is forward-declared only.
    // All boundary-condition logic removed for initial port validation.
    if (obc != nullptr) {
       AMREX_ABORT_LOC("OBC pointer provided but not yet implemented");
    }
    /*
    bool local_OBC = false;
    if (obc != nullptr) local_OBC = obc->open_u_BCs_exist_globally;
    */

    const bool use_vhbt  = (vhbt.p != nullptr);
    const bool use_vh_3d = (vh_3d.p != nullptr);

    const Real tol_vel   = CS.tol_vel;
    const int  max_itts  = 20;

    const int kmin = bxC.smallEnd(2);
    const int kmax = bxC.bigEnd(2);

    // Iteration box for v-point (V-grid) fields: grown by 1 at the lower y-extent
    Box bxV = growLo(bxC, 1, 1);
    Box bx2d(IntVect(bxV.smallEnd(0), bxV.smallEnd(1), 0),
             IntVect(bxV.bigEnd(0),   bxV.bigEnd(1),   0));

    ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        bool do_I = (do_I_in(i,j,0) != 0);
        Real dv_val = 0.0_rt;
        Real dv_max = dv_max_CFL(i,j,0);
        Real dv_min = dv_min_CFL(i,j,0);
        Real vh_err = vh_tot_0(i,j,0);
        if (use_vhbt) vh_err -= vhbt(i,j,0);
        Real dvhdv_tot = dvhdv_tot_0(i,j,0);
        Real vh_err_best = amrex::Math::abs(vh_err);

        for (int itt = 1; itt <= max_itts; ++itt) {
            Real tol_eta;
            if (itt <= 1) {
                tol_eta = 1.0e-6_rt * CS.tol_eta;
            } else if (itt == 2) {
                tol_eta = 1.0e-4_rt * CS.tol_eta;
            } else if (itt == 3) {
                tol_eta = 1.0e-2_rt * CS.tol_eta;
            } else {
                tol_eta = CS.tol_eta;
            }

            // Unmasked: runs every iteration regardless of do_I.
            if (vh_err > 0.0_rt) {
                dv_max = dv_val;
            } else if (vh_err < 0.0_rt) {
                dv_min = dv_val;
            } else {
                do_I = false;
            }

            // do_I-masked: Newton step / bisection.
            if (do_I) {
                if ((dt * amrex::min(IareaT(i,j,0), IareaT(i,j+1,0)) * amrex::Math::abs(vh_err) > tol_eta) ||
                    (CS.better_iter &&
                     ((amrex::Math::abs(vh_err) > tol_vel * dvhdv_tot) ||
                      (amrex::Math::abs(vh_err) > vh_err_best)))) {
                    // Use Newton's method, provided it stays bounded. Otherwise bisect
                    // the value with the appropriate bound.
                    Real const ddv = -vh_err / dvhdv_tot;
                    Real const dv_prev = dv_val;
                    dv_val = dv_val + ddv;
                    if (amrex::Math::abs(ddv) < 1.0e-15_rt * amrex::Math::abs(dv_val)) {
                        do_I = false; // ddv is small enough to quit.
                    } else if (ddv > 0.0_rt) {
                        if (dv_val >= dv_max) {
                            dv_val = 0.5_rt * (dv_prev + dv_max);
                            if (dv_max - dv_prev < 1.0e-15_rt * amrex::Math::abs(dv_val)) do_I = false;
                        }
                    } else { // ddv < 0.0
                        if (dv_val <= dv_min) {
                            dv_val = 0.5_rt * (dv_prev + dv_min);
                            if (dv_prev - dv_min < 1.0e-15_rt * amrex::Math::abs(dv_val)) do_I = false;
                        }
                    }
                } else {
                    do_I = false;
                }
            }

            if ((itt < max_itts) || use_vh_3d) {
                vh_err = 0.0_rt;
                dvhdv_tot = 0.0_rt;
                if (use_vhbt) vh_err = -vhbt(i,j,0);
                if (do_I) {
                    for (int k = kmin; k <= kmax; ++k) {
                        Real const v_new = v(i,j,k) + dv_val * visc_rem(i,j,k);
                        Real vh_val, dvhdv_val;
                        flux_elem_point(v_new, h_in(i,j,k), h_in(i,j+1,k), h_S(i,j,k), h_S(i,j+1,k),
                                        h_N(i,j,k), h_N(i,j+1,k), vh_val, dvhdv_val, visc_rem(i,j,k),
                                        dx_Cv(i,j,0), IareaT(i,j,0), IareaT(i,j+1,0), IdyT(i,j,0), IdyT(i,j+1,0),
                                        dt, CS.vol_CFL, por_face_areaV(i,j,k));
                        /*
                        if (local_OBC) {
                            int const l_seg = obc->segnum_v(i,j,k);
                            if (l_seg != 0 && obc->segment[amrex::Math::abs(l_seg)].open) {
                                if (l_seg > 0) {
                                    vh_val    = (dx_Cv(i,j,0) * por_face_areaV(i,j,k)) * v_new * h_in(i,j,k);
                                    dvhdv_val = (dx_Cv(i,j,0) * por_face_areaV(i,j,k)) * h_in(i,j,k) * visc_rem(i,j,k);
                                } else {
                                    vh_val    = (dx_Cv(i,j,0) * por_face_areaV(i,j,k)) * v_new * h_in(i,j+1,k);
                                    dvhdv_val = (dx_Cv(i,j,0) * por_face_areaV(i,j,k)) * h_in(i,j+1,k) * visc_rem(i,j,k);
                                }
                            }
                        }
                        */
                        if (use_vh_3d) vh_3d(i,j,k) = vh_val;
                        vh_err += vh_val;
                        dvhdv_tot += dvhdv_val;
                    }
                    vh_err_best = amrex::min(vh_err_best, amrex::Math::abs(vh_err));
                }
            }
        }

        dv(i,j,0) = dv_val;
    });
}

//> Accumulates the vertically-summed zonal barotropic mass/volume transport
//  across the water column, for use as the barotropic solver's target
//  transport in the transport-adjustment iteration.
void zonal_BT_mass_flux(
    const Box& bxC,                          //!< Iteration box for continuity solver
    Array4<const Real> const& u,             //!< Zonal velocity [L T-1 ~> m s-1]
    Array4<const Real> const& h_in,          //!< Layer thickness used to calculate fluxes [H ~> m or kg m-2]
    Array4<const Real> const& h_W,           //!< Western edge thickness in the PPM reconstruction
                                              //!< [H ~> m or kg m-2]
    Array4<const Real> const& h_E,           //!< Eastern edge thickness in the PPM reconstruction
                                              //!< [H ~> m or kg m-2]
    Array4<Real> const& uhbt,                //!< Summed volume flux through zonal faces
                                              //!< [H L2 T-1 ~> m3 s-1 or kg s-1]
    Real dt,                                 //!< Time increment [T ~> s]
    Array4<const Real> const& dy_Cu,         //!< The grid cell's unblocked lengths of the u-faces
                                              //!< of the h-cell [L ~> m]
    Array4<const Real> const& IareaT,        //!< The grid cell's 1/areaT [L-2 ~> m-2]
    Array4<const Real> const& IdxT,          //!< The grid cell's 1/dxT [L-1 ~> m-1]
    const transport_adjust_CS_C& CS,         //!< Options controlling the transport adjustment
                                              //!< and barotropic-consistency iteration
    OceanOBC* obc,                           //!< Open boundary control structure
    Array4<const Real> const& por_face_areaU)//!< Fractional open area of U-faces [nondim]
{
    BL_PROFILE("zonal_BT_mass_flux");

    // NOTE: OBC support temporarily disabled.
    // OceanOBC is forward-declared only.
    // All boundary-condition logic removed for initial port validation.
    if (obc != nullptr) {
       AMREX_ABORT_LOC("OBC pointer provided but not yet implemented");
    }
    /*
    bool local_specified_BC = false;
    if (obc != nullptr) {
        if (obc->OBC_pe) {
            local_specified_BC = obc->specified_v_BCs_exist_globally;
        }
    }
    */

    const int kmin = bxC.smallEnd(2);
    const int kmax = bxC.bigEnd(2);

    // Iteration box for u-point (U-grid) fields: grown by 1 at the lower x-extent
    Box bxU = growLo(bxC, 0, 1);
    Box bx2d(IntVect(bxU.smallEnd(0), bxU.smallEnd(1), 0),
             IntVect(bxU.bigEnd(0),   bxU.bigEnd(1),   0));

    ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        Real uhbt_val = 0.0_rt;
        for (int k = kmin; k <= kmax; ++k) {
            Real uh_val, duhdu_val;
            flux_elem_point(u(i,j,k), h_in(i,j,k), h_in(i+1,j,k), h_W(i,j,k), h_W(i+1,j,k),
                            h_E(i,j,k), h_E(i+1,j,k), uh_val, duhdu_val, 1.0_rt,
                            dy_Cu(i,j,0), IareaT(i,j,0), IareaT(i+1,j,0), IdxT(i,j,0), IdxT(i+1,j,0),
                            dt, CS.vol_CFL, por_face_areaU(i,j,k));
            /*
            // untested (Fortran source's own comment)!
            if (local_specified_BC) {
                flux_elem_obc_point(u(i,j,k), h_in(i,j,k), h_in(i+1,j,k), uh_val, duhdu_val, 1.0_rt,
                                    por_face_areaU(i,j,k), dy_Cu(i,j,0), obc, obc->segnum_u(i,j,0));
            }
            // Second pass in the Fortran source: if a specified-transport OBC
            // segment is active in this row, uh(I,j,k) is replaced by the
            // segment's normal transport before being summed into uhbt.
            if (local_specified_BC && obc->segnum_u(i,j,0) != 0) {
                int const l_seg = amrex::Math::abs(obc->segnum_u(i,j,0));
                if (obc->segment[l_seg].specified) uh_val = obc->segment[l_seg].normal_trans(i,j,k);
            }
            */
            uhbt_val += uh_val;
        }
        uhbt(i,j,0) = uhbt_val;
    });
}

//> Accumulates the vertically-summed meridional barotropic mass/volume
//  transport across the water column, for use as the barotropic solver's
//  target transport in the transport-adjustment iteration.
void meridional_BT_mass_flux(
    const Box& bxC,                          //!< Iteration box for continuity solver
    Array4<const Real> const& v,             //!< Meridional velocity [L T-1 ~> m s-1]
    Array4<const Real> const& h_in,          //!< Layer thickness used to calculate fluxes [H ~> m or kg m-2]
    Array4<const Real> const& h_S,           //!< Southern edge thickness in the PPM reconstruction
                                              //!< [H ~> m or kg m-2]
    Array4<const Real> const& h_N,           //!< Northern edge thickness in the PPM reconstruction
                                              //!< [H ~> m or kg m-2]
    Array4<Real> const& vhbt,                //!< Summed volume flux through meridional faces
                                              //!< [H L2 T-1 ~> m3 s-1 or kg s-1]
    Real dt,                                 //!< Time increment [T ~> s]
    Array4<const Real> const& dx_Cv,         //!< The grid cell's unblocked lengths of the v-faces
                                              //!< of the h-cell [L ~> m]
    Array4<const Real> const& IareaT,        //!< The grid cell's 1/areaT [L-2 ~> m-2]
    Array4<const Real> const& IdyT,          //!< The grid cell's 1/dyT [L-1 ~> m-1]
    const transport_adjust_CS_C& CS,         //!< Options controlling the transport adjustment
                                              //!< and barotropic-consistency iteration
    OceanOBC* obc,                           //!< Open boundary control structure
    Array4<const Real> const& por_face_areaV)//!< Fractional open area of V-faces [nondim]
{
    BL_PROFILE("meridional_BT_mass_flux");

    // NOTE: OBC support temporarily disabled.
    // OceanOBC is forward-declared only.
    // All boundary-condition logic removed for initial port validation.
    if (obc != nullptr) {
       AMREX_ABORT_LOC("OBC pointer provided but not yet implemented");
    }
    /*
    bool local_specified_BC = false;
    if (obc != nullptr) {
        if (obc->OBC_pe) {
            local_specified_BC = obc->specified_v_BCs_exist_globally;
        }
    }
    */

    const int kmin = bxC.smallEnd(2);
    const int kmax = bxC.bigEnd(2);

    // Iteration box for v-point (V-grid) fields: grown by 1 at the lower y-extent
    Box bxV = growLo(bxC, 1, 1);
    Box bx2d(IntVect(bxV.smallEnd(0), bxV.smallEnd(1), 0),
             IntVect(bxV.bigEnd(0),   bxV.bigEnd(1),   0));

    ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        Real vhbt_val = 0.0_rt;
        for (int k = kmin; k <= kmax; ++k) {
            Real vh_val, dvhdv_val;
            flux_elem_point(v(i,j,k), h_in(i,j,k), h_in(i,j+1,k), h_S(i,j,k), h_S(i,j+1,k),
                            h_N(i,j,k), h_N(i,j+1,k), vh_val, dvhdv_val, 1.0_rt,
                            dx_Cv(i,j,0), IareaT(i,j,0), IareaT(i,j+1,0), IdyT(i,j,0), IdyT(i,j+1,0),
                            dt, CS.vol_CFL, por_face_areaV(i,j,k));
            /*
            // untested (Fortran source's own comment)!
            if (local_specified_BC) {
                flux_elem_obc_point(v(i,j,k), h_in(i,j,k), h_in(i,j+1,k), vh_val, dvhdv_val, 1.0_rt,
                                    por_face_areaV(i,j,k), dx_Cv(i,j,0), obc, obc->segnum_v(i,j,0));
            }
            // Second pass in the Fortran source: if a specified-transport OBC
            // segment is active in this row, vh(i,J,k) is replaced by the
            // segment's normal transport before being summed into vhbt.
            if (local_specified_BC && obc->segnum_v(i,j,0) != 0) {
                int const l_seg = amrex::Math::abs(obc->segnum_v(i,j,0));
                if (obc->segment[l_seg].specified) vh_val = obc->segment[l_seg].normal_trans(i,j,k);
            }
            */
            vhbt_val += vh_val;
        }
        vhbt(i,j,0) = vhbt_val;
    });
}

//> Calculates the mass or volume fluxes through the zonal faces, and other
//  related quantities -- including, optionally, the barotropic mass-flux
//  correction (u_cor/du_cor) and the barotropic-consistency
//  face-area/velocity-correction diagnostics (FA_u_*/uBT_*).
void zonal_mass_flux(
    const Box& bxC,                          //!< Iteration box for continuity solver
    Array4<const Real> const& u,             //!< Zonal velocity [L T-1 ~> m s-1]
    Array4<const Real> const& h_in,          //!< Layer thickness used to calculate fluxes [H ~> m or kg m-2]
    Array4<const Real> const& h_W,           //!< West edge thickness in the reconstruction [H ~> m or kg m-2]
    Array4<const Real> const& h_E,           //!< East edge thickness in the reconstruction [H ~> m or kg m-2]
    Array4<Real> const& uh,                  //!< Volume flux through zonal faces, u*h*dy
                                              //!< [H L2 T-1 ~> m3 s-1 or kg s-1]
    Real dt,                                 //!< Time increment [T ~> s]
    Array4<const Real> const& dy_Cu,         //!< The grid cell's unblocked lengths of the u-faces
                                              //!< of the h-cell [L ~> m]
    Array4<const Real> const& IareaT,        //!< The grid cell's 1/areaT [L-2 ~> m-2]
    Array4<const Real> const& IdxT,          //!< The grid cell's 1/dxT [L-1 ~> m-1]
    Array4<const Real> const& areaT,         //!< The area of the h-cell [L2 ~> m2]
    Array4<const Real> const& dxT,           //!< The x-extent of the h-cell [L ~> m]
    Array4<const Real> const& mask2dCu,      //!< 0 for land points, 1 for ocean points at u-locations [nondim]
    Array4<const Real> const& dxCu,          //!< The grid cell's u-point x-extent [L ~> m]
    Real H_subroundoff,                      //!< A negligibly small thickness used to avoid division
                                              //!< by zero [H ~> m or kg m-2]
    const transport_adjust_CS_C& CS,         //!< Options controlling the transport adjustment
                                              //!< and barotropic-consistency iteration
    OceanOBC* obc,                           //!< Open boundary control structure
    Array4<const Real> const& por_face_areaU,//!< Fractional open area of U-faces [nondim]
    Array4<const Real> const& uhbt,          //!< Summed volume flux through zonal faces
                                              //!< [H L2 T-1 ~> m3 s-1 or kg s-1]; may be absent (.p == nullptr)
    Array4<const Real> const& visc_rem_u,    //!< Fraction of momentum/barotropic acceleration remaining
                                              //!< after viscosity [nondim]; may be absent (.p == nullptr)
    Array4<Real> const& u_cor,               //!< Zonal velocity with barotropic correction [L T-1 ~> m s-1];
                                              //!< may be absent (.p == nullptr)
    Array4<Real> const& FA_u_W0,             //!< Effective open face area, west, 0 transport;
                                              //!< may be absent (.p == nullptr)
    Array4<Real> const& FA_u_E0,             //!< Effective open face area, east, 0 transport;
                                              //!< may be absent (.p == nullptr)
    Array4<Real> const& FA_u_WW,             //!< Effective open face area, westerly test velocity;
                                              //!< may be absent (.p == nullptr)
    Array4<Real> const& FA_u_EE,             //!< Effective open face area, easterly test velocity;
                                              //!< may be absent (.p == nullptr)
    Array4<Real> const& uBT_WW,              //!< Westerly correction to the barotropic velocity;
                                              //!< may be absent (.p == nullptr)
    Array4<Real> const& uBT_EE,              //!< Easterly correction to the barotropic velocity;
                                              //!< may be absent (.p == nullptr)
    Array4<Real> const& du_cor)              //!< Zonal velocity increment from u that gives uhbt as the
                                              //!< depth-integrated transport [L T-1 ~> m s-1];
                                              //!< may be absent (.p == nullptr)
{
    BL_PROFILE("zonal_mass_flux");

    // NOTE: OBC support temporarily disabled.
    // OceanOBC is forward-declared only.
    // All boundary-condition logic removed for initial port validation.
    if (obc != nullptr) {
       AMREX_ABORT_LOC("OBC pointer provided but not yet implemented");
    }
    /*
    bool local_specified_BC = false, local_Flather_OBC = false, local_open_BC = false;
    if (obc != nullptr) {
        if (obc->OBC_pe) {
            local_specified_BC = obc->specified_u_BCs_exist_globally;
            local_Flather_OBC  = obc->Flather_u_BCs_exist_globally;
            local_open_BC      = obc->open_u_BCs_exist_globally;
        }
    }
    */

    const bool use_visc_rem = (visc_rem_u.p != nullptr);
    const bool set_BT_cont  = (FA_u_W0.p != nullptr);
    const bool need_bt      = (uhbt.p != nullptr) || set_BT_cont;
    const bool has_u_cor    = (u_cor.p != nullptr);
    const bool has_du_cor   = (du_cor.p != nullptr);

    const int kmin = bxC.smallEnd(2);
    const int kmax = bxC.bigEnd(2);

    Real CFL_dt = CS.CFL_limit_adjust / dt;
    const Real I_dt = 1.0_rt / dt;
    if (CS.aggress_adjust) CFL_dt = I_dt;

    // Iteration box for u-point (U-grid) fields: grown by 1 at the lower x-extent
    Box bxU = growLo(bxC, 0, 1);
    Box bx2d(IntVect(bxU.smallEnd(0), bxU.smallEnd(1), 0),
             IntVect(bxU.bigEnd(0),   bxU.bigEnd(1),   0));

    // Scratch, internal to this orchestrator -- never crosses the Fortran/C
    // boundary. Mirrors the Fortran source's own local allocView()/free()
    // calls (see lessons.md's backlog note on an alternative, scratch-free
    // design deferred for later).
    FArrayBox visc_rem_u_tmp_fab(bxU, 1, amrex::The_Arena());
    FArrayBox uh_tot_0_fab(bx2d, 1, amrex::The_Arena());
    FArrayBox duhdu_tot_0_fab(bx2d, 1, amrex::The_Arena());
    FArrayBox du_max_CFL_fab(bx2d, 1, amrex::The_Arena());
    FArrayBox du_min_CFL_fab(bx2d, 1, amrex::The_Arena());
    FArrayBox visc_rem_max_fab(bx2d, 1, amrex::The_Arena());
    FArrayBox du_fab(bx2d, 1, amrex::The_Arena());
    IArrayBox do_I_fab(bx2d, 1, amrex::The_Arena());

    Array4<Real> const visc_rem_u_tmp = visc_rem_u_tmp_fab.array();
    Array4<Real> const uh_tot_0       = uh_tot_0_fab.array();
    Array4<Real> const duhdu_tot_0    = duhdu_tot_0_fab.array();
    Array4<Real> const du_max_CFL     = du_max_CFL_fab.array();
    Array4<Real> const du_min_CFL     = du_min_CFL_fab.array();
    Array4<Real> const visc_rem_max   = visc_rem_max_fab.array();
    Array4<Real> const du             = du_fab.array();
    Array4<int>  const do_I           = do_I_fab.array();

    ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        if (has_du_cor) du_cor(i,j,0) = 0.0_rt;

        // Set uh (and, only when needed below, the running sums of its
        // velocity derivative and of visc_rem_u_tmp's running max).
        Real uh_tot_0_val = 0.0_rt, duhdu_tot_0_val = 0.0_rt;
        Real visc_rem_max_running = 0.0_rt;
        for (int k = kmin; k <= kmax; ++k) {
            Real const vrt = use_visc_rem ? visc_rem_u(i,j,k) : 1.0_rt;
            visc_rem_u_tmp(i,j,k) = vrt;

            Real duhdu_val;
            flux_elem_point(u(i,j,k), h_in(i,j,k), h_in(i+1,j,k), h_W(i,j,k), h_W(i+1,j,k),
                            h_E(i,j,k), h_E(i+1,j,k), uh(i,j,k), duhdu_val, vrt,
                            dy_Cu(i,j,0), IareaT(i,j,0), IareaT(i+1,j,0), IdxT(i,j,0), IdxT(i+1,j,0),
                            dt, CS.vol_CFL, por_face_areaU(i,j,k));
            /*
            if (local_open_BC) {
                int const l_seg = obc->segnum_u(i,j,0);
                flux_elem_obc_point(u(i,j,k), h_in(i,j,k), h_in(i+1,j,k), uh(i,j,k), duhdu_val,
                                    vrt, por_face_areaU(i,j,k), dy_Cu(i,j,0), obc, l_seg);
            }
            // untested (Fortran source's own comment)!
            if (local_specified_BC && obc->segnum_u(i,j,0) != 0) {
                int const l_seg = amrex::Math::abs(obc->segnum_u(i,j,0));
                if (obc->segment[l_seg].specified) uh(i,j,k) = obc->segment[l_seg].normal_trans(i,j,k);
            }
            */

            if (need_bt) {
                uh_tot_0_val += uh(i,j,k);
                duhdu_tot_0_val += duhdu_val;
                if (use_visc_rem && CS.use_visc_rem_max) {
                    visc_rem_max_running = (k == kmin) ? vrt : amrex::max(visc_rem_max_running, vrt);
                }
            }
        }

        if (need_bt) {
            Real const visc_rem_max_val = (use_visc_rem && CS.use_visc_rem_max) ? visc_rem_max_running : 1.0_rt;
            Real const I_vrm = (visc_rem_max_val > 0.0_rt) ? 1.0_rt / visc_rem_max_val : 0.0_rt;

            //   Set limits on du that will keep the CFL number between -1 and 1.
            // This should be adequate to keep the root bracketed in all cases.
            Real dx_W, dx_E;
            if (CS.vol_CFL) {
                dx_W = ratio_max_point(areaT(i,j,0), dy_Cu(i,j,0), 1000.0_rt*dxT(i,j,0));
                dx_E = ratio_max_point(areaT(i+1,j,0), dy_Cu(i,j,0), 1000.0_rt*dxT(i+1,j,0));
            } else {
                dx_W = dxT(i,j,0);
                dx_E = dxT(i+1,j,0);
            }
            Real du_max_CFL_val = 2.0_rt * (CFL_dt * dx_W) * I_vrm;
            Real du_min_CFL_val = -2.0_rt * (CFL_dt * dx_E) * I_vrm;

            if (use_visc_rem) {
                if (CS.aggress_adjust) {
                    // untested (Fortran source's own comment)!
                    for (int k = kmin; k <= kmax; ++k) {
                        Real const vrt = visc_rem_u_tmp(i,j,k);
                        Real const du_lim_max = 0.499_rt * ((dx_W*I_dt - u(i,j,k)) + amrex::min(0.0_rt, u(i-1,j,k)));
                        if (du_max_CFL_val * vrt > du_lim_max) du_max_CFL_val = du_lim_max / vrt;
                        Real const du_lim_min = 0.499_rt * ((-dx_E*I_dt - u(i,j,k)) + amrex::max(0.0_rt, u(i+1,j,k)));
                        if (du_min_CFL_val * vrt < du_lim_min) du_min_CFL_val = du_lim_min / vrt;
                    }
                } else {
                    for (int k = kmin; k <= kmax; ++k) {
                        Real const vrt = visc_rem_u_tmp(i,j,k);
                        if (du_max_CFL_val * vrt > dx_W*CFL_dt - u(i,j,k)*mask2dCu(i,j,0))
                            du_max_CFL_val = (dx_W*CFL_dt - u(i,j,k)) / vrt;
                        if (du_min_CFL_val * vrt < -dx_E*CFL_dt - u(i,j,k)*mask2dCu(i,j,0))
                            du_min_CFL_val = -(dx_E*CFL_dt + u(i,j,k)) / vrt;
                    }
                }
            } else {
                // untested (Fortran source's own comment)!
                if (CS.aggress_adjust) {
                    for (int k = kmin; k <= kmax; ++k) {
                        du_max_CFL_val = amrex::min(du_max_CFL_val,
                            0.499_rt*((dx_W*I_dt - u(i,j,k)) + amrex::min(0.0_rt, u(i-1,j,k))));
                        du_min_CFL_val = amrex::max(du_min_CFL_val,
                            0.499_rt*((-dx_E*I_dt - u(i,j,k)) + amrex::max(0.0_rt, u(i+1,j,k))));
                    }
                } else {
                    for (int k = kmin; k <= kmax; ++k) {
                        du_max_CFL_val = amrex::min(du_max_CFL_val, dx_W*CFL_dt - u(i,j,k));
                        du_min_CFL_val = amrex::max(du_min_CFL_val, -(dx_E*CFL_dt + u(i,j,k)));
                    }
                }
            }
            du_max_CFL_val = amrex::max(du_max_CFL_val, 0.0_rt);
            du_min_CFL_val = amrex::min(du_min_CFL_val, 0.0_rt);

            uh_tot_0(i,j,0)    = uh_tot_0_val;
            duhdu_tot_0(i,j,0) = duhdu_tot_0_val;
            du_max_CFL(i,j,0)  = du_max_CFL_val;
            du_min_CFL(i,j,0)  = du_min_CFL_val;
            visc_rem_max(i,j,0) = visc_rem_max_val;
            // do_I is only ever masked false by OBC (specified/Flather segment)
            // logic, out of scope until OceanOBC is implemented in C++.
            do_I(i,j,0) = 1;
        }
    });

    if (need_bt) {
        if (uhbt.p != nullptr) {
            // Find du and uh.
            MOM::zonal_flux_adjust(bxC, u, h_in, h_W, h_E,
                                   uh_tot_0, duhdu_tot_0, du,
                                   du_max_CFL, du_min_CFL, dt,
                                   dy_Cu, IareaT, IdxT, CS,
                                   visc_rem_u_tmp, do_I,
                                   por_face_areaU, uhbt, uh, obc);

            if (has_u_cor || has_du_cor) {
                ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
                {
                    if (has_u_cor) {
                        for (int k = kmin; k <= kmax; ++k) {
                            u_cor(i,j,k) = u(i,j,k) + du(i,j,0) * visc_rem_u_tmp(i,j,k);
                        }
                        /*
                        // untested (Fortran source's own comment)!
                        if (any_simple_OBC && simple_OBC_pt(i,j)) {
                            for (int k = kmin; k <= kmax; ++k) {
                                u_cor(i,j,k) = obc->segment[amrex::Math::abs(obc->segnum_u(i,j,0))]
                                                  .normal_vel(i,j,k);
                            }
                        }
                        */
                    }
                    if (has_du_cor) {
                        du_cor(i,j,0) = du(i,j,0);
                    }
                });
            }
        }

        if (set_BT_cont) {
            // Diagnose the zero-transport correction, du0.
            MOM::zonal_flux_adjust(bxC, u, h_in, h_W, h_E,
                                   uh_tot_0, duhdu_tot_0, du,
                                   du_max_CFL, du_min_CFL, dt,
                                   dy_Cu, IareaT, IdxT, CS,
                                   visc_rem_u_tmp, do_I,
                                   por_face_areaU, Array4<const Real>{}, Array4<Real>{}, obc);

            MOM::set_zonal_BT_cont(bxC, u, h_in, h_W, h_E,
                                   FA_u_W0, FA_u_E0, FA_u_WW, FA_u_EE, uBT_WW, uBT_EE,
                                   du, uh_tot_0, duhdu_tot_0,
                                   du_max_CFL, du_min_CFL, dt,
                                   dxCu, dy_Cu, IareaT, IdxT, CS,
                                   visc_rem_u_tmp, visc_rem_max,
                                   do_I, por_face_areaU);

            /*
            // untested (Fortran source's own comment)!
            if (any_simple_OBC) {
                ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
                {
                    if (simple_OBC_pt(i,j)) {
                        Real FAuI = H_subroundoff * dy_Cu(i,j,0);
                        for (int k = kmin; k <= kmax; ++k) {
                            int const l_seg = amrex::Math::abs(obc->segnum_u(i,j,0));
                            if ((amrex::Math::abs(obc->segment[l_seg].normal_vel(i,j,k)) > 0.0_rt) &&
                                obc->segment[l_seg].specified) {
                                FAuI += obc->segment[l_seg].normal_trans(i,j,k) /
                                        obc->segment[l_seg].normal_vel(i,j,k);
                            }
                        }
                        FA_u_W0(i,j,0) = FAuI; FA_u_E0(i,j,0) = FAuI;
                        FA_u_WW(i,j,0) = FAuI; FA_u_EE(i,j,0) = FAuI;
                        uBT_WW(i,j,0) = 0.0_rt; uBT_EE(i,j,0) = 0.0_rt;
                    }
                });
            }
            */
        }
    }

    /*
    // untested (Fortran source's own comment)!
    if (local_open_BC && set_BT_cont) {
        for (int n = 0; n < obc->number_of_segments; ++n) {
            if (obc->segment[n].open && obc->segment[n].is_E_or_W) {
                int const i = obc->segment[n].HI.IsdB;
                bool const is_E = (obc->segment[n].direction == OBC_DIRECTION_E);
                ParallelFor(Box(IntVect(i, obc->segment[n].HI.Jsd, kmin),
                                IntVect(i, obc->segment[n].HI.Jed, kmax)),
                    [=] AMREX_GPU_DEVICE (int ii, int j, int) noexcept
                {
                    Real FA_u = 0.0_rt;
                    for (int k = kmin; k <= kmax; ++k) {
                        FA_u += (is_E ? h_in(ii,j,k) : h_in(ii+1,j,k)) * (dy_Cu(ii,j,0) * por_face_areaU(ii,j,k));
                    }
                    FA_u_W0(ii,j,0) = FA_u; FA_u_E0(ii,j,0) = FA_u;
                    FA_u_WW(ii,j,0) = FA_u; FA_u_EE(ii,j,0) = FA_u;
                    uBT_WW(ii,j,0) = 0.0_rt; uBT_EE(ii,j,0) = 0.0_rt;
                });
            }
        }
    }
    */

    // NOTE: BT_cont%h_u (a zonal_flux_thickness call) is out of scope here --
    // turbotmp_zonal_mass_flux_bridge does not currently expose an h_u
    // channel, so there is nothing on the C side to write that result into.
}

//> Calculates the mass or volume fluxes through the meridional faces, and
//  other related quantities -- including, optionally, the barotropic
//  mass-flux correction (v_cor/dv_cor) and the barotropic-consistency
//  face-area/velocity-correction diagnostics (FA_v_*/vBT_*).
void meridional_mass_flux(
    const Box& bxC,                          //!< Iteration box for continuity solver
    Array4<const Real> const& v,             //!< Meridional velocity [L T-1 ~> m s-1]
    Array4<const Real> const& h_in,          //!< Layer thickness used to calculate fluxes [H ~> m or kg m-2]
    Array4<const Real> const& h_S,           //!< South edge thickness in the reconstruction [H ~> m or kg m-2]
    Array4<const Real> const& h_N,           //!< North edge thickness in the reconstruction [H ~> m or kg m-2]
    Array4<Real> const& vh,                  //!< Volume flux through meridional faces, v*h*dx
                                              //!< [H L2 T-1 ~> m3 s-1 or kg s-1]
    Real dt,                                 //!< Time increment [T ~> s]
    Array4<const Real> const& dx_Cv,         //!< The grid cell's unblocked lengths of the v-faces
                                              //!< of the h-cell [L ~> m]
    Array4<const Real> const& IareaT,        //!< The grid cell's 1/areaT [L-2 ~> m-2]
    Array4<const Real> const& IdyT,          //!< The grid cell's 1/dyT [L-1 ~> m-1]
    Array4<const Real> const& areaT,         //!< The area of the h-cell [L2 ~> m2]
    Array4<const Real> const& dyT,           //!< The y-extent of the h-cell [L ~> m]
    Array4<const Real> const& mask2dCv,      //!< 0 for land points, 1 for ocean points at v-locations [nondim]
    Array4<const Real> const& dyCv,          //!< The grid cell's v-point y-extent [L ~> m]
    int isd,                                 //!< Start i-index of the data domain (0-based)
    int ied,                                 //!< End i-index of the data domain (0-based)
    Real H_subroundoff,                      //!< A negligibly small thickness used to avoid division
                                              //!< by zero [H ~> m or kg m-2]
    const transport_adjust_CS_C& CS,         //!< Options controlling the transport adjustment
                                              //!< and barotropic-consistency iteration
    OceanOBC* obc,                           //!< Open boundary control structure
    Array4<const Real> const& por_face_areaV,//!< Fractional open area of V-faces [nondim]
    Array4<const Real> const& vhbt,          //!< Summed volume flux through meridional faces
                                              //!< [H L2 T-1 ~> m3 s-1 or kg s-1]; may be absent (.p == nullptr)
    Array4<const Real> const& visc_rem_v,    //!< Fraction of momentum/barotropic acceleration remaining
                                              //!< after viscosity [nondim]; may be absent (.p == nullptr)
    Array4<Real> const& v_cor,               //!< Meridional velocity with barotropic correction [L T-1 ~> m s-1];
                                              //!< may be absent (.p == nullptr)
    Array4<Real> const& FA_v_S0,             //!< Effective open face area, south, 0 transport;
                                              //!< may be absent (.p == nullptr)
    Array4<Real> const& FA_v_N0,             //!< Effective open face area, north, 0 transport;
                                              //!< may be absent (.p == nullptr)
    Array4<Real> const& FA_v_SS,             //!< Effective open face area, southerly test velocity;
                                              //!< may be absent (.p == nullptr)
    Array4<Real> const& FA_v_NN,             //!< Effective open face area, northerly test velocity;
                                              //!< may be absent (.p == nullptr)
    Array4<Real> const& vBT_SS,              //!< Southerly correction to the barotropic velocity;
                                              //!< may be absent (.p == nullptr)
    Array4<Real> const& vBT_NN,              //!< Northerly correction to the barotropic velocity;
                                              //!< may be absent (.p == nullptr)
    Array4<Real> const& dv_cor)              //!< Meridional velocity increment from v that gives vhbt as the
                                              //!< depth-integrated transport [L T-1 ~> m s-1];
                                              //!< may be absent (.p == nullptr)
{
    BL_PROFILE("meridional_mass_flux");

    // NOTE: OBC support temporarily disabled.
    // OceanOBC is forward-declared only.
    // All boundary-condition logic removed for initial port validation.
    if (obc != nullptr) {
       AMREX_ABORT_LOC("OBC pointer provided but not yet implemented");
    }
    /*
    bool local_specified_BC = false, local_Flather_OBC = false, local_open_BC = false;
    if (obc != nullptr) {
        if (obc->OBC_pe) {
            local_specified_BC = obc->specified_v_BCs_exist_globally;
            local_Flather_OBC  = obc->Flather_v_BCs_exist_globally;
            local_open_BC      = obc->open_v_BCs_exist_globally;
        }
    }
    */

    const bool use_visc_rem = (visc_rem_v.p != nullptr);
    const bool set_BT_cont  = (FA_v_S0.p != nullptr);
    const bool need_bt      = (vhbt.p != nullptr) || set_BT_cont;
    const bool has_v_cor    = (v_cor.p != nullptr);
    const bool has_dv_cor   = (dv_cor.p != nullptr);

    const int kmin = bxC.smallEnd(2);
    const int kmax = bxC.bigEnd(2);

    Real CFL_dt = CS.CFL_limit_adjust / dt;
    const Real I_dt = 1.0_rt / dt;
    if (CS.aggress_adjust) CFL_dt = I_dt;

    // Iteration box for v-point (V-grid) fields: grown by 1 at the lower y-extent
    Box bxV = growLo(bxC, 1, 1);
    Box bx2d(IntVect(bxV.smallEnd(0), bxV.smallEnd(1), 0),
             IntVect(bxV.bigEnd(0),   bxV.bigEnd(1),   0));

    // visc_rem_v_tmp is filled over the full data-domain i-range [isd,ied]
    // (matching the Fortran source), wider than this box's own i-range.
    Box bx_visc_wide(IntVect(isd,               bxV.smallEnd(1), kmin),
                      IntVect(ied,               bxV.bigEnd(1),   kmax));

    // Scratch, internal to this orchestrator -- never crosses the Fortran/C
    // boundary. Mirrors the Fortran source's own local allocView()/free()
    // calls (see lessons.md's backlog note on an alternative, scratch-free
    // design deferred for later).
    FArrayBox visc_rem_v_tmp_fab(bx_visc_wide, 1, amrex::The_Arena());
    FArrayBox vh_tot_0_fab(bx2d, 1, amrex::The_Arena());
    FArrayBox dvhdv_tot_0_fab(bx2d, 1, amrex::The_Arena());
    FArrayBox dv_max_CFL_fab(bx2d, 1, amrex::The_Arena());
    FArrayBox dv_min_CFL_fab(bx2d, 1, amrex::The_Arena());
    FArrayBox visc_rem_max_fab(bx2d, 1, amrex::The_Arena());
    FArrayBox dv_fab(bx2d, 1, amrex::The_Arena());
    IArrayBox do_I_fab(bx2d, 1, amrex::The_Arena());

    Array4<Real> const visc_rem_v_tmp = visc_rem_v_tmp_fab.array();
    Array4<Real> const vh_tot_0        = vh_tot_0_fab.array();
    Array4<Real> const dvhdv_tot_0     = dvhdv_tot_0_fab.array();
    Array4<Real> const dv_max_CFL      = dv_max_CFL_fab.array();
    Array4<Real> const dv_min_CFL      = dv_min_CFL_fab.array();
    Array4<Real> const visc_rem_max    = visc_rem_max_fab.array();
    Array4<Real> const dv              = dv_fab.array();
    Array4<int>  const do_I            = do_I_fab.array();

    ParallelFor(bx_visc_wide, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        visc_rem_v_tmp(i,j,k) = use_visc_rem ? visc_rem_v(i,j,k) : 1.0_rt;
    });

    ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        if (has_dv_cor) dv_cor(i,j,0) = 0.0_rt;

        // Set vh (and, only when needed below, the running sums of its
        // velocity derivative and of visc_rem_v_tmp's running max).
        Real vh_tot_0_val = 0.0_rt, dvhdv_tot_0_val = 0.0_rt;
        Real visc_rem_max_running = 0.0_rt;
        for (int k = kmin; k <= kmax; ++k) {
            Real const vrt = visc_rem_v_tmp(i,j,k);

            Real dvhdv_val;
            flux_elem_point(v(i,j,k), h_in(i,j,k), h_in(i,j+1,k), h_S(i,j,k), h_S(i,j+1,k),
                            h_N(i,j,k), h_N(i,j+1,k), vh(i,j,k), dvhdv_val, vrt,
                            dx_Cv(i,j,0), IareaT(i,j,0), IareaT(i,j+1,0), IdyT(i,j,0), IdyT(i,j+1,0),
                            dt, CS.vol_CFL, por_face_areaV(i,j,k));
            /*
            // untested (Fortran source's own comment)!
            if (local_open_BC) {
                int const l_seg = obc->segnum_v(i,j,0);
                flux_elem_obc_point(v(i,j,k), h_in(i,j,k), h_in(i,j+1,k), vh(i,j,k), dvhdv_val,
                                    vrt, por_face_areaV(i,j,k), dx_Cv(i,j,0), obc, l_seg);
            }
            // untested (Fortran source's own comment)!
            if (local_specified_BC && obc->segnum_v(i,j,0) != 0) {
                int const l_seg = amrex::Math::abs(obc->segnum_v(i,j,0));
                if (obc->segment[l_seg].specified) vh(i,j,k) = obc->segment[l_seg].normal_trans(i,j,k);
            }
            */

            if (need_bt) {
                vh_tot_0_val += vh(i,j,k);
                dvhdv_tot_0_val += dvhdv_val;
                if (use_visc_rem && CS.use_visc_rem_max) {
                    visc_rem_max_running = amrex::max(visc_rem_max_running, vrt);
                }
            }
        }

        if (need_bt) {
            Real const visc_rem_max_val = (use_visc_rem && CS.use_visc_rem_max) ? visc_rem_max_running : 1.0_rt;
            Real const I_vrm = (visc_rem_max_val > 0.0_rt) ? 1.0_rt / visc_rem_max_val : 0.0_rt;

            //   Set limits on dv that will keep the CFL number between -1 and 1.
            // This should be adequate to keep the root bracketed in all cases.
            Real dy_S, dy_N;
            if (CS.vol_CFL) {
                dy_S = ratio_max_point(areaT(i,j,0), dx_Cv(i,j,0), 1000.0_rt*dyT(i,j,0));
                dy_N = ratio_max_point(areaT(i,j+1,0), dx_Cv(i,j,0), 1000.0_rt*dyT(i,j+1,0));
            } else {
                dy_S = dyT(i,j,0);
                dy_N = dyT(i,j+1,0);
            }
            Real dv_max_CFL_val = 2.0_rt * (CFL_dt * dy_S) * I_vrm;
            Real dv_min_CFL_val = -2.0_rt * (CFL_dt * dy_N) * I_vrm;

            if (use_visc_rem) {
                if (CS.aggress_adjust) {
                    // untested (Fortran source's own comment)!
                    for (int k = kmin; k <= kmax; ++k) {
                        Real const vrt = visc_rem_v_tmp(i,j,k);
                        Real const dv_lim_max = 0.499_rt * ((dy_S*I_dt - v(i,j,k)) + amrex::min(0.0_rt, v(i,j-1,k)));
                        if (dv_max_CFL_val * vrt > dv_lim_max) dv_max_CFL_val = dv_lim_max / vrt;
                        // Fortran source uses CFL_dt (not I_dt) on this line; the two are
                        // guaranteed equal here since CS.aggress_adjust forces CFL_dt = I_dt.
                        Real const dv_lim_min = 0.499_rt * ((-dy_N*CFL_dt - v(i,j,k)) + amrex::max(0.0_rt, v(i,j+1,k)));
                        if (dv_min_CFL_val * vrt < dv_lim_min) dv_min_CFL_val = dv_lim_min / vrt;
                    }
                } else {
                    for (int k = kmin; k <= kmax; ++k) {
                        Real const vrt = visc_rem_v_tmp(i,j,k);
                        if (dv_max_CFL_val * vrt > dy_S*CFL_dt - v(i,j,k)*mask2dCv(i,j,0))
                            dv_max_CFL_val = (dy_S*CFL_dt - v(i,j,k)) / vrt;
                        if (dv_min_CFL_val * vrt < -dy_N*CFL_dt - v(i,j,k)*mask2dCv(i,j,0))
                            dv_min_CFL_val = -(dy_N*CFL_dt + v(i,j,k)) / vrt;
                    }
                }
            } else {
                if (CS.aggress_adjust) {
                    // untested (Fortran source's own comment)!
                    for (int k = kmin; k <= kmax; ++k) {
                        dv_max_CFL_val = amrex::min(dv_max_CFL_val,
                            0.499_rt*((dy_S*I_dt - v(i,j,k)) + amrex::min(0.0_rt, v(i,j-1,k))));
                        dv_min_CFL_val = amrex::max(dv_min_CFL_val,
                            0.499_rt*((-dy_N*I_dt - v(i,j,k)) + amrex::max(0.0_rt, v(i,j+1,k))));
                    }
                } else {
                    for (int k = kmin; k <= kmax; ++k) {
                        dv_max_CFL_val = amrex::min(dv_max_CFL_val, dy_S*CFL_dt - v(i,j,k));
                        dv_min_CFL_val = amrex::max(dv_min_CFL_val, -(dy_N*CFL_dt + v(i,j,k)));
                    }
                }
            }
            dv_max_CFL_val = amrex::max(dv_max_CFL_val, 0.0_rt);
            dv_min_CFL_val = amrex::min(dv_min_CFL_val, 0.0_rt);

            vh_tot_0(i,j,0)    = vh_tot_0_val;
            dvhdv_tot_0(i,j,0) = dvhdv_tot_0_val;
            dv_max_CFL(i,j,0)  = dv_max_CFL_val;
            dv_min_CFL(i,j,0)  = dv_min_CFL_val;
            visc_rem_max(i,j,0) = visc_rem_max_val;
            // do_I is only ever masked false by OBC (specified/Flather segment)
            // logic, out of scope until OceanOBC is implemented in C++.
            do_I(i,j,0) = 1;
        }
    });

    if (need_bt) {
        if (vhbt.p != nullptr) {
            // Find dv and vh.
            MOM::meridional_flux_adjust(bxC, v, h_in, h_S, h_N,
                                        vh_tot_0, dvhdv_tot_0, dv,
                                        dv_max_CFL, dv_min_CFL, dt,
                                        dx_Cv, IareaT, IdyT, CS,
                                        visc_rem_v_tmp, do_I,
                                        por_face_areaV, vhbt, vh, obc);

            if (has_v_cor || has_dv_cor) {
                ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
                {
                    if (has_v_cor) {
                        for (int k = kmin; k <= kmax; ++k) {
                            v_cor(i,j,k) = v(i,j,k) + dv(i,j,0) * visc_rem_v_tmp(i,j,k);
                        }
                        /*
                        // untested (Fortran source's own comment)!
                        if (any_simple_OBC && simple_OBC_pt(i,j)) {
                            for (int k = kmin; k <= kmax; ++k) {
                                v_cor(i,j,k) = obc->segment[amrex::Math::abs(obc->segnum_v(i,j,0))]
                                                  .normal_vel(i,j,k);
                            }
                        }
                        */
                    }
                    if (has_dv_cor) {
                        dv_cor(i,j,0) = dv(i,j,0);
                    }
                });
            }
        }

        if (set_BT_cont) {
            // Diagnose the zero-transport correction, dv0.
            MOM::meridional_flux_adjust(bxC, v, h_in, h_S, h_N,
                                        vh_tot_0, dvhdv_tot_0, dv,
                                        dv_max_CFL, dv_min_CFL, dt,
                                        dx_Cv, IareaT, IdyT, CS,
                                        visc_rem_v_tmp, do_I,
                                        por_face_areaV, Array4<const Real>{}, Array4<Real>{}, obc);

            MOM::set_merid_BT_cont(bxC, v, h_in, h_S, h_N,
                                   FA_v_S0, FA_v_N0, FA_v_SS, FA_v_NN, vBT_SS, vBT_NN,
                                   dv, vh_tot_0, dvhdv_tot_0,
                                   dv_max_CFL, dv_min_CFL, dt,
                                   dyCv, dx_Cv, IareaT, IdyT, CS,
                                   visc_rem_v_tmp, visc_rem_max,
                                   do_I, por_face_areaV);

            /*
            // untested (Fortran source's own comment)!
            // NOTE: the Fortran source's own i-range here reads i=ish:jeh (not
            // ieh) -- ported verbatim; likely a pre-existing Fortran-side typo.
            if (any_simple_OBC) {
                ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
                {
                    if (simple_OBC_pt(i,j)) {
                        Real FAvi = H_subroundoff * dx_Cv(i,j,0);
                        for (int k = kmin; k <= kmax; ++k) {
                            int const l_seg = amrex::Math::abs(obc->segnum_v(i,j,0));
                            if ((amrex::Math::abs(obc->segment[l_seg].normal_vel(i,j,k)) > 0.0_rt) &&
                                obc->segment[l_seg].specified) {
                                FAvi += obc->segment[l_seg].normal_trans(i,j,k) /
                                        obc->segment[l_seg].normal_vel(i,j,k);
                            }
                        }
                        FA_v_S0(i,j,0) = FAvi; FA_v_N0(i,j,0) = FAvi;
                        FA_v_SS(i,j,0) = FAvi; FA_v_NN(i,j,0) = FAvi;
                        vBT_SS(i,j,0) = 0.0_rt; vBT_NN(i,j,0) = 0.0_rt;
                    }
                });
            }
            */
        }
    }

    /*
    // untested (Fortran source's own comment)!
    if (local_open_BC && set_BT_cont) {
        for (int n = 0; n < obc->number_of_segments; ++n) {
            if (obc->segment[n].open && obc->segment[n].is_N_or_S) {
                int const j = obc->segment[n].HI.JsdB;
                bool const is_N = (obc->segment[n].direction == OBC_DIRECTION_N);
                ParallelFor(Box(IntVect(obc->segment[n].HI.Isd, j, kmin),
                                IntVect(obc->segment[n].HI.Ied, j, kmax)),
                    [=] AMREX_GPU_DEVICE (int i, int jj, int) noexcept
                {
                    Real FA_v = 0.0_rt;
                    for (int k = kmin; k <= kmax; ++k) {
                        FA_v += (is_N ? h_in(i,j,k) : h_in(i,j+1,k)) * (dx_Cv(i,j,0) * por_face_areaV(i,j,k));
                    }
                    FA_v_S0(i,j,0) = FA_v; FA_v_N0(i,j,0) = FA_v;
                    FA_v_SS(i,j,0) = FA_v; FA_v_NN(i,j,0) = FA_v;
                    vBT_SS(i,j,0) = 0.0_rt; vBT_NN(i,j,0) = 0.0_rt;
                });
            }
        }
    }
    */

    // NOTE: BT_cont%h_v (a meridional_flux_thickness call) is out of scope
    // here -- turbotmp_meridional_mass_flux_bridge does not currently expose
    // an h_v channel, so there is nothing on the C side to write that
    // result into.
}
}
