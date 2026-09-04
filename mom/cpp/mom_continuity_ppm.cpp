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
}
