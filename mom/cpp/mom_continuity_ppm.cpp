// mom_continuity_ppm.cpp
/**
 * @file mom_continuity_ppm.cpp
 * @brief Box-level AMReX kernel implementations for MOM6 PPM continuity.
 */
/// @brief Abort with @p msg annotated by the source file and line number.
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
    const Box& bxC,
    Array4<const Real> const& u,
    Array4<const Real> const& h_in,
    Array4<const Real> const& h_W,
    Array4<const Real> const& h_E,
    Array4<const Real> const& uh_tot_0,
    Array4<const Real> const& duhdu_tot_0,
    Array4<Real> const& du,
    Array4<const Real> const& du_max_CFL,
    Array4<const Real> const& du_min_CFL,
    Real dt,
    Array4<const Real> const& dy_Cu,
    Array4<const Real> const& IareaT,
    Array4<const Real> const& IdxT,
    const transport_adjust_CS_C& CS,
    Array4<const Real> const& visc_rem,
    Array4<const int> const& do_I_in,
    Array4<const Real> const& por_face_areaU,
    Array4<const Real> const& uhbt,
    Array4<Real> const& uh_3d,
    OceanOBC* obc)
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
    const Box& bxC,
    Array4<const Real> const& v,
    Array4<const Real> const& h_in,
    Array4<const Real> const& h_S,
    Array4<const Real> const& h_N,
    Array4<const Real> const& vh_tot_0,
    Array4<const Real> const& dvhdv_tot_0,
    Array4<Real> const& dv,
    Array4<const Real> const& dv_max_CFL,
    Array4<const Real> const& dv_min_CFL,
    Real dt,
    Array4<const Real> const& dx_Cv,
    Array4<const Real> const& IareaT,
    Array4<const Real> const& IdyT,
    const transport_adjust_CS_C& CS,
    Array4<const Real> const& visc_rem,
    Array4<const int> const& do_I_in,
    Array4<const Real> const& por_face_areaV,
    Array4<const Real> const& vhbt,
    Array4<Real> const& vh_3d,
    OceanOBC* obc)
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
}
