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
}
