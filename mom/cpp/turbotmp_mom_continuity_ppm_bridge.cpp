/**
 * @file turbotmp_mom_continuity_ppm_bridge.cpp
 * @brief Bridge that moves data (host to device, Fortran to C++ array order, and
 *        Box setup) between the MOM6 Fortran shim and the AMReX PPM continuity kernels.
 */
#include "mom_continuity_ppm.hpp"
#include "turbotmp_helper.hpp"
#include "turbotmp_mom_continuity_ppm_bridge.h"
#include <AMReX_Print.H>
#include <fstream>
#include <string>

using namespace amrex;


namespace {
bool verbose = false;
}

/**
 * @brief Bridge for the function PPM_limit_pos function
 *
 * This function acts as a bridge between a Fortran interface
 * and an AMReX C++ implementation. It also provides the ability
 * to either capture the input, or output or execute the AMReX C++ 
 * implementation based on the setting of the @p mode parameter.
 *
 * @param bx_HOST   Box over which to iterate
 * @param h_in_HOST Layer thickness [H → m or kg m^-2]
 * 	on the host in Fortran order
 * @param h_L_HOST Left thickness of the reconstruction {host, Fortran order} 
 * 	[H → m or kg m^-2]
 * @param h_R_HOST Right thickness in the reconstruction {host, Fortran order} 
 * 	[H → m or kg m^-2] 
 * @param h_min Minimum thickness allowed by the parabolic fit (host, Fortran order) 
 * 	[H → m or kg m^-2]
 *
 * @note On return, @p h_L_HOST and @p h_R_HOST hold the modified thickness values.
 */
void turbotmp_ppm_limit_pos_bridge(const Box_C* bx_HOST,
		                   const RealArray_C* h_in_HOST,
			           RealArray_C* h_L_HOST,
			           RealArray_C* h_R_HOST,
	                           const double h_min)
{ 
    /// Define Active domain (kernel launch only on real cells)
    Box bx(IntVect(bx_HOST->idxS[0]-1, bx_HOST->idxS[1]-1, bx_HOST->idxS[2]-1),
	   IntVect(bx_HOST->idxE[0]-1, bx_HOST->idxE[1]-1, bx_HOST->idxE[2]-1));

    /// Create A4 containers for the Fortran arrays
    auto h_in_DEV = turbotmp::make_array4(h_in_HOST->shape[0], h_in_HOST->shape[1], h_in_HOST->shape[2], 1, h_in_HOST->lb[0], h_in_HOST->lb[1], h_in_HOST->lb[2]);
    auto h_L_DEV  = turbotmp::make_array4(h_L_HOST->shape[0], h_L_HOST->shape[1], h_L_HOST->shape[2], 1, h_L_HOST->lb[0], h_L_HOST->lb[1], h_L_HOST->lb[2]);
    auto h_R_DEV  = turbotmp::make_array4(h_R_HOST->shape[0], h_R_HOST->shape[1], h_R_HOST->shape[2], 1, h_R_HOST->lb[0], h_R_HOST->lb[1], h_R_HOST->lb[2]);

    /// Copy from Fortran arrays to A4 container
    turbotmp::copy_FortranHost_to_array4(h_in_HOST->data, h_in_DEV);
    turbotmp::copy_FortranHost_to_array4(h_L_HOST->data, h_L_DEV);
    turbotmp::copy_FortranHost_to_array4(h_R_HOST->data, h_R_DEV);

    if(verbose) amrex::Print() << "Entered turbotmp_ppm_limit_pos_bridge\n";
    ///-------------------------------------------------
    ///  Execute kernel
    ///-------------------------------------------------
    MOM::ppm_limit_pos(bx,h_in_DEV.arr, h_L_DEV.arr, h_R_DEV.arr, h_min);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device → host
    turbotmp::copy_array4_to_FortranHost(h_L_DEV, h_L_HOST->data);
    turbotmp::copy_array4_to_FortranHost(h_R_DEV, h_R_HOST->data);

    /// Free a4 container
    turbotmp::free_array4(h_in_DEV);
    turbotmp::free_array4(h_R_DEV);
    turbotmp::free_array4(h_L_DEV);
}

/**
 * @brief Bridge for the function PPM_limit_cw84 function
 *
 * This function acts as a bridge between a Fortran interface
 * and an AMReX C++ implementation. It also provides the ability
 * to either capture the input, or output or execute the AMReX C++
 * implementation based on the setting of the @p mode parameter.
 *
 * @param bx_HOST   Box over which to iterate 
 * @param h_in_HOST Layer thickness [H → m or kg m^-2]
 *      on the host in Fortran order
 * @param h_L_HOST Left thickness of the reconstruction {host, Fortran order}
 *      [H → m or kg m^-2]
 * @param h_R_HOST Right thickness in the reconstruction {host, Fortran order}
 *      [H → m or kg m^-2]
 *
 * @note On return, @p h_L_HOST and @p h_R_HOST hold the modified thickness values.
 */
void turbotmp_ppm_limit_cw84_bridge(const Box_C* bx_HOST,
	                  const RealArray_C* h_in_HOST,
                          RealArray_C* h_L_HOST,
                          RealArray_C* h_R_HOST)
{
    /// Define Active domain (kernel launch only on real cells)
    Box bx(IntVect(bx_HOST->idxS[0]-1, bx_HOST->idxS[1]-1, bx_HOST->idxS[2]-1),
           IntVect(bx_HOST->idxE[0]-1, bx_HOST->idxE[1]-1, bx_HOST->idxE[2]-1));

    /// Create A4 containers for the Fortran arrays
    auto h_in_DEV = turbotmp::make_array4(h_in_HOST->shape[0], h_in_HOST->shape[1], h_in_HOST->shape[2], 1, h_in_HOST->lb[0], h_in_HOST->lb[1], h_in_HOST->lb[2]);
    auto h_L_DEV  = turbotmp::make_array4(h_L_HOST->shape[0], h_L_HOST->shape[1], h_L_HOST->shape[2], 1, h_L_HOST->lb[0], h_L_HOST->lb[1], h_L_HOST->lb[2]);
    auto h_R_DEV  = turbotmp::make_array4(h_R_HOST->shape[0], h_R_HOST->shape[1], h_R_HOST->shape[2], 1, h_R_HOST->lb[0], h_R_HOST->lb[1], h_R_HOST->lb[2]);

    /// Copy from Fortran arrays to A4 container
    turbotmp::copy_FortranHost_to_array4(h_in_HOST->data, h_in_DEV);
    turbotmp::copy_FortranHost_to_array4(h_L_HOST->data, h_L_DEV);
    turbotmp::copy_FortranHost_to_array4(h_R_HOST->data, h_R_DEV);

    if(verbose) amrex::Print() << "Entered turbotmp_ppm_limit_cw84_bridge\n";
    ///-------------------------------------------------
    ///  Execute kernel
    ///-------------------------------------------------
    MOM::ppm_limit_cw84(bx,h_in_DEV.arr, h_L_DEV.arr, h_R_DEV.arr);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device → host
    turbotmp::copy_array4_to_FortranHost(h_L_DEV, h_L_HOST->data);
    turbotmp::copy_array4_to_FortranHost(h_R_DEV, h_R_HOST->data);

    /// Free memory from a4 containers
    turbotmp::free_array4(h_in_DEV);
    turbotmp::free_array4(h_R_DEV);
    turbotmp::free_array4(h_L_DEV);
}

/**
 * @brief Bridge for the function PPM_reconstruction_y
 *
 * This function acts as a bridge between a Fortran interface
 * and an AMReX C++ implementation. It also provides the ability
 * to either capture the input, or output or execute the AMReX C++
 * implementation based on the setting of the @p mode parameter.
 *
 * @param bx_HOST        Box over which to iterate
 * @param h_in_HOST      Layer thickness [H → m or kg m^-2] (host, Fortran order)
 * @param h_S_HOST       South edge thickness (host, Fortran order)
 * @param h_N_HOST       North edge thickness (host, Fortran order)
 * @param mask2dT_HOST   Mask (0 land, 1 ocean) (host, Fortran order)
 * @param h_min       Minimum thickness
 * @param monotonic   Use CW84 limiter if true
 * @param simple_2nd  Use simple 2nd order scheme if true
 * @param obc         Open boundary control structure
 *
 * @note On return, @p h_S_HOST and @p h_N_HOST hold the modified thickness values.
 */
void turbotmp_ppm_reconstruction_y_bridge(const Box_C* bx_HOST,
                                          const RealArray_C* h_in_HOST,
                                          RealArray_C* h_S_HOST,
                                          RealArray_C* h_N_HOST,
                                          const RealArray_C* mask2dT_HOST,
                                          const double h_min,
					  const bool monotonic,
                                          const bool simple_2nd,
				          OceanOBC* obc)
{
    /// Define Active domain (kernel launch only on real cells)
    amrex::Box bx(amrex::IntVect(bx_HOST->idxS[0]-1, bx_HOST->idxS[1]-1, bx_HOST->idxS[2]-1),
                  amrex::IntVect(bx_HOST->idxE[0]-1, bx_HOST->idxE[1]-1, bx_HOST->idxE[2]-1));

    /// Create A4 containers for the Fortran arrays
    auto h_in_DEV    = turbotmp::make_array4(h_in_HOST->shape[0], h_in_HOST->shape[1], h_in_HOST->shape[2], 1, h_in_HOST->lb[0], h_in_HOST->lb[1], h_in_HOST->lb[2]);
    auto h_S_DEV     = turbotmp::make_array4(h_S_HOST->shape[0], h_S_HOST->shape[1], h_S_HOST->shape[2], 1, h_S_HOST->lb[0], h_S_HOST->lb[1], h_S_HOST->lb[2]);
    auto h_N_DEV     = turbotmp::make_array4(h_N_HOST->shape[0], h_N_HOST->shape[1], h_N_HOST->shape[2], 1, h_N_HOST->lb[0], h_N_HOST->lb[1], h_N_HOST->lb[2]);
    auto mask2dT_DEV = turbotmp::make_array4(mask2dT_HOST->shape[0], mask2dT_HOST->shape[1], 1, 1, mask2dT_HOST->lb[0], mask2dT_HOST->lb[1], 1);

    /// Copy from Fortran arrays to A4 container
    turbotmp::copy_FortranHost_to_array4(h_in_HOST->data,    h_in_DEV);
    turbotmp::copy_FortranHost_to_array4(h_S_HOST->data,     h_S_DEV);
    turbotmp::copy_FortranHost_to_array4(h_N_HOST->data,     h_N_DEV);
    turbotmp::copy_FortranHost_to_array4(mask2dT_HOST->data, mask2dT_DEV);

    if(verbose) amrex::Print() << "Entered turbotmp_ppm_reconstruction_y_bridge\n";
    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------

    MOM::PPM_reconstruction_y(bx,
                         h_in_DEV.arr,
                         h_S_DEV.arr,
                         h_N_DEV.arr,
                         mask2dT_DEV.arr,
                         h_min,
                         monotonic,
                         simple_2nd,
                         obc);

    /// Ensure kernel is done before copying back
    amrex::Gpu::synchronize();

    /// Copy device → host
    turbotmp::copy_array4_to_FortranHost(h_S_DEV, h_S_HOST->data);
    turbotmp::copy_array4_to_FortranHost(h_N_DEV, h_N_HOST->data);

    /// Free memory from a4 containers
    turbotmp::free_array4(h_in_DEV);
    turbotmp::free_array4(h_S_DEV);
    turbotmp::free_array4(h_N_DEV);
    turbotmp::free_array4(mask2dT_DEV);
}

/**
 * @brief Bridge for the function PPM_reconstruction_x
 *
 * @param bx_HOST        Box over which to iterate
 * @param h_in_HOST      Layer thickness [H → m or kg m^-2] (host, Fortran order)
 * @param h_W_HOST       West edge thickness (host, Fortran order)
 * @param h_E_HOST       East edge thickness (host, Fortran order)
 * @param mask2dT_HOST   Mask (0 land, 1 ocean) (host, Fortran order)
 * @param h_min       Minimum thickness
 * @param monotonic   Use CW84 limiter if true
 * @param simple_2nd  Use simple 2nd order scheme if true
 * @param obc         Open boundary control structure
 *
 * @note On return, @p h_W_HOST and @p h_E_HOST hold the modified thickness values.
 */
void turbotmp_ppm_reconstruction_x_bridge(const Box_C* bx_HOST,
                                          const RealArray_C* h_in_HOST,
                                          RealArray_C* h_W_HOST,
                                          RealArray_C* h_E_HOST,
                                          const RealArray_C* mask2dT_HOST,
                                          const double h_min,
                                          const bool monotonic,
                                          const bool simple_2nd,
                                          OceanOBC* obc)
{
    /// Define Active domain (kernel launch only on real cells)
    amrex::Box bx(amrex::IntVect(bx_HOST->idxS[0]-1, bx_HOST->idxS[1]-1, bx_HOST->idxS[2]-1),
                  amrex::IntVect(bx_HOST->idxE[0]-1, bx_HOST->idxE[1]-1, bx_HOST->idxE[2]-1));

    /// Create A4 containers for the Fortran arrays
    auto h_in_DEV    = turbotmp::make_array4(h_in_HOST->shape[0], h_in_HOST->shape[1], h_in_HOST->shape[2], 1, h_in_HOST->lb[0], h_in_HOST->lb[1], h_in_HOST->lb[2]);
    auto h_W_DEV     = turbotmp::make_array4(h_W_HOST->shape[0], h_W_HOST->shape[1], h_W_HOST->shape[2], 1, h_W_HOST->lb[0], h_W_HOST->lb[1], h_W_HOST->lb[2]);
    auto h_E_DEV     = turbotmp::make_array4(h_E_HOST->shape[0], h_E_HOST->shape[1], h_E_HOST->shape[2], 1, h_E_HOST->lb[0], h_E_HOST->lb[1], h_E_HOST->lb[2]);
    auto mask2dT_DEV = turbotmp::make_array4(mask2dT_HOST->shape[0], mask2dT_HOST->shape[1], 1, 1, mask2dT_HOST->lb[0], mask2dT_HOST->lb[1], 1);

    /// Copy host → device (h_W and h_E are inout: copy in before kernel)
    turbotmp::copy_FortranHost_to_array4(h_in_HOST->data,    h_in_DEV);
    turbotmp::copy_FortranHost_to_array4(h_W_HOST->data,     h_W_DEV);
    turbotmp::copy_FortranHost_to_array4(h_E_HOST->data,     h_E_DEV);
    turbotmp::copy_FortranHost_to_array4(mask2dT_HOST->data, mask2dT_DEV);

    if(verbose) amrex::Print() << "Entered turbotmp_ppm_reconstruction_x_bridge\n";
    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::PPM_reconstruction_x(bx,
                         h_in_DEV.arr,
                         h_W_DEV.arr,
                         h_E_DEV.arr,
                         mask2dT_DEV.arr,
                         h_min,
                         monotonic,
                         simple_2nd,
                         obc);

    /// Ensure kernel is done before copying back
    amrex::Gpu::synchronize();

    /// Copy device → host (outputs only)
    turbotmp::copy_array4_to_FortranHost(h_W_DEV, h_W_HOST->data);
    turbotmp::copy_array4_to_FortranHost(h_E_DEV, h_E_HOST->data);

    /// Free memory from a4 containers
    turbotmp::free_array4(h_in_DEV);
    turbotmp::free_array4(h_W_DEV);
    turbotmp::free_array4(h_E_DEV);
    turbotmp::free_array4(mask2dT_DEV);
}
/**
 * @brief Bridge for zonal_edge_thickness
 *
 * @param bx_HOST        Box over which to iterate
 * @param h_in_HOST      Layer thickness (host, Fortran order)
 * @param h_W_HOST       West edge thickness (host, Fortran order)
 * @param h_E_HOST       East edge thickness (host, Fortran order)
 * @param mask2dT_HOST   Mask (0 land, 1 ocean) (host, Fortran order)
 * @param h_min       Minimum thickness
 * @param upwind_1st  If true, use 1st-order upwind reconstruction
 * @param monotonic   Use CW84 limiter if true
 * @param simple_2nd  Use simple 2nd order scheme if true
 * @param obc         Open boundary control structure
 *
 * @note On return, @p h_W_HOST and @p h_E_HOST hold the modified thickness values.
 */
void turbotmp_zonal_edge_thickness_bridge(const Box_C* bx_HOST,
                                          const RealArray_C* h_in_HOST,
                                          RealArray_C* h_W_HOST,
                                          RealArray_C* h_E_HOST,
                                          const RealArray_C* mask2dT_HOST,
                                          const double h_min,
                                          const bool upwind_1st,
                                          const bool monotonic,
                                          const bool simple_2nd,
                                          OceanOBC* obc)
{
    /// Define Active domain (kernel launch only on real cells)
    amrex::Box bx(amrex::IntVect(bx_HOST->idxS[0]-1, bx_HOST->idxS[1]-1, bx_HOST->idxS[2]-1),
                  amrex::IntVect(bx_HOST->idxE[0]-1, bx_HOST->idxE[1]-1, bx_HOST->idxE[2]-1));

    /// Create A4 containers for the Fortran arrays
    auto h_in_DEV    = turbotmp::make_array4(h_in_HOST->shape[0], h_in_HOST->shape[1], h_in_HOST->shape[2], 1, h_in_HOST->lb[0], h_in_HOST->lb[1], h_in_HOST->lb[2]);
    auto h_W_DEV     = turbotmp::make_array4(h_W_HOST->shape[0], h_W_HOST->shape[1], h_W_HOST->shape[2], 1, h_W_HOST->lb[0], h_W_HOST->lb[1], h_W_HOST->lb[2]);
    auto h_E_DEV     = turbotmp::make_array4(h_E_HOST->shape[0], h_E_HOST->shape[1], h_E_HOST->shape[2], 1, h_E_HOST->lb[0], h_E_HOST->lb[1], h_E_HOST->lb[2]);
    auto mask2dT_DEV = turbotmp::make_array4(mask2dT_HOST->shape[0], mask2dT_HOST->shape[1], 1, 1, mask2dT_HOST->lb[0], mask2dT_HOST->lb[1], 1);

    /// Copy host → device (h_W and h_E are inout: copy in before kernel)
    turbotmp::copy_FortranHost_to_array4(h_in_HOST->data,    h_in_DEV);
    turbotmp::copy_FortranHost_to_array4(h_W_HOST->data,     h_W_DEV);
    turbotmp::copy_FortranHost_to_array4(h_E_HOST->data,     h_E_DEV);
    turbotmp::copy_FortranHost_to_array4(mask2dT_HOST->data, mask2dT_DEV);

    if(verbose) amrex::Print() << "Entered turbotmp_zonal_edge_thickness_bridge\n";
    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::zonal_edge_thickness(bx,
                              h_in_DEV.arr,
                              h_W_DEV.arr,
                              h_E_DEV.arr,
                              mask2dT_DEV.arr,
                              h_min,
                              upwind_1st,
                              monotonic,
                              simple_2nd,
                              obc);

    /// Ensure kernel is done before copying back
    amrex::Gpu::synchronize();

    /// Copy device → host (outputs only)
    turbotmp::copy_array4_to_FortranHost(h_W_DEV, h_W_HOST->data);
    turbotmp::copy_array4_to_FortranHost(h_E_DEV, h_E_HOST->data);

    /// Free memory from a4 containers
    turbotmp::free_array4(h_in_DEV);
    turbotmp::free_array4(h_W_DEV);
    turbotmp::free_array4(h_E_DEV);
    turbotmp::free_array4(mask2dT_DEV);
}

/**
 * @brief Bridge for meridional_edge_thickness
 *
 * @param bx_HOST        Box over which to iterate
 * @param h_in_HOST      Layer thickness (host, Fortran order)
 * @param h_S_HOST       South edge thickness (host, Fortran order)
 * @param h_N_HOST       North edge thickness (host, Fortran order)
 * @param mask2dT_HOST   Mask (0 land, 1 ocean) (host, Fortran order)
 * @param h_min       Minimum thickness
 * @param upwind_1st  If true, use 1st-order upwind reconstruction
 * @param monotonic   Use CW84 limiter if true
 * @param simple_2nd  Use simple 2nd order scheme if true
 * @param obc         Open boundary control structure
 *
 * @note On return, @p h_S_HOST and @p h_N_HOST hold the modified thickness values.
 */
void turbotmp_meridional_edge_thickness_bridge(const Box_C* bx_HOST,
                                               const RealArray_C* h_in_HOST,
                                               RealArray_C* h_S_HOST,
                                               RealArray_C* h_N_HOST,
                                               const RealArray_C* mask2dT_HOST,
                                               const double h_min,
                                               const bool upwind_1st,
                                               const bool monotonic,
                                               const bool simple_2nd,
                                               OceanOBC* obc)
{
    /// Define Active domain (kernel launch only on real cells)
    amrex::Box bx(amrex::IntVect(bx_HOST->idxS[0]-1, bx_HOST->idxS[1]-1, bx_HOST->idxS[2]-1),
                  amrex::IntVect(bx_HOST->idxE[0]-1, bx_HOST->idxE[1]-1, bx_HOST->idxE[2]-1));

    /// Create A4 containers for the Fortran arrays
    auto h_in_DEV    = turbotmp::make_array4(h_in_HOST->shape[0], h_in_HOST->shape[1], h_in_HOST->shape[2], 1, h_in_HOST->lb[0], h_in_HOST->lb[1], h_in_HOST->lb[2]);
    auto h_S_DEV     = turbotmp::make_array4(h_S_HOST->shape[0], h_S_HOST->shape[1], h_S_HOST->shape[2], 1, h_S_HOST->lb[0], h_S_HOST->lb[1], h_S_HOST->lb[2]);
    auto h_N_DEV     = turbotmp::make_array4(h_N_HOST->shape[0], h_N_HOST->shape[1], h_N_HOST->shape[2], 1, h_N_HOST->lb[0], h_N_HOST->lb[1], h_N_HOST->lb[2]);
    auto mask2dT_DEV = turbotmp::make_array4(mask2dT_HOST->shape[0], mask2dT_HOST->shape[1], 1, 1, mask2dT_HOST->lb[0], mask2dT_HOST->lb[1], 1);

    /// Copy host → device (h_S and h_N are inout: copy in before kernel)
    turbotmp::copy_FortranHost_to_array4(h_in_HOST->data,    h_in_DEV);
    turbotmp::copy_FortranHost_to_array4(h_S_HOST->data,     h_S_DEV);
    turbotmp::copy_FortranHost_to_array4(h_N_HOST->data,     h_N_DEV);
    turbotmp::copy_FortranHost_to_array4(mask2dT_HOST->data, mask2dT_DEV);

    if(verbose) amrex::Print() << "Entered: turbotmp_meridional_edge_thickness_bridge\n";
    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::meridional_edge_thickness(bx,
                                   h_in_DEV.arr,
                                   h_S_DEV.arr,
                                   h_N_DEV.arr,
                                   mask2dT_DEV.arr,
                                   h_min,
                                   upwind_1st,
                                   monotonic,
                                   simple_2nd,
                                   obc);

    /// Ensure kernel is done before copying back
    amrex::Gpu::synchronize();

    /// Copy device → host (outputs only)
    turbotmp::copy_array4_to_FortranHost(h_S_DEV, h_S_HOST->data);
    turbotmp::copy_array4_to_FortranHost(h_N_DEV, h_N_HOST->data);

    /// Free memory from a4 containers
    turbotmp::free_array4(h_in_DEV);
    turbotmp::free_array4(h_S_DEV);
    turbotmp::free_array4(h_N_DEV);
    turbotmp::free_array4(mask2dT_DEV);
}

/**
 * @brief Bridge for set_zonal_BT_cont
 *
 * @param bxC_HOST            Box over which to iterate
 * @param u_HOST              Zonal velocity (host, Fortran order)
 * @param h_in_HOST           Layer thickness used to calculate fluxes (host, Fortran order)
 * @param h_W_HOST            West edge thickness in the reconstruction (host, Fortran order)
 * @param h_E_HOST            East edge thickness in the reconstruction (host, Fortran order)
 * @param FA_u_W0_HOST        Effective open face area, west, 0 transport (host, Fortran order)
 * @param FA_u_E0_HOST        Effective open face area, east, 0 transport (host, Fortran order)
 * @param FA_u_WW_HOST        Effective open face area, westerly test velocity (host, Fortran order)
 * @param FA_u_EE_HOST        Effective open face area, easterly test velocity (host, Fortran order)
 * @param uBT_WW_HOST         Westerly correction to the barotropic velocity (host, Fortran order)
 * @param uBT_EE_HOST         Easterly correction to the barotropic velocity (host, Fortran order)
 * @param du0_HOST            Barotropic velocity increment giving 0 transport (host, Fortran order)
 * @param uh_tot_0_HOST       Summed transport with 0 adjustment (host, Fortran order)
 * @param duhdu_tot_0_HOST    Partial derivative of du_err with du at 0 adjustment (host, Fortran order)
 * @param du_max_CFL_HOST     Maximum acceptable value of du (host, Fortran order)
 * @param du_min_CFL_HOST     Minimum acceptable value of du (host, Fortran order)
 * @param dt                  Time increment
 * @param dxCu_HOST           The grid cell's u-point x-extent (host, Fortran order)
 * @param dy_Cu_HOST          Unblocked u-face length, 2D (host, Fortran order)
 * @param IareaT_HOST         1/areaT, 2D (host, Fortran order)
 * @param IdxT_HOST           1/dxT, 2D (host, Fortran order)
 * @param CS_HOST             Transport-adjustment and barotropic-consistency options
 * @param visc_rem_HOST       Fraction of momentum/barotropic acceleration remaining after
 *                            viscosity (host, Fortran order) [nondim]
 * @param visc_rem_max_HOST   Maximum allowable viscosity remnant (host, Fortran order)
 * @param do_I_HOST           Logical flag indicating which I values to work on (host, Fortran order)
 * @param por_face_areaU_HOST Fractional open area of U-faces (host, Fortran order)
 *
 * @note On return, @p FA_u_W0_HOST, @p FA_u_E0_HOST, @p FA_u_WW_HOST, @p FA_u_EE_HOST,
 *       @p uBT_WW_HOST, and @p uBT_EE_HOST hold the modified values.
 */
void turbotmp_set_zonal_bt_cont_bridge(const Box_C* bxC_HOST,
                                       const RealArray_C* u_HOST,
                                       const RealArray_C* h_in_HOST,
                                       const RealArray_C* h_W_HOST,
                                       const RealArray_C* h_E_HOST,
                                       RealArray_C* FA_u_W0_HOST,
                                       RealArray_C* FA_u_E0_HOST,
                                       RealArray_C* FA_u_WW_HOST,
                                       RealArray_C* FA_u_EE_HOST,
                                       RealArray_C* uBT_WW_HOST,
                                       RealArray_C* uBT_EE_HOST,
                                       const RealArray_C* du0_HOST,
                                       const RealArray_C* uh_tot_0_HOST,
                                       const RealArray_C* duhdu_tot_0_HOST,
                                       const RealArray_C* du_max_CFL_HOST,
                                       const RealArray_C* du_min_CFL_HOST,
                                       const double dt,
                                       const RealArray_C* dxCu_HOST,
                                       const RealArray_C* dy_Cu_HOST,
                                       const RealArray_C* IareaT_HOST,
                                       const RealArray_C* IdxT_HOST,
                                       const transport_adjust_CS_C* CS_HOST,
                                       const RealArray_C* visc_rem_HOST,
                                       const RealArray_C* visc_rem_max_HOST,
                                       const LogicalArray_C* do_I_HOST,
                                       const RealArray_C* por_face_areaU_HOST)
{
    /// Define Active domain (kernel launch only on real cells)
    amrex::Box bx(amrex::IntVect(bxC_HOST->idxS[0]-1, bxC_HOST->idxS[1]-1, bxC_HOST->idxS[2]-1),
                  amrex::IntVect(bxC_HOST->idxE[0]-1, bxC_HOST->idxE[1]-1, bxC_HOST->idxE[2]-1));

    /// Create A4 containers for the Fortran arrays (2D fields: nz=1)
    auto u_DEV               = turbotmp::make_array4(u_HOST->shape[0], u_HOST->shape[1], u_HOST->shape[2], 1, u_HOST->lb[0], u_HOST->lb[1], u_HOST->lb[2]);
    auto h_in_DEV            = turbotmp::make_array4(h_in_HOST->shape[0], h_in_HOST->shape[1], h_in_HOST->shape[2], 1, h_in_HOST->lb[0], h_in_HOST->lb[1], h_in_HOST->lb[2]);
    auto h_W_DEV             = turbotmp::make_array4(h_W_HOST->shape[0], h_W_HOST->shape[1], h_W_HOST->shape[2], 1, h_W_HOST->lb[0], h_W_HOST->lb[1], h_W_HOST->lb[2]);
    auto h_E_DEV             = turbotmp::make_array4(h_E_HOST->shape[0], h_E_HOST->shape[1], h_E_HOST->shape[2], 1, h_E_HOST->lb[0], h_E_HOST->lb[1], h_E_HOST->lb[2]);
    auto FA_u_W0_DEV         = turbotmp::make_array4(FA_u_W0_HOST->shape[0], FA_u_W0_HOST->shape[1], 1, 1, FA_u_W0_HOST->lb[0], FA_u_W0_HOST->lb[1], 1);
    auto FA_u_E0_DEV         = turbotmp::make_array4(FA_u_E0_HOST->shape[0], FA_u_E0_HOST->shape[1], 1, 1, FA_u_E0_HOST->lb[0], FA_u_E0_HOST->lb[1], 1);
    auto FA_u_WW_DEV         = turbotmp::make_array4(FA_u_WW_HOST->shape[0], FA_u_WW_HOST->shape[1], 1, 1, FA_u_WW_HOST->lb[0], FA_u_WW_HOST->lb[1], 1);
    auto FA_u_EE_DEV         = turbotmp::make_array4(FA_u_EE_HOST->shape[0], FA_u_EE_HOST->shape[1], 1, 1, FA_u_EE_HOST->lb[0], FA_u_EE_HOST->lb[1], 1);
    auto uBT_WW_DEV          = turbotmp::make_array4(uBT_WW_HOST->shape[0], uBT_WW_HOST->shape[1], 1, 1, uBT_WW_HOST->lb[0], uBT_WW_HOST->lb[1], 1);
    auto uBT_EE_DEV          = turbotmp::make_array4(uBT_EE_HOST->shape[0], uBT_EE_HOST->shape[1], 1, 1, uBT_EE_HOST->lb[0], uBT_EE_HOST->lb[1], 1);
    auto du0_DEV             = turbotmp::make_array4(du0_HOST->shape[0], du0_HOST->shape[1], 1, 1, du0_HOST->lb[0], du0_HOST->lb[1], 1);
    auto uh_tot_0_DEV        = turbotmp::make_array4(uh_tot_0_HOST->shape[0], uh_tot_0_HOST->shape[1], 1, 1, uh_tot_0_HOST->lb[0], uh_tot_0_HOST->lb[1], 1);
    auto duhdu_tot_0_DEV     = turbotmp::make_array4(duhdu_tot_0_HOST->shape[0], duhdu_tot_0_HOST->shape[1], 1, 1, duhdu_tot_0_HOST->lb[0], duhdu_tot_0_HOST->lb[1], 1);
    auto du_max_CFL_DEV      = turbotmp::make_array4(du_max_CFL_HOST->shape[0], du_max_CFL_HOST->shape[1], 1, 1, du_max_CFL_HOST->lb[0], du_max_CFL_HOST->lb[1], 1);
    auto du_min_CFL_DEV      = turbotmp::make_array4(du_min_CFL_HOST->shape[0], du_min_CFL_HOST->shape[1], 1, 1, du_min_CFL_HOST->lb[0], du_min_CFL_HOST->lb[1], 1);
    auto dxCu_DEV            = turbotmp::make_array4(dxCu_HOST->shape[0], dxCu_HOST->shape[1], 1, 1, dxCu_HOST->lb[0], dxCu_HOST->lb[1], 1);
    auto dy_Cu_DEV           = turbotmp::make_array4(dy_Cu_HOST->shape[0], dy_Cu_HOST->shape[1], 1, 1, dy_Cu_HOST->lb[0], dy_Cu_HOST->lb[1], 1);
    auto IareaT_DEV          = turbotmp::make_array4(IareaT_HOST->shape[0], IareaT_HOST->shape[1], 1, 1, IareaT_HOST->lb[0], IareaT_HOST->lb[1], 1);
    auto IdxT_DEV            = turbotmp::make_array4(IdxT_HOST->shape[0], IdxT_HOST->shape[1], 1, 1, IdxT_HOST->lb[0], IdxT_HOST->lb[1], 1);
    auto visc_rem_DEV        = turbotmp::make_array4(visc_rem_HOST->shape[0], visc_rem_HOST->shape[1], visc_rem_HOST->shape[2], 1, visc_rem_HOST->lb[0], visc_rem_HOST->lb[1], visc_rem_HOST->lb[2]);
    auto visc_rem_max_DEV    = turbotmp::make_array4(visc_rem_max_HOST->shape[0], visc_rem_max_HOST->shape[1], 1, 1, visc_rem_max_HOST->lb[0], visc_rem_max_HOST->lb[1], 1);
    auto do_I_DEV            = turbotmp::make_int_array4(do_I_HOST->shape[0], do_I_HOST->shape[1], 1, 1, do_I_HOST->lb[0], do_I_HOST->lb[1], 1);
    auto por_face_areaU_DEV  = turbotmp::make_array4(por_face_areaU_HOST->shape[0], por_face_areaU_HOST->shape[1], por_face_areaU_HOST->shape[2], 1, por_face_areaU_HOST->lb[0], por_face_areaU_HOST->lb[1], por_face_areaU_HOST->lb[2]);

    /// Copy host → device (FA_u_*/uBT_* are inout: copy in before kernel)
    turbotmp::copy_FortranHost_to_array4(u_HOST->data,               u_DEV);
    turbotmp::copy_FortranHost_to_array4(h_in_HOST->data,            h_in_DEV);
    turbotmp::copy_FortranHost_to_array4(h_W_HOST->data,             h_W_DEV);
    turbotmp::copy_FortranHost_to_array4(h_E_HOST->data,             h_E_DEV);
    turbotmp::copy_FortranHost_to_array4(FA_u_W0_HOST->data,         FA_u_W0_DEV);
    turbotmp::copy_FortranHost_to_array4(FA_u_E0_HOST->data,         FA_u_E0_DEV);
    turbotmp::copy_FortranHost_to_array4(FA_u_WW_HOST->data,         FA_u_WW_DEV);
    turbotmp::copy_FortranHost_to_array4(FA_u_EE_HOST->data,         FA_u_EE_DEV);
    turbotmp::copy_FortranHost_to_array4(uBT_WW_HOST->data,          uBT_WW_DEV);
    turbotmp::copy_FortranHost_to_array4(uBT_EE_HOST->data,          uBT_EE_DEV);
    turbotmp::copy_FortranHost_to_array4(du0_HOST->data,             du0_DEV);
    turbotmp::copy_FortranHost_to_array4(uh_tot_0_HOST->data,        uh_tot_0_DEV);
    turbotmp::copy_FortranHost_to_array4(duhdu_tot_0_HOST->data,     duhdu_tot_0_DEV);
    turbotmp::copy_FortranHost_to_array4(du_max_CFL_HOST->data,      du_max_CFL_DEV);
    turbotmp::copy_FortranHost_to_array4(du_min_CFL_HOST->data,      du_min_CFL_DEV);
    turbotmp::copy_FortranHost_to_array4(dxCu_HOST->data,            dxCu_DEV);
    turbotmp::copy_FortranHost_to_array4(dy_Cu_HOST->data,           dy_Cu_DEV);
    turbotmp::copy_FortranHost_to_array4(IareaT_HOST->data,          IareaT_DEV);
    turbotmp::copy_FortranHost_to_array4(IdxT_HOST->data,            IdxT_DEV);
    turbotmp::copy_FortranHost_to_array4(visc_rem_HOST->data,        visc_rem_DEV);
    turbotmp::copy_FortranHost_to_array4(visc_rem_max_HOST->data,    visc_rem_max_DEV);
    turbotmp::copy_FortranHost_to_int_array4(do_I_HOST->data,        do_I_DEV);
    turbotmp::copy_FortranHost_to_array4(por_face_areaU_HOST->data,  por_face_areaU_DEV);

    if(verbose) amrex::Print() << "Entered: turbotmp_set_zonal_bt_cont_bridge\n";
    ///-------------------------------------------------
    /// Execute kernel
    MOM::set_zonal_BT_cont(bx,
                           u_DEV.arr,
                           h_in_DEV.arr,
                           h_W_DEV.arr,
                           h_E_DEV.arr,
                           FA_u_W0_DEV.arr,
                           FA_u_E0_DEV.arr,
                           FA_u_WW_DEV.arr,
                           FA_u_EE_DEV.arr,
                           uBT_WW_DEV.arr,
                           uBT_EE_DEV.arr,
                           du0_DEV.arr,
                           uh_tot_0_DEV.arr,
                           duhdu_tot_0_DEV.arr,
                           du_max_CFL_DEV.arr,
                           du_min_CFL_DEV.arr,
                           dt,
                           dxCu_DEV.arr,
                           dy_Cu_DEV.arr,
                           IareaT_DEV.arr,
                           IdxT_DEV.arr,
                           *CS_HOST,
                           visc_rem_DEV.arr,
                           visc_rem_max_DEV.arr,
                           do_I_DEV.arr,
                           por_face_areaU_DEV.arr);

    /// Ensure kernel is done before copying back
    amrex::Gpu::synchronize();

    /// Copy device → host (FA_u_*/uBT_* are the only outputs)
    turbotmp::copy_array4_to_FortranHost(FA_u_W0_DEV, FA_u_W0_HOST->data);
    turbotmp::copy_array4_to_FortranHost(FA_u_E0_DEV, FA_u_E0_HOST->data);
    turbotmp::copy_array4_to_FortranHost(FA_u_WW_DEV, FA_u_WW_HOST->data);
    turbotmp::copy_array4_to_FortranHost(FA_u_EE_DEV, FA_u_EE_HOST->data);
    turbotmp::copy_array4_to_FortranHost(uBT_WW_DEV,  uBT_WW_HOST->data);
    turbotmp::copy_array4_to_FortranHost(uBT_EE_DEV,  uBT_EE_HOST->data);

    /// Free memory from a4 containers
    turbotmp::free_array4(u_DEV);
    turbotmp::free_array4(h_in_DEV);
    turbotmp::free_array4(h_W_DEV);
    turbotmp::free_array4(h_E_DEV);
    turbotmp::free_array4(FA_u_W0_DEV);
    turbotmp::free_array4(FA_u_E0_DEV);
    turbotmp::free_array4(FA_u_WW_DEV);
    turbotmp::free_array4(FA_u_EE_DEV);
    turbotmp::free_array4(uBT_WW_DEV);
    turbotmp::free_array4(uBT_EE_DEV);
    turbotmp::free_array4(du0_DEV);
    turbotmp::free_array4(uh_tot_0_DEV);
    turbotmp::free_array4(duhdu_tot_0_DEV);
    turbotmp::free_array4(du_max_CFL_DEV);
    turbotmp::free_array4(du_min_CFL_DEV);
    turbotmp::free_array4(dxCu_DEV);
    turbotmp::free_array4(dy_Cu_DEV);
    turbotmp::free_array4(IareaT_DEV);
    turbotmp::free_array4(IdxT_DEV);
    turbotmp::free_array4(visc_rem_DEV);
    turbotmp::free_array4(visc_rem_max_DEV);
    turbotmp::free_int_array4(do_I_DEV);
    turbotmp::free_array4(por_face_areaU_DEV);
}

/**
 * @brief Bridge for set_merid_BT_cont
 *
 * @param bxC_HOST            Box over which to iterate
 * @param v_HOST              Meridional velocity (host, Fortran order)
 * @param h_in_HOST           Layer thickness used to calculate fluxes (host, Fortran order)
 * @param h_S_HOST            South edge thickness in the reconstruction (host, Fortran order)
 * @param h_N_HOST            North edge thickness in the reconstruction (host, Fortran order)
 * @param FA_v_S0_HOST        Effective open face area, south, 0 transport (host, Fortran order)
 * @param FA_v_N0_HOST        Effective open face area, north, 0 transport (host, Fortran order)
 * @param FA_v_SS_HOST        Effective open face area, southerly test velocity (host, Fortran order)
 * @param FA_v_NN_HOST        Effective open face area, northerly test velocity (host, Fortran order)
 * @param vBT_SS_HOST         Southerly correction to the barotropic velocity (host, Fortran order)
 * @param vBT_NN_HOST         Northerly correction to the barotropic velocity (host, Fortran order)
 * @param dv0_HOST            Barotropic velocity increment giving 0 transport (host, Fortran order)
 * @param vh_tot_0_HOST       Summed transport with 0 adjustment (host, Fortran order)
 * @param dvhdv_tot_0_HOST    Partial derivative of du_err with dv at 0 adjustment (host, Fortran order)
 * @param dv_max_CFL_HOST     Maximum acceptable value of dv (host, Fortran order)
 * @param dv_min_CFL_HOST     Minimum acceptable value of dv (host, Fortran order)
 * @param dt                  Time increment
 * @param dyCv_HOST           The grid cell's v-point y-extent (host, Fortran order)
 * @param dx_Cv_HOST          Unblocked v-face length, 2D (host, Fortran order)
 * @param IareaT_HOST         1/areaT, 2D (host, Fortran order)
 * @param IdyT_HOST           1/dyT, 2D (host, Fortran order)
 * @param CS_HOST             Transport-adjustment and barotropic-consistency options
 * @param visc_rem_HOST       Fraction of momentum/barotropic acceleration remaining after
 *                            viscosity (host, Fortran order) [nondim]
 * @param visc_rem_max_HOST   Maximum allowable viscosity remnant (host, Fortran order)
 * @param do_I_HOST           Logical flag indicating which I values to work on (host, Fortran order)
 * @param por_face_areaV_HOST Fractional open area of V-faces (host, Fortran order)
 *
 * @note On return, @p FA_v_S0_HOST, @p FA_v_N0_HOST, @p FA_v_SS_HOST, @p FA_v_NN_HOST,
 *       @p vBT_SS_HOST, and @p vBT_NN_HOST hold the modified values.
 */
void turbotmp_set_merid_bt_cont_bridge(const Box_C* bxC_HOST,
                                       const RealArray_C* v_HOST,
                                       const RealArray_C* h_in_HOST,
                                       const RealArray_C* h_S_HOST,
                                       const RealArray_C* h_N_HOST,
                                       RealArray_C* FA_v_S0_HOST,
                                       RealArray_C* FA_v_N0_HOST,
                                       RealArray_C* FA_v_SS_HOST,
                                       RealArray_C* FA_v_NN_HOST,
                                       RealArray_C* vBT_SS_HOST,
                                       RealArray_C* vBT_NN_HOST,
                                       const RealArray_C* dv0_HOST,
                                       const RealArray_C* vh_tot_0_HOST,
                                       const RealArray_C* dvhdv_tot_0_HOST,
                                       const RealArray_C* dv_max_CFL_HOST,
                                       const RealArray_C* dv_min_CFL_HOST,
                                       const double dt,
                                       const RealArray_C* dyCv_HOST,
                                       const RealArray_C* dx_Cv_HOST,
                                       const RealArray_C* IareaT_HOST,
                                       const RealArray_C* IdyT_HOST,
                                       const transport_adjust_CS_C* CS_HOST,
                                       const RealArray_C* visc_rem_HOST,
                                       const RealArray_C* visc_rem_max_HOST,
                                       const LogicalArray_C* do_I_HOST,
                                       const RealArray_C* por_face_areaV_HOST)
{
    /// Define Active domain (kernel launch only on real cells)
    amrex::Box bx(amrex::IntVect(bxC_HOST->idxS[0]-1, bxC_HOST->idxS[1]-1, bxC_HOST->idxS[2]-1),
                  amrex::IntVect(bxC_HOST->idxE[0]-1, bxC_HOST->idxE[1]-1, bxC_HOST->idxE[2]-1));

    /// Create A4 containers for the Fortran arrays (2D fields: nz=1)
    auto v_DEV               = turbotmp::make_array4(v_HOST->shape[0], v_HOST->shape[1], v_HOST->shape[2], 1, v_HOST->lb[0], v_HOST->lb[1], v_HOST->lb[2]);
    auto h_in_DEV            = turbotmp::make_array4(h_in_HOST->shape[0], h_in_HOST->shape[1], h_in_HOST->shape[2], 1, h_in_HOST->lb[0], h_in_HOST->lb[1], h_in_HOST->lb[2]);
    auto h_S_DEV             = turbotmp::make_array4(h_S_HOST->shape[0], h_S_HOST->shape[1], h_S_HOST->shape[2], 1, h_S_HOST->lb[0], h_S_HOST->lb[1], h_S_HOST->lb[2]);
    auto h_N_DEV             = turbotmp::make_array4(h_N_HOST->shape[0], h_N_HOST->shape[1], h_N_HOST->shape[2], 1, h_N_HOST->lb[0], h_N_HOST->lb[1], h_N_HOST->lb[2]);
    auto FA_v_S0_DEV         = turbotmp::make_array4(FA_v_S0_HOST->shape[0], FA_v_S0_HOST->shape[1], 1, 1, FA_v_S0_HOST->lb[0], FA_v_S0_HOST->lb[1], 1);
    auto FA_v_N0_DEV         = turbotmp::make_array4(FA_v_N0_HOST->shape[0], FA_v_N0_HOST->shape[1], 1, 1, FA_v_N0_HOST->lb[0], FA_v_N0_HOST->lb[1], 1);
    auto FA_v_SS_DEV         = turbotmp::make_array4(FA_v_SS_HOST->shape[0], FA_v_SS_HOST->shape[1], 1, 1, FA_v_SS_HOST->lb[0], FA_v_SS_HOST->lb[1], 1);
    auto FA_v_NN_DEV         = turbotmp::make_array4(FA_v_NN_HOST->shape[0], FA_v_NN_HOST->shape[1], 1, 1, FA_v_NN_HOST->lb[0], FA_v_NN_HOST->lb[1], 1);
    auto vBT_SS_DEV          = turbotmp::make_array4(vBT_SS_HOST->shape[0], vBT_SS_HOST->shape[1], 1, 1, vBT_SS_HOST->lb[0], vBT_SS_HOST->lb[1], 1);
    auto vBT_NN_DEV          = turbotmp::make_array4(vBT_NN_HOST->shape[0], vBT_NN_HOST->shape[1], 1, 1, vBT_NN_HOST->lb[0], vBT_NN_HOST->lb[1], 1);
    auto dv0_DEV             = turbotmp::make_array4(dv0_HOST->shape[0], dv0_HOST->shape[1], 1, 1, dv0_HOST->lb[0], dv0_HOST->lb[1], 1);
    auto vh_tot_0_DEV        = turbotmp::make_array4(vh_tot_0_HOST->shape[0], vh_tot_0_HOST->shape[1], 1, 1, vh_tot_0_HOST->lb[0], vh_tot_0_HOST->lb[1], 1);
    auto dvhdv_tot_0_DEV     = turbotmp::make_array4(dvhdv_tot_0_HOST->shape[0], dvhdv_tot_0_HOST->shape[1], 1, 1, dvhdv_tot_0_HOST->lb[0], dvhdv_tot_0_HOST->lb[1], 1);
    auto dv_max_CFL_DEV      = turbotmp::make_array4(dv_max_CFL_HOST->shape[0], dv_max_CFL_HOST->shape[1], 1, 1, dv_max_CFL_HOST->lb[0], dv_max_CFL_HOST->lb[1], 1);
    auto dv_min_CFL_DEV      = turbotmp::make_array4(dv_min_CFL_HOST->shape[0], dv_min_CFL_HOST->shape[1], 1, 1, dv_min_CFL_HOST->lb[0], dv_min_CFL_HOST->lb[1], 1);
    auto dyCv_DEV            = turbotmp::make_array4(dyCv_HOST->shape[0], dyCv_HOST->shape[1], 1, 1, dyCv_HOST->lb[0], dyCv_HOST->lb[1], 1);
    auto dx_Cv_DEV           = turbotmp::make_array4(dx_Cv_HOST->shape[0], dx_Cv_HOST->shape[1], 1, 1, dx_Cv_HOST->lb[0], dx_Cv_HOST->lb[1], 1);
    auto IareaT_DEV          = turbotmp::make_array4(IareaT_HOST->shape[0], IareaT_HOST->shape[1], 1, 1, IareaT_HOST->lb[0], IareaT_HOST->lb[1], 1);
    auto IdyT_DEV            = turbotmp::make_array4(IdyT_HOST->shape[0], IdyT_HOST->shape[1], 1, 1, IdyT_HOST->lb[0], IdyT_HOST->lb[1], 1);
    auto visc_rem_DEV        = turbotmp::make_array4(visc_rem_HOST->shape[0], visc_rem_HOST->shape[1], visc_rem_HOST->shape[2], 1, visc_rem_HOST->lb[0], visc_rem_HOST->lb[1], visc_rem_HOST->lb[2]);
    auto visc_rem_max_DEV    = turbotmp::make_array4(visc_rem_max_HOST->shape[0], visc_rem_max_HOST->shape[1], 1, 1, visc_rem_max_HOST->lb[0], visc_rem_max_HOST->lb[1], 1);
    auto do_I_DEV            = turbotmp::make_int_array4(do_I_HOST->shape[0], do_I_HOST->shape[1], 1, 1, do_I_HOST->lb[0], do_I_HOST->lb[1], 1);
    auto por_face_areaV_DEV  = turbotmp::make_array4(por_face_areaV_HOST->shape[0], por_face_areaV_HOST->shape[1], por_face_areaV_HOST->shape[2], 1, por_face_areaV_HOST->lb[0], por_face_areaV_HOST->lb[1], por_face_areaV_HOST->lb[2]);

    /// Copy host → device (FA_v_*/vBT_* are inout: copy in before kernel)
    turbotmp::copy_FortranHost_to_array4(v_HOST->data,               v_DEV);
    turbotmp::copy_FortranHost_to_array4(h_in_HOST->data,            h_in_DEV);
    turbotmp::copy_FortranHost_to_array4(h_S_HOST->data,             h_S_DEV);
    turbotmp::copy_FortranHost_to_array4(h_N_HOST->data,             h_N_DEV);
    turbotmp::copy_FortranHost_to_array4(FA_v_S0_HOST->data,         FA_v_S0_DEV);
    turbotmp::copy_FortranHost_to_array4(FA_v_N0_HOST->data,         FA_v_N0_DEV);
    turbotmp::copy_FortranHost_to_array4(FA_v_SS_HOST->data,         FA_v_SS_DEV);
    turbotmp::copy_FortranHost_to_array4(FA_v_NN_HOST->data,         FA_v_NN_DEV);
    turbotmp::copy_FortranHost_to_array4(vBT_SS_HOST->data,          vBT_SS_DEV);
    turbotmp::copy_FortranHost_to_array4(vBT_NN_HOST->data,          vBT_NN_DEV);
    turbotmp::copy_FortranHost_to_array4(dv0_HOST->data,             dv0_DEV);
    turbotmp::copy_FortranHost_to_array4(vh_tot_0_HOST->data,        vh_tot_0_DEV);
    turbotmp::copy_FortranHost_to_array4(dvhdv_tot_0_HOST->data,     dvhdv_tot_0_DEV);
    turbotmp::copy_FortranHost_to_array4(dv_max_CFL_HOST->data,      dv_max_CFL_DEV);
    turbotmp::copy_FortranHost_to_array4(dv_min_CFL_HOST->data,      dv_min_CFL_DEV);
    turbotmp::copy_FortranHost_to_array4(dyCv_HOST->data,            dyCv_DEV);
    turbotmp::copy_FortranHost_to_array4(dx_Cv_HOST->data,           dx_Cv_DEV);
    turbotmp::copy_FortranHost_to_array4(IareaT_HOST->data,          IareaT_DEV);
    turbotmp::copy_FortranHost_to_array4(IdyT_HOST->data,            IdyT_DEV);
    turbotmp::copy_FortranHost_to_array4(visc_rem_HOST->data,        visc_rem_DEV);
    turbotmp::copy_FortranHost_to_array4(visc_rem_max_HOST->data,    visc_rem_max_DEV);
    turbotmp::copy_FortranHost_to_int_array4(do_I_HOST->data,        do_I_DEV);
    turbotmp::copy_FortranHost_to_array4(por_face_areaV_HOST->data,  por_face_areaV_DEV);

    if(verbose) amrex::Print() << "Entered: turbotmp_set_merid_bt_cont_bridge\n";
    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::set_merid_BT_cont(bx,
                           v_DEV.arr,
                           h_in_DEV.arr,
                           h_S_DEV.arr,
                           h_N_DEV.arr,
                           FA_v_S0_DEV.arr,
                           FA_v_N0_DEV.arr,
                           FA_v_SS_DEV.arr,
                           FA_v_NN_DEV.arr,
                           vBT_SS_DEV.arr,
                           vBT_NN_DEV.arr,
                           dv0_DEV.arr,
                           vh_tot_0_DEV.arr,
                           dvhdv_tot_0_DEV.arr,
                           dv_max_CFL_DEV.arr,
                           dv_min_CFL_DEV.arr,
                           dt,
                           dyCv_DEV.arr,
                           dx_Cv_DEV.arr,
                           IareaT_DEV.arr,
                           IdyT_DEV.arr,
                           *CS_HOST,
                           visc_rem_DEV.arr,
                           visc_rem_max_DEV.arr,
                           do_I_DEV.arr,
                           por_face_areaV_DEV.arr);

    /// Ensure kernel is done before copying back
    amrex::Gpu::synchronize();

    /// Copy device → host (FA_v_*/vBT_* are the only outputs)
    turbotmp::copy_array4_to_FortranHost(FA_v_S0_DEV, FA_v_S0_HOST->data);
    turbotmp::copy_array4_to_FortranHost(FA_v_N0_DEV, FA_v_N0_HOST->data);
    turbotmp::copy_array4_to_FortranHost(FA_v_SS_DEV, FA_v_SS_HOST->data);
    turbotmp::copy_array4_to_FortranHost(FA_v_NN_DEV, FA_v_NN_HOST->data);
    turbotmp::copy_array4_to_FortranHost(vBT_SS_DEV,  vBT_SS_HOST->data);
    turbotmp::copy_array4_to_FortranHost(vBT_NN_DEV,  vBT_NN_HOST->data);

    /// Free memory from a4 containers
    turbotmp::free_array4(v_DEV);
    turbotmp::free_array4(h_in_DEV);
    turbotmp::free_array4(h_S_DEV);
    turbotmp::free_array4(h_N_DEV);
    turbotmp::free_array4(FA_v_S0_DEV);
    turbotmp::free_array4(FA_v_N0_DEV);
    turbotmp::free_array4(FA_v_SS_DEV);
    turbotmp::free_array4(FA_v_NN_DEV);
    turbotmp::free_array4(vBT_SS_DEV);
    turbotmp::free_array4(vBT_NN_DEV);
    turbotmp::free_array4(dv0_DEV);
    turbotmp::free_array4(vh_tot_0_DEV);
    turbotmp::free_array4(dvhdv_tot_0_DEV);
    turbotmp::free_array4(dv_max_CFL_DEV);
    turbotmp::free_array4(dv_min_CFL_DEV);
    turbotmp::free_array4(dyCv_DEV);
    turbotmp::free_array4(dx_Cv_DEV);
    turbotmp::free_array4(IareaT_DEV);
    turbotmp::free_array4(IdyT_DEV);
    turbotmp::free_array4(visc_rem_DEV);
    turbotmp::free_array4(visc_rem_max_DEV);
    turbotmp::free_int_array4(do_I_DEV);
    turbotmp::free_array4(por_face_areaV_DEV);
}
