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
 * @brief Bridge for zonal_flux_thickness
 *
 * @param bxC_HOST            Box over which to iterate
 * @param u_HOST              Zonal velocity (host, Fortran order)
 * @param h_HOST              Layer thickness used to calculate fluxes (host, Fortran order)
 * @param h_W_HOST            West edge thickness in the reconstruction (host, Fortran order)
 * @param h_E_HOST            East edge thickness in the reconstruction (host, Fortran order)
 * @param h_u_HOST            Effective thickness at zonal faces (host, Fortran order)
 * @param dt                  Time increment
 * @param dy_Cu_HOST          Unblocked u-face length, 2D (host, Fortran order)
 * @param IareaT_HOST         1/areaT, 2D (host, Fortran order)
 * @param IdxT_HOST           1/dxT, 2D (host, Fortran order)
 * @param vol_CFL             If true, rescale face/cell area ratio for CFL
 * @param marginal            If true, report marginal (not transport-averaged) thickness
 * @param obc                 Open boundary control structure
 * @param por_face_areaU_HOST Fractional open area of U-faces (host, Fortran order)
 * @param visc_rem_u_HOST     Both the fraction of the momentum originally in a layer
 *                            that remains after a time-step of viscosity, and the
 *                            fraction of a time-step's worth of a barotropic
 *                            acceleration that a layer experiences after viscosity is
 *                            applied (host, Fortran order) [nondim]. Between 0 (at the
 *                            bottom) and 1 (far above the bottom).
 *
 * @note On return, @p h_u_HOST holds the modified thickness values.
 */
void turbotmp_zonal_flux_thickness_bridge(const Box_C* bxC_HOST,
                                          const RealArray_C* u_HOST,
                                          const RealArray_C* h_HOST,
                                          const RealArray_C* h_W_HOST,
                                          const RealArray_C* h_E_HOST,
                                          RealArray_C* h_u_HOST,
                                          const double dt,
                                          const RealArray_C* dy_Cu_HOST,
                                          const RealArray_C* IareaT_HOST,
                                          const RealArray_C* IdxT_HOST,
                                          const bool vol_CFL,
                                          const bool marginal,
                                          OceanOBC* obc,
                                          const RealArray_C* por_face_areaU_HOST,
                                          const RealArray_C* visc_rem_u_HOST)
{
    /// Define Active domain (kernel launch only on real cells)
    amrex::Box bx(amrex::IntVect(bxC_HOST->idxS[0]-1, bxC_HOST->idxS[1]-1, bxC_HOST->idxS[2]-1),
                  amrex::IntVect(bxC_HOST->idxE[0]-1, bxC_HOST->idxE[1]-1, bxC_HOST->idxE[2]-1));

    /// Create A4 containers for the Fortran arrays (dy_Cu/IareaT/IdxT are 2D: nz=1)
    auto u_DEV      = turbotmp::make_array4(u_HOST->shape[0], u_HOST->shape[1], u_HOST->shape[2], 1,
                                            u_HOST->lb[0], u_HOST->lb[1], u_HOST->lb[2]);
    auto h_DEV      = turbotmp::make_array4(h_HOST->shape[0], h_HOST->shape[1], h_HOST->shape[2], 1,
                                            h_HOST->lb[0], h_HOST->lb[1], h_HOST->lb[2]);
    auto h_W_DEV    = turbotmp::make_array4(h_W_HOST->shape[0], h_W_HOST->shape[1], h_W_HOST->shape[2], 1,
                                            h_W_HOST->lb[0], h_W_HOST->lb[1], h_W_HOST->lb[2]);
    auto h_E_DEV    = turbotmp::make_array4(h_E_HOST->shape[0], h_E_HOST->shape[1], h_E_HOST->shape[2], 1,
                                            h_E_HOST->lb[0], h_E_HOST->lb[1], h_E_HOST->lb[2]);
    auto h_u_DEV    = turbotmp::make_array4(h_u_HOST->shape[0], h_u_HOST->shape[1], h_u_HOST->shape[2], 1,
                                            h_u_HOST->lb[0], h_u_HOST->lb[1], h_u_HOST->lb[2]);
    auto dy_Cu_DEV  = turbotmp::make_array4(dy_Cu_HOST->shape[0], dy_Cu_HOST->shape[1], 1, 1,
                                            dy_Cu_HOST->lb[0], dy_Cu_HOST->lb[1], 1);
    auto IareaT_DEV = turbotmp::make_array4(IareaT_HOST->shape[0], IareaT_HOST->shape[1], 1, 1,
                                            IareaT_HOST->lb[0], IareaT_HOST->lb[1], 1);
    auto IdxT_DEV   = turbotmp::make_array4(IdxT_HOST->shape[0], IdxT_HOST->shape[1], 1, 1,
                                            IdxT_HOST->lb[0], IdxT_HOST->lb[1], 1);
    auto por_face_areaU_DEV = turbotmp::make_array4(por_face_areaU_HOST->shape[0], por_face_areaU_HOST->shape[1],
                                                    por_face_areaU_HOST->shape[2], 1,
                                                    por_face_areaU_HOST->lb[0], por_face_areaU_HOST->lb[1],
                                                    por_face_areaU_HOST->lb[2]);

    /// visc_rem_u_HOST may be absent (data == nullptr); only allocate/copy it when present.
    const bool has_visc_rem_u = (visc_rem_u_HOST->data != nullptr);
    turbotmp::A4Box visc_rem_u_DEV{};
    if (has_visc_rem_u) {
        visc_rem_u_DEV = turbotmp::make_array4(visc_rem_u_HOST->shape[0], visc_rem_u_HOST->shape[1],
                                               visc_rem_u_HOST->shape[2], 1,
                                               visc_rem_u_HOST->lb[0], visc_rem_u_HOST->lb[1],
                                               visc_rem_u_HOST->lb[2]);
    }

    /// Copy host → device (h_u is inout: copy in before kernel)
    turbotmp::copy_FortranHost_to_array4(u_HOST->data, u_DEV);
    turbotmp::copy_FortranHost_to_array4(h_HOST->data, h_DEV);
    turbotmp::copy_FortranHost_to_array4(h_W_HOST->data, h_W_DEV);
    turbotmp::copy_FortranHost_to_array4(h_E_HOST->data, h_E_DEV);
    turbotmp::copy_FortranHost_to_array4(h_u_HOST->data, h_u_DEV);
    turbotmp::copy_FortranHost_to_array4(dy_Cu_HOST->data, dy_Cu_DEV);
    turbotmp::copy_FortranHost_to_array4(IareaT_HOST->data, IareaT_DEV);
    turbotmp::copy_FortranHost_to_array4(IdxT_HOST->data, IdxT_DEV);
    turbotmp::copy_FortranHost_to_array4(por_face_areaU_HOST->data, por_face_areaU_DEV);
    if (has_visc_rem_u) {
        turbotmp::copy_FortranHost_to_array4(visc_rem_u_HOST->data, visc_rem_u_DEV);
    }

    if(verbose) amrex::Print() << "Entered: turbotmp_zonal_flux_thickness_bridge\n";
    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::zonal_flux_thickness(bx, u_DEV.arr, h_DEV.arr, h_W_DEV.arr, h_E_DEV.arr, h_u_DEV.arr, dt,
                              dy_Cu_DEV.arr, IareaT_DEV.arr, IdxT_DEV.arr, vol_CFL, marginal, obc,
                              por_face_areaU_DEV.arr, visc_rem_u_DEV.arr);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device → host (h_u is the only output)
    turbotmp::copy_array4_to_FortranHost(h_u_DEV, h_u_HOST->data);

    /// Free memory from a4 containers (free_array4 on a never-allocated
    /// visc_rem_u_DEV is a safe no-op -- its .data/.data_f are both null)
    turbotmp::free_array4(u_DEV);
    turbotmp::free_array4(h_DEV);
    turbotmp::free_array4(h_W_DEV);
    turbotmp::free_array4(h_E_DEV);
    turbotmp::free_array4(h_u_DEV);
    turbotmp::free_array4(dy_Cu_DEV);
    turbotmp::free_array4(IareaT_DEV);
    turbotmp::free_array4(IdxT_DEV);
    turbotmp::free_array4(por_face_areaU_DEV);
    turbotmp::free_array4(visc_rem_u_DEV);
}

/**
 * @brief Bridge for meridional_flux_thickness
 *
 * @param bxC_HOST            Box over which to iterate
 * @param v_HOST              Meridional velocity (host, Fortran order)
 * @param h_HOST              Layer thickness used to calculate fluxes (host, Fortran order)
 * @param h_S_HOST            South edge thickness in the reconstruction (host, Fortran order)
 * @param h_N_HOST            North edge thickness in the reconstruction (host, Fortran order)
 * @param h_v_HOST            Effective thickness at meridional faces (host, Fortran order)
 * @param dt                  Time increment
 * @param dx_Cv_HOST          Unblocked v-face length, 2D (host, Fortran order)
 * @param IareaT_HOST         1/areaT, 2D (host, Fortran order)
 * @param IdyT_HOST           1/dyT, 2D (host, Fortran order)
 * @param vol_CFL             If true, rescale face/cell area ratio for CFL
 * @param marginal            If true, report marginal (not transport-averaged) thickness
 * @param obc                 Open boundary control structure
 * @param por_face_areaV_HOST Fractional open area of V-faces (host, Fortran order)
 * @param visc_rem_v_HOST     Both the fraction of the momentum originally in a layer
 *                            that remains after a time-step of viscosity, and the
 *                            fraction of a time-step's worth of a barotropic
 *                            acceleration that a layer experiences after viscosity is
 *                            applied (host, Fortran order) [nondim]. Between 0 (at the
 *                            bottom) and 1 (far above the bottom).
 *
 * @note On return, @p h_v_HOST holds the modified thickness values.
 */
void turbotmp_meridional_flux_thickness_bridge(const Box_C* bxC_HOST,
                                               const RealArray_C* v_HOST,
                                               const RealArray_C* h_HOST,
                                               const RealArray_C* h_S_HOST,
                                               const RealArray_C* h_N_HOST,
                                               RealArray_C* h_v_HOST,
                                               const double dt,
                                               const RealArray_C* dx_Cv_HOST,
                                               const RealArray_C* IareaT_HOST,
                                               const RealArray_C* IdyT_HOST,
                                               const bool vol_CFL,
                                               const bool marginal,
                                               OceanOBC* obc,
                                               const RealArray_C* por_face_areaV_HOST,
                                               const RealArray_C* visc_rem_v_HOST)
{
    /// Define Active domain (kernel launch only on real cells)
    amrex::Box bx(amrex::IntVect(bxC_HOST->idxS[0]-1, bxC_HOST->idxS[1]-1, bxC_HOST->idxS[2]-1),
                  amrex::IntVect(bxC_HOST->idxE[0]-1, bxC_HOST->idxE[1]-1, bxC_HOST->idxE[2]-1));

    /// Create A4 containers for the Fortran arrays (dx_Cv/IareaT/IdyT are 2D: nz=1)
    auto v_DEV      = turbotmp::make_array4(v_HOST->shape[0], v_HOST->shape[1], v_HOST->shape[2], 1,
                                            v_HOST->lb[0], v_HOST->lb[1], v_HOST->lb[2]);
    auto h_DEV      = turbotmp::make_array4(h_HOST->shape[0], h_HOST->shape[1], h_HOST->shape[2], 1,
                                            h_HOST->lb[0], h_HOST->lb[1], h_HOST->lb[2]);
    auto h_S_DEV    = turbotmp::make_array4(h_S_HOST->shape[0], h_S_HOST->shape[1], h_S_HOST->shape[2], 1,
                                            h_S_HOST->lb[0], h_S_HOST->lb[1], h_S_HOST->lb[2]);
    auto h_N_DEV    = turbotmp::make_array4(h_N_HOST->shape[0], h_N_HOST->shape[1], h_N_HOST->shape[2], 1,
                                            h_N_HOST->lb[0], h_N_HOST->lb[1], h_N_HOST->lb[2]);
    auto h_v_DEV    = turbotmp::make_array4(h_v_HOST->shape[0], h_v_HOST->shape[1], h_v_HOST->shape[2], 1,
                                            h_v_HOST->lb[0], h_v_HOST->lb[1], h_v_HOST->lb[2]);
    auto dx_Cv_DEV  = turbotmp::make_array4(dx_Cv_HOST->shape[0], dx_Cv_HOST->shape[1], 1, 1,
                                            dx_Cv_HOST->lb[0], dx_Cv_HOST->lb[1], 1);
    auto IareaT_DEV = turbotmp::make_array4(IareaT_HOST->shape[0], IareaT_HOST->shape[1], 1, 1,
                                            IareaT_HOST->lb[0], IareaT_HOST->lb[1], 1);
    auto IdyT_DEV   = turbotmp::make_array4(IdyT_HOST->shape[0], IdyT_HOST->shape[1], 1, 1,
                                            IdyT_HOST->lb[0], IdyT_HOST->lb[1], 1);
    auto por_face_areaV_DEV = turbotmp::make_array4(por_face_areaV_HOST->shape[0], por_face_areaV_HOST->shape[1],
                                                    por_face_areaV_HOST->shape[2], 1,
                                                    por_face_areaV_HOST->lb[0], por_face_areaV_HOST->lb[1],
                                                    por_face_areaV_HOST->lb[2]);

    /// visc_rem_v_HOST may be absent (data == nullptr); only allocate/copy it when present.
    const bool has_visc_rem_v = (visc_rem_v_HOST->data != nullptr);
    turbotmp::A4Box visc_rem_v_DEV{};
    if (has_visc_rem_v) {
        visc_rem_v_DEV = turbotmp::make_array4(visc_rem_v_HOST->shape[0], visc_rem_v_HOST->shape[1],
                                               visc_rem_v_HOST->shape[2], 1,
                                               visc_rem_v_HOST->lb[0], visc_rem_v_HOST->lb[1],
                                               visc_rem_v_HOST->lb[2]);
    }

    /// Copy host → device (h_v is inout: copy in before kernel)
    turbotmp::copy_FortranHost_to_array4(v_HOST->data, v_DEV);
    turbotmp::copy_FortranHost_to_array4(h_HOST->data, h_DEV);
    turbotmp::copy_FortranHost_to_array4(h_S_HOST->data, h_S_DEV);
    turbotmp::copy_FortranHost_to_array4(h_N_HOST->data, h_N_DEV);
    turbotmp::copy_FortranHost_to_array4(h_v_HOST->data, h_v_DEV);
    turbotmp::copy_FortranHost_to_array4(dx_Cv_HOST->data, dx_Cv_DEV);
    turbotmp::copy_FortranHost_to_array4(IareaT_HOST->data, IareaT_DEV);
    turbotmp::copy_FortranHost_to_array4(IdyT_HOST->data, IdyT_DEV);
    turbotmp::copy_FortranHost_to_array4(por_face_areaV_HOST->data, por_face_areaV_DEV);
    if (has_visc_rem_v) {
        turbotmp::copy_FortranHost_to_array4(visc_rem_v_HOST->data, visc_rem_v_DEV);
    }

    if(verbose) amrex::Print() << "Entered: turbotmp_meridional_flux_thickness_bridge\n";
    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::meridional_flux_thickness(bx, v_DEV.arr, h_DEV.arr, h_S_DEV.arr, h_N_DEV.arr, h_v_DEV.arr, dt,
                                   dx_Cv_DEV.arr, IareaT_DEV.arr, IdyT_DEV.arr, vol_CFL, marginal, obc,
                                   por_face_areaV_DEV.arr, visc_rem_v_DEV.arr);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device → host (h_v is the only output)
    turbotmp::copy_array4_to_FortranHost(h_v_DEV, h_v_HOST->data);

    /// Free memory from a4 containers (free_array4 on a never-allocated
    /// visc_rem_v_DEV is a safe no-op -- its .data/.data_f are both null)
    turbotmp::free_array4(v_DEV);
    turbotmp::free_array4(h_DEV);
    turbotmp::free_array4(h_S_DEV);
    turbotmp::free_array4(h_N_DEV);
    turbotmp::free_array4(h_v_DEV);
    turbotmp::free_array4(dx_Cv_DEV);
    turbotmp::free_array4(IareaT_DEV);
    turbotmp::free_array4(IdyT_DEV);
    turbotmp::free_array4(por_face_areaV_DEV);
    turbotmp::free_array4(visc_rem_v_DEV);
}
