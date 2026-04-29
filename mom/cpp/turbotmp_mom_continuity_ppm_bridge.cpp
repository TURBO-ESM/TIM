
#include "mom_continuity_ppm.hpp"
#include "turbotmp_helper.hpp"
#include "turbotmp_mom_continuity_ppm_bridge.h"
#include <fstream>
#include <string>

using namespace amrex;

/**
 * @brief Bridge for the function PPM_limit_pos function
 *
 * This function acts as a bridge between a Fortran interface
 * and an AMReX C++ implementation. It also provides the ability
 * to either capture the input, or output or execute the AMReX C++ 
 * implementation based on the setting of the @p mode parameter.
 *
 * @param h_in_HOST Layer thickness [H → m or kg m^-2] 
 * 	on the host in Fortran order
 * @param h_L_HOST, Left thickness of the reconstruction {host, Fortran order} 
 * 	[H → m or kg m^-2]
 * @param h_R_HOST, Right thickness in the reconstruction {host, Fortran order} 
 * 	[H → m or kg m^-2] 
 * @param hmin Minimum thickness allowed by the parabolic fit (host, Fortran order) 
 * 	[H → m or kg m^-2]
 *
 * @return Modified thickness values @p h_L_HOST and @p h_R_HOST
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
    auto h_in_DEV = turbotmp::make_array4(h_in_HOST->shape[0], h_in_HOST->shape[1], h_in_HOST->shape[2], 1);
    auto h_L_DEV  = turbotmp::make_array4(h_L_HOST->shape[0],  h_L_HOST->shape[1],  h_L_HOST->shape[2], 1);
    auto h_R_DEV  = turbotmp::make_array4(h_R_HOST->shape[0],  h_R_HOST->shape[1],  h_R_HOST->shape[2], 1);

    /// Copy from Fortran arrays to A4 container
    turbotmp::copy_FortranHost_to_array4(h_in_HOST->data, h_in_DEV);
    turbotmp::copy_FortranHost_to_array4(h_L_HOST->data, h_L_DEV);
    turbotmp::copy_FortranHost_to_array4(h_R_HOST->data, h_R_DEV);

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
 * @param h_L_HOST, Left thickness of the reconstruction {host, Fortran order}
 *      [H → m or kg m^-2]
 * @param h_R_HOST, Right thickness in the reconstruction {host, Fortran order}
 *      [H → m or kg m^-2]
 *
 * @return Modified thickness values @p h_L_HOST and @p h_R_HOST
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
    auto h_in_DEV = turbotmp::make_array4(h_in_HOST->shape[0], h_in_HOST->shape[1], h_in_HOST->shape[2], 1);
    auto h_L_DEV  = turbotmp::make_array4(h_L_HOST->shape[0],  h_L_HOST->shape[1],  h_L_HOST->shape[2], 1);
    auto h_R_DEV  = turbotmp::make_array4(h_R_HOST->shape[0],  h_R_HOST->shape[1],  h_R_HOST->shape[2], 1);

    /// Copy from Fortran arrays to A4 container
    turbotmp::copy_FortranHost_to_array4(h_in_HOST->data, h_in_DEV);
    turbotmp::copy_FortranHost_to_array4(h_L_HOST->data, h_L_DEV);
    turbotmp::copy_FortranHost_to_array4(h_R_HOST->data, h_R_DEV);

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
 *
 * @return Modified thickness values @p h_S_HOST and @p h_N_HOST
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
    auto h_in_DEV    = turbotmp::make_array4(h_in_HOST->shape[0], h_in_HOST->shape[1], h_in_HOST->shape[2],    1);
    auto h_S_DEV     = turbotmp::make_array4(h_S_HOST->shape[0], h_S_HOST->shape[1],  h_S_HOST->shape[2],     1);
    auto h_N_DEV     = turbotmp::make_array4(h_N_HOST->shape[0], h_N_HOST->shape[1],  h_N_HOST->shape[2],     1);
    auto mask2dT_DEV = turbotmp::make_array4(mask2dT_HOST->shape[0], mask2dT_HOST->shape[1], 1,               1);

    /// Copy from Fortran arrays to A4 container
    turbotmp::copy_FortranHost_to_array4(h_in_HOST->data,    h_in_DEV);
    turbotmp::copy_FortranHost_to_array4(h_S_HOST->data,     h_S_DEV);
    turbotmp::copy_FortranHost_to_array4(h_N_HOST->data,     h_N_DEV);
    turbotmp::copy_FortranHost_to_array4(mask2dT_HOST->data, mask2dT_DEV);

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
