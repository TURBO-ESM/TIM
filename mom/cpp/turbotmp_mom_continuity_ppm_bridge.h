#pragma once
/**
 * @file turbotmp_mom_continuity_ppm_bridge.h
 * @brief extern "C" bridge declarations between the MOM6 Fortran
 *        shim and the AMReX PPM continuity kernels (temporary turbotmp layer).
 */

#include "turbotmp_bridge_c_types.h"

struct OceanOBC;    // Undefined at the moment

#ifdef __cplusplus
extern "C" {
#endif

void turbotmp_ppm_limit_pos_bridge(const Box_C* bx_HOST, const RealArray_C* h_in_HOST,
			RealArray_C* h_L_HOST, RealArray_C* h_R_HOST, const double h_min);
void turbotmp_ppm_limit_cw84_bridge(const Box_C* bx_HOST, const RealArray_C* h_in_HOST,
                        RealArray_C* h_L_HOST, RealArray_C* h_R_HOST);
void turbotmp_ppm_reconstruction_y_bridge(const Box_C* bx_HOST, const RealArray_C* h_in_HOST,
                        RealArray_C* h_S_HOST, RealArray_C* h_N_HOST, const RealArray_C* mask2dT_HOST,
                        const double h_min, const bool monotonic, const bool simple_2nd, OceanOBC* obc);
void turbotmp_ppm_reconstruction_x_bridge(const Box_C* bx_HOST, const RealArray_C* h_in_HOST,
                        RealArray_C* h_W_HOST, RealArray_C* h_E_HOST, const RealArray_C* mask2dT_HOST,
                        const double h_min, const bool monotonic, const bool simple_2nd, OceanOBC* obc);
void turbotmp_zonal_edge_thickness_bridge(const Box_C* bx_HOST, const RealArray_C* h_in_HOST,
                        RealArray_C* h_W_HOST, RealArray_C* h_E_HOST, const RealArray_C* mask2dT_HOST,
                        const double h_min, const bool upwind_1st, const bool monotonic,
                        const bool simple_2nd, OceanOBC* obc);
void turbotmp_meridional_edge_thickness_bridge(const Box_C* bx_HOST, const RealArray_C* h_in_HOST,
                        RealArray_C* h_S_HOST, RealArray_C* h_N_HOST, const RealArray_C* mask2dT_HOST,
                        const double h_min, const bool upwind_1st, const bool monotonic,
                        const bool simple_2nd, OceanOBC* obc);

#ifdef __cplusplus
}
#endif
