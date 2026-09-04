// SKILLS: 0.3.1
#pragma once
/**
 * @file turbotmp_mom_continuity_ppm_bridge.h
 * @brief extern "C" bridge declarations between the MOM6 Fortran
 *        shim and the AMReX PPM continuity kernels (temporary turbotmp layer).
 */

#include "turbotmp_bridge_c_types.h"

struct OceanOBC;    // Undefined at the moment
struct transport_adjust_CS_C;   // Defined in mom_continuity_ppm.hpp -- field-for-field
                                 // mirror of the Fortran bind(C) type of the same name.
                                 // A forward declaration suffices here: every prototype
                                 // below only takes a pointer to it.

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
void turbotmp_zonal_bt_mass_flux_bridge(const Box_C* bxC_HOST, const RealArray_C* u_HOST,
                        const RealArray_C* h_in_HOST, const RealArray_C* h_W_HOST, const RealArray_C* h_E_HOST,
                        RealArray_C* uhbt_HOST, const double dt, const RealArray_C* dy_Cu_HOST,
                        const RealArray_C* IareaT_HOST, const RealArray_C* IdxT_HOST,
                        const transport_adjust_CS_C* CS_HOST, OceanOBC* obc,
                        const RealArray_C* por_face_areaU_HOST);
void turbotmp_meridional_bt_mass_flux_bridge(const Box_C* bxC_HOST, const RealArray_C* v_HOST,
                        const RealArray_C* h_in_HOST, const RealArray_C* h_S_HOST, const RealArray_C* h_N_HOST,
                        RealArray_C* vhbt_HOST, const double dt, const RealArray_C* dx_Cv_HOST,
                        const RealArray_C* IareaT_HOST, const RealArray_C* IdyT_HOST,
                        const transport_adjust_CS_C* CS_HOST, OceanOBC* obc,
                        const RealArray_C* por_face_areaV_HOST);

void turbotmp_continuity_zonal_convergence_bridge(const Box_C* bxC_HOST, RealArray_C* h_HOST,
                        const RealArray_C* uh_HOST, const double dt, const RealArray_C* IareaT_HOST,
                        const RealArray_C* hin_HOST, const double h_min);
void turbotmp_continuity_meridional_convergence_bridge(const Box_C* bxC_HOST, RealArray_C* h_HOST,
                        const RealArray_C* vh_HOST, const double dt, const RealArray_C* IareaT_HOST,
                        const RealArray_C* hin_HOST, const double h_min);
void turbotmp_zonal_flux_thickness_bridge(const Box_C* bxC_HOST, const RealArray_C* u_HOST,
                        const RealArray_C* h_HOST, const RealArray_C* h_W_HOST, const RealArray_C* h_E_HOST,
                        RealArray_C* h_u_HOST, const double dt, const RealArray_C* dy_Cu_HOST,
                        const RealArray_C* IareaT_HOST, const RealArray_C* IdxT_HOST,
                        const bool vol_CFL, const bool marginal, OceanOBC* obc,
                        const RealArray_C* por_face_areaU_HOST, const RealArray_C* visc_rem_u_HOST);

void turbotmp_meridional_flux_thickness_bridge(const Box_C* bxC_HOST, const RealArray_C* v_HOST,
                        const RealArray_C* h_HOST, const RealArray_C* h_S_HOST, const RealArray_C* h_N_HOST,
                        RealArray_C* h_v_HOST, const double dt, const RealArray_C* dx_Cv_HOST,
                        const RealArray_C* IareaT_HOST, const RealArray_C* IdyT_HOST,
                        const bool vol_CFL, const bool marginal, OceanOBC* obc,
                        const RealArray_C* por_face_areaV_HOST, const RealArray_C* visc_rem_v_HOST);
void turbotmp_set_zonal_bt_cont_bridge(const Box_C* bxC_HOST, const RealArray_C* u_HOST,
                        const RealArray_C* h_in_HOST, const RealArray_C* h_W_HOST, const RealArray_C* h_E_HOST,
                        RealArray_C* FA_u_W0_HOST, RealArray_C* FA_u_E0_HOST, RealArray_C* FA_u_WW_HOST,
                        RealArray_C* FA_u_EE_HOST, RealArray_C* uBT_WW_HOST, RealArray_C* uBT_EE_HOST,
                        const RealArray_C* du0_HOST, const RealArray_C* uh_tot_0_HOST,
                        const RealArray_C* duhdu_tot_0_HOST, const RealArray_C* du_max_CFL_HOST,
                        const RealArray_C* du_min_CFL_HOST, const double dt, const RealArray_C* dxCu_HOST,
                        const RealArray_C* dy_Cu_HOST, const RealArray_C* IareaT_HOST,
                        const RealArray_C* IdxT_HOST, const transport_adjust_CS_C* CS_HOST,
                        const RealArray_C* visc_rem_HOST, const RealArray_C* visc_rem_max_HOST,
                        const LogicalArray_C* do_I_HOST, const RealArray_C* por_face_areaU_HOST);

void turbotmp_set_merid_bt_cont_bridge(const Box_C* bxC_HOST, const RealArray_C* v_HOST,
                        const RealArray_C* h_in_HOST, const RealArray_C* h_S_HOST, const RealArray_C* h_N_HOST,
                        RealArray_C* FA_v_S0_HOST, RealArray_C* FA_v_N0_HOST, RealArray_C* FA_v_SS_HOST,
                        RealArray_C* FA_v_NN_HOST, RealArray_C* vBT_SS_HOST, RealArray_C* vBT_NN_HOST,
                        const RealArray_C* dv0_HOST, const RealArray_C* vh_tot_0_HOST,
                        const RealArray_C* dvhdv_tot_0_HOST, const RealArray_C* dv_max_CFL_HOST,
                        const RealArray_C* dv_min_CFL_HOST, const double dt, const RealArray_C* dyCv_HOST,
                        const RealArray_C* dx_Cv_HOST, const RealArray_C* IareaT_HOST,
                        const RealArray_C* IdyT_HOST, const transport_adjust_CS_C* CS_HOST,
                        const RealArray_C* visc_rem_HOST, const RealArray_C* visc_rem_max_HOST,
                        const LogicalArray_C* do_I_HOST, const RealArray_C* por_face_areaV_HOST);

void turbotmp_zonal_flux_adjust_bridge(const Box_C* bxC_HOST, const RealArray_C* u_HOST,
                        const RealArray_C* h_in_HOST, const RealArray_C* h_W_HOST, const RealArray_C* h_E_HOST,
                        const RealArray_C* uh_tot_0_HOST, const RealArray_C* duhdu_tot_0_HOST,
                        RealArray_C* du_HOST, const RealArray_C* du_max_CFL_HOST,
                        const RealArray_C* du_min_CFL_HOST, const double dt, const RealArray_C* dy_Cu_HOST,
                        const RealArray_C* IareaT_HOST, const RealArray_C* IdxT_HOST,
                        const transport_adjust_CS_C* CS_HOST, const RealArray_C* visc_rem_HOST,
                        const LogicalArray_C* do_I_in_HOST, const RealArray_C* por_face_areaU_HOST,
                        const RealArray_C* uhbt_HOST, RealArray_C* uh_3d_HOST, OceanOBC* obc);
void turbotmp_meridional_flux_adjust_bridge(const Box_C* bxC_HOST, const RealArray_C* v_HOST,
                        const RealArray_C* h_in_HOST, const RealArray_C* h_S_HOST, const RealArray_C* h_N_HOST,
                        const RealArray_C* vh_tot_0_HOST, const RealArray_C* dvhdv_tot_0_HOST,
                        RealArray_C* dv_HOST, const RealArray_C* dv_max_CFL_HOST,
                        const RealArray_C* dv_min_CFL_HOST, const double dt, const RealArray_C* dx_Cv_HOST,
                        const RealArray_C* IareaT_HOST, const RealArray_C* IdyT_HOST,
                        const transport_adjust_CS_C* CS_HOST, const RealArray_C* visc_rem_HOST,
                        const LogicalArray_C* do_I_in_HOST, const RealArray_C* por_face_areaV_HOST,
                        const RealArray_C* vhbt_HOST, RealArray_C* vh_3d_HOST, OceanOBC* obc);

#ifdef __cplusplus
}
#endif
