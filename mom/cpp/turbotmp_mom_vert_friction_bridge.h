// SKILLS: 0.3.3
#pragma once
/// @file turbotmp_mom_vert_friction_bridge.h
/// @brief extern "C" bridge declarations between the MOM6 Fortran shims in
///        MOM_vert_friction and the AMReX kernels (temporary turbotmp layer).

#include <stdbool.h>

#include "turbotmp_bridge_c_types.h"

/// Defined in mom_vert_friction.hpp -- field-for-field mirror of the Fortran bind(C) type
/// of the same name.
struct vertvisc_opts_CS_C;
/// @brief Alias letting C callers spell the options struct without the `struct` keyword.
typedef struct vertvisc_opts_CS_C vertvisc_opts_CS_C;
/// Defined in mom_vert_friction.hpp -- field-for-field mirror of the Fortran bind(C) type
/// of the same name.
struct coef_opts_CS_C;
/// @brief Alias letting C callers spell the coupling coefficient options without the `struct`
/// keyword.
typedef struct coef_opts_CS_C coef_opts_CS_C;
/// Opaque: the OBC structure is never dereferenced on this side.
struct OceanOBC;
/// @brief Alias letting C callers spell the opaque OBC handle without the `struct` keyword.
typedef struct OceanOBC OceanOBC;

/// @brief bind(C) mirror of vgrid_consts_CS, field-for-field, same order.
struct vgrid_consts_C {
    int    ke;                     ///< The number of layers/levels in the vertical
    bool   Boussinesq;             ///< If true, make the Boussinesq approximation.
    int    nkml;                   ///< The number of layers at the top that should be treated
                                   ///  as parts of a homogeneous region.
    double dZ_subroundoff;         ///< A thickness in height units that is so small that it can be
                                   ///  added to a vertical distance of Angstrom_Z or 1e-17 m without
                                   ///  changing it at the bit level [Z ~> m].  This is the height
                                   ///  equivalent of H_subroundoff.
    double H_subroundoff;          ///< A thickness that is so small that it can be added to a
                                   ///  thickness of Angstrom or larger without changing it at the bit
                                   ///  level [H ~> m or kg m-2]. If Angstrom is 0 or exceedingly
                                   ///  small, this is negligible compared to 1e-17 m.
    double H_to_m;                 ///< A constant that translates distances in the units of thickness
                                   ///  to m [m H-1 ~> 1 or m3 kg-1].
    double m_to_H;                 ///< A constant that translates distances in m to the units of
                                   ///  thickness [H m-1 ~> 1 or kg m-3].
    double RZ_to_H;                ///< A constant that translates mass per unit area units to
                                   ///  thickness units [H R-1 Z-1 ~> m3 kg-2 or 1].
    double Z_to_H;                 ///< A constant that translates depth units to thickness units
                                   ///  depth [H Z-1 ~> 1 or kg m-3].
    double T_to_s;                 ///< A constant that translates the units of time to seconds
                                   ///  [s T-1 ~> 1]
    double Angstrom_H;             ///< A one-Angstrom thickness in the model thickness units
                                   ///  [H ~> m or kg m-2].
    double H_to_RZ;                ///< A constant that translates thickness units to the units of
                                   ///  mass per unit area [R Z H-1 ~> kg m-3 or 1].
    double m_s_to_L_T;             ///< Convert lateral velocities from m s-1 to L T-1
                                   ///  [L s T-1 m-1 ~> 1]
    struct RealArray_C g_prime_a;  ///< The reduced gravity at each interface [L2 Z-1 T-2 ~> m s-2].
};
/// @brief Alias letting C callers spell the vertical grid constants without the `struct` keyword.
typedef struct vgrid_consts_C vgrid_consts_C;

/// @brief bind(C) mirror of coupling_state_type, field-for-field, same order.
struct coupling_state_C {
    struct RealArray_C a_u;        ///< The u-drag coefficient across an interface
                                   ///  [H T-1 ~> m s-1 or Pa s m-1]
    struct RealArray_C a_u_gl90;   ///< The u-drag coefficient associated with GL90 across an
                                   ///  interface [H T-1 ~> m s-1 or Pa s m-1]
    struct RealArray_C h_u;        ///< The effective layer thickness at u-points [H ~> m or kg m-2].
    struct RealArray_C a_v;        ///< The v-drag coefficient across an interface
                                   ///  [H T-1 ~> m s-1 or Pa s m-1]
    struct RealArray_C a_v_gl90;   ///< The v-drag coefficient associated with GL90 across an
                                   ///  interface [H T-1 ~> m s-1 or Pa s m-1]
    struct RealArray_C h_v;        ///< The effective layer thickness at v-points [H ~> m or kg m-2].
    struct RealArray_C a1_shelf_u; ///< The u-momentum coupling coefficient under ice shelves
                                   ///  [H T-1 ~> m s-1 or Pa s m-1]. Retained to determine stress
                                   ///  under shelves.
    struct RealArray_C a1_shelf_v; ///< The v-momentum coupling coefficient under ice shelves
                                   ///  [H T-1 ~> m s-1 or Pa s m-1]. Retained to determine stress
                                   ///  under shelves.
};
/// @brief Alias letting C callers spell the coupling state mirror without the `struct` keyword.
typedef struct coupling_state_C coupling_state_C;

/// @brief bind(C) mirror of limit_vel_report_type, field-for-field, same order.
struct limit_vel_report_C {
    struct RealArray_C    u_old;        ///< The previous u-velocity [L T-1 ~> m s-1]
    struct RealArray_C    v_old;        ///< The previous v-velocity [L T-1 ~> m s-1]
    struct RealArray_C    vel_report_u; ///< The u-velocity to report [L T-1 ~> m s-1]
    struct RealArray_C    vel_report_v; ///< The v-velocity to report [L T-1 ~> m s-1]
    struct LogicalArray_C dowrite_u;    ///< True for the u-point columns to report
    struct LogicalArray_C dowrite_v;    ///< True for the v-point columns to report
};
/// @brief Alias letting C callers spell the truncation report mirror without the `struct` keyword.
typedef struct limit_vel_report_C limit_vel_report_C;

/// @brief bind(C) mirror of vertvisc_container_type, field-for-field, same order.
struct vertvisc_container_C {
    struct RealArray_C Ray_u;             ///< The Rayleigh drag velocity to be applied to each layer
                                          ///  at u-points [H T-1 ~> m s-1 or Pa s m-1].
    struct RealArray_C Ray_v;             ///< The Rayleigh drag velocity to be applied to each layer
                                          ///  at v-points [H T-1 ~> m s-1 or Pa s m-1].
    struct RealArray_C taux_shelf;        ///< The zonal stresses on the ocean under shelves
                                          ///  [R Z L T-2 ~> Pa].
    struct RealArray_C tauy_shelf;        ///< The meridional stresses on the ocean under shelves
                                          ///  [R Z L T-2 ~> Pa].
    struct RealArray_C kv_bbl_u;          ///< The bottom boundary layer viscosity at the u-points
                                          ///  [H Z T-1 ~> m2 s-1 or Pa s]
    struct RealArray_C kv_bbl_v;          ///< The bottom boundary layer viscosity at the v-points
                                          ///  [H Z T-1 ~> m2 s-1 or Pa s]
    struct RealArray_C bbl_thick_u;       ///< The bottom boundary layer thickness at the u-points
                                          ///  [Z ~> m].
    struct RealArray_C bbl_thick_v;       ///< The bottom boundary layer thickness at the v-points
                                          ///  [Z ~> m].
    struct RealArray_C tbl_thick_shelf_u; ///< Thickness of the viscous top boundary layer under ice
                                          ///  shelves at u-points [Z ~> m].
    struct RealArray_C tbl_thick_shelf_v; ///< Thickness of the viscous top boundary layer under ice
                                          ///  shelves at v-points [Z ~> m].
    struct RealArray_C kv_tbl_shelf_u;    ///< Viscosity in the viscous top boundary layer under ice
                                          ///  shelves at u-points [H Z T-1 ~> m2 s-1 or Pa s]
    struct RealArray_C kv_tbl_shelf_v;    ///< Viscosity in the viscous top boundary layer under ice
                                          ///  shelves at v-points [H Z T-1 ~> m2 s-1 or Pa s]
    struct RealArray_C nkml_visc_u;       ///< The number of layers in the viscous surface mixed
                                          ///  layer at u-points [nondim]. This is not an integer
                                          ///  because there may be fractional layers, and it is
                                          ///  stored in terms of layers, not depth, to facilitate
                                          ///  the movement of the viscous boundary layer with the
                                          ///  flow.
    struct RealArray_C nkml_visc_v;       ///< The number of layers in the viscous surface mixed
                                          ///  layer at v-points [nondim].
    struct RealArray_C Kv_shear;          ///< The shear-driven turbulent vertical viscosity at the
                                          ///  interfaces between layers in tracer columns
                                          ///  [H Z T-1 ~> m2 s-1 or Pa s]
    struct RealArray_C Kv_shear_Bu;       ///< The shear-driven turbulent vertical viscosity at the
                                          ///  interfaces between layers in corner columns
                                          ///  [H Z T-1 ~> m2 s-1 or Pa s]
};
/// @brief Alias letting C callers spell the viscosity container mirror without the `struct` keyword.
typedef struct vertvisc_container_C vertvisc_container_C;

/// @brief bind(C) mirror of mech_forcing_container_type, field-for-field, same order.
struct mech_forcing_container_C {
    struct RealArray_C taux;         ///< zonal wind stress [R L Z T-2 ~> Pa]
    struct RealArray_C tauy;         ///< meridional wind stress [R L Z T-2 ~> Pa]
    struct RealArray_C frac_shelf_u; ///< Fractional ice shelf coverage of u-cells, nondimensional
                                     ///  from 0 to 1 [nondim]. This is only associated if ice
                                     ///  shelves are enabled, and is exactly 0 away from shelves or
                                     ///  on land.
    struct RealArray_C frac_shelf_v; ///< Fractional ice shelf coverage of v-cells, nondimensional
                                     ///  from 0 to 1 [nondim]. This is only associated if ice
                                     ///  shelves are enabled, and is exactly 0 away from shelves or
                                     ///  on land.
};
/// @brief Alias letting C callers spell the mechanical forcing mirror without the `struct` keyword.
typedef struct mech_forcing_container_C mech_forcing_container_C;

/// @brief bind(C) mirror of FPmix_diag_type, field-for-field, same order.
struct FPmix_diag_C {
    struct RealArray_C uE_u;   ///< zonal Eulerian u-pts [L T-1 ~> m s-1]
    struct RealArray_C vE_v;   ///< merid Eulerian v-pts [L T-1 ~> m s-1]
    struct RealArray_C uInc_u; ///< zonal Eulerian u-pts [L T-1 ~> m s-1]
    struct RealArray_C vInc_v; ///< merid Eulerian v-pts [L T-1 ~> m s-1]
};
/// @brief Alias letting C callers spell the FPmix diagnostics mirror without the `struct` keyword.
typedef struct FPmix_diag_C FPmix_diag_C;

/// @brief bind(C) mirror of vertvisc_coef_diag_type, field-for-field, same order.
struct vertvisc_coef_diag_C {
    struct RealArray_C hML_u;     ///< Diagnostic of the mixed layer depth at u points [Z ~> m].
    struct RealArray_C hML_v;     ///< Diagnostic of the mixed layer depth at v points [Z ~> m].
    struct RealArray_C Kv_u;      ///< Total vertical viscosity at u-points in
                                  ///  thickness-based units [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1].
    struct RealArray_C Kv_v;      ///< Total vertical viscosity at v-points in
                                  ///  thickness-based units [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1].
    struct RealArray_C Kv_gl90_u; ///< GL90 vertical viscosity at u-points in
                                  ///  thickness-based units [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1].
    struct RealArray_C Kv_gl90_v; ///< GL90 vertical viscosity at v-points in
                                  ///  thickness-based units [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1].
};
/// @brief Alias letting C callers spell the coupling coefficient diagnostics mirror without the
/// `struct` keyword.
typedef struct vertvisc_coef_diag_C vertvisc_coef_diag_C;

/// @brief bind(C) mirror of vertvisc_AD_type, field-for-field, same order.
struct vertvisc_AD_C {
    struct RealArray_C du_dt_visc;      ///< Zonal acceleration due to vertical viscosity
                                        ///  [L T-2 ~> m s-2]
    struct RealArray_C dv_dt_visc;      ///< Meridional acceleration due to vertical viscosity
                                        ///  [L T-2 ~> m s-2]
    struct RealArray_C du_dt_visc_gl90; ///< Zonal acceleration due to GL90 vertical viscosity
                                        ///  (is included in du_dt_visc) [L T-2 ~> m s-2]
    struct RealArray_C dv_dt_visc_gl90; ///< Meridional acceleration due to GL90 vertical viscosity
                                        ///  (is included in dv_dt_visc) [L T-2 ~> m s-2]
    struct RealArray_C du_dt_str;       ///< Zonal acceleration due to the surface stress (included
                                        ///  in du_dt_visc) [L T-2 ~> m s-2]
    struct RealArray_C dv_dt_str;       ///< Meridional acceleration due to the surface stress
                                        ///  (included in dv_dt_visc) [L T-2 ~> m s-2]
};
/// @brief Alias letting C callers spell the acceleration diagnostics mirror without the `struct`
/// keyword.
typedef struct vertvisc_AD_C vertvisc_AD_C;

#ifdef __cplusplus
extern "C" {
#endif

void turbotmp_vertvisc_limit_vel_bridge(const Box_C* bxU0_HOST,
                                        const Box_C* bxV0_HOST,
                                        RealArray_C* u_HOST,
                                        RealArray_C* v_HOST,
                                        const double dt,
                                        const grid_core_C* Gcore_HOST,
                                        const vgrid_consts_C* vgrid_HOST,
                                        const vertvisc_opts_CS_C* CS_HOST,
                                        const coupling_state_C* coupling_state_HOST,
                                        limit_vel_report_C* report_HOST,
                                        int* ntrunc_sum);

void turbotmp_vertvisc_remnant_tr_bridge(const Box_C* bxU0_HOST,
                                         const Box_C* bxV0_HOST,
                                         const vertvisc_container_C* visc_HOST,
                                         RealArray_C* visc_rem_u_HOST,
                                         RealArray_C* visc_rem_v_HOST,
                                         const double dt,
                                         const grid_core_C* Gcore_HOST,
                                         const int ke,
                                         const coupling_state_C* coupling_state_HOST);

void turbotmp_vertfpmix_tr_bridge(const Box_C* bxU0_HOST,
                                  const Box_C* bxV0_HOST,
                                  RealArray_C* ui_HOST,
                                  RealArray_C* vi_HOST,
                                  const RealArray_C* hbl_h_HOST,
                                  const mech_forcing_container_C* forces_HOST,
                                  const RealArray_C* Us_x_HOST,
                                  const RealArray_C* Us_y_HOST,
                                  const double dt,
                                  const double Cemp_NL,
                                  const double Irho0,
                                  const grid_core_C* Gcore_HOST,
                                  const coupling_state_C* coupling_state_HOST,
                                  FPmix_diag_C* FPdiag_HOST);

void turbotmp_vertvisc_coef_tr_bridge(const Box_C* bxU0_HOST,
                                      const Box_C* bxV0_HOST,
                                      const RealArray_C* u_HOST,
                                      const RealArray_C* v_HOST,
                                      const RealArray_C* h_HOST,
                                      const RealArray_C* dz_HOST,
                                      const mech_forcing_container_C* forces_HOST,
                                      const RealArray_C* Ustar_2d_HOST,
                                      const vertvisc_container_C* visc_HOST,
                                      const RealArray_C* SpV_avg_HOST,
                                      const double dt,
                                      const grid_core_C* Gcore_HOST,
                                      const vgrid_consts_C* vgrid_HOST,
                                      const coef_opts_CS_C* CS_HOST,
                                      coupling_state_C* coupling_state_HOST,
                                      OceanOBC* OBC,
                                      const RealArray_C* kdgl90_struct_HOST,
                                      vertvisc_coef_diag_C* coef_diag_HOST);

void turbotmp_vertvisc_tr_bridge(const Box_C* bxU0_HOST,
                                 const Box_C* bxV0_HOST,
                                 RealArray_C* u_HOST,
                                 RealArray_C* v_HOST,
                                 const RealArray_C* h_HOST,
                                 const mech_forcing_container_C* forces_HOST,
                                 vertvisc_container_C* visc_HOST,
                                 const double dt,
                                 OceanOBC* OBC,
                                 vertvisc_AD_C* AD_HOST,
                                 const grid_core_C* Gcore_HOST,
                                 const vgrid_consts_C* vgrid_HOST,
                                 const vertvisc_opts_CS_C* CS_HOST,
                                 const coupling_state_C* coupling_state_HOST,
                                 RealArray_C* taux_bot_HOST,
                                 RealArray_C* tauy_bot_HOST,
                                 RealArray_C* KE_u_HOST,
                                 RealArray_C* KE_v_HOST,
                                 limit_vel_report_C* report_HOST,
                                 int* ntrunc_sum,
                                 const bool fpmix,
                                 const RealArray_C* Us_x_HOST,
                                 const RealArray_C* Us_y_HOST);

#ifdef __cplusplus
}
#endif
