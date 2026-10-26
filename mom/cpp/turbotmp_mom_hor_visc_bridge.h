// SKILLS: 0.3.3
#pragma once
/// @file turbotmp_mom_hor_visc_bridge.h
/// @brief extern "C" bridge declarations between the MOM6 Fortran shims in
///        MOM_hor_visc and the AMReX kernels (temporary turbotmp layer).

#include <stdbool.h>

#include "turbotmp_bridge_c_types.h"

/// Defined in mom_hor_visc.hpp -- field-for-field mirror of the Fortran bind(C) type
/// of the same name.
struct hv_opts_CS_C;
/// @brief Alias letting C callers spell the options struct without the `struct` keyword.
typedef struct hv_opts_CS_C hv_opts_CS_C;

/// Opaque: the OBC structure is never dereferenced on this side.
struct OceanOBC;
/// @brief Alias letting C callers spell the opaque OBC handle without the `struct` keyword.
typedef struct OceanOBC OceanOBC;

/// @brief bind(C) mirror of hor_visc_MEKE_type, field-for-field, same order.
struct hor_visc_MEKE_C {
    struct RealArray_C Ku;         ///< The MEKE-derived lateral viscosity coefficient
                                   ///  [L2 T-1 ~> m2 s-1]. This viscosity can be negative when representing
                                   ///  backscatter from unresolved eddies (see Jansen and Held, 2014).
    struct RealArray_C Au;         ///< The MEKE-derived lateral biharmonic viscosity
                                   ///  coefficient [L4 T-1 ~> m4 s-1].
    struct RealArray_C mom_src;    ///< MEKE source from lateral friction in the
                                   ///  momentum equations [R Z L2 T-3 ~> W m-2].
    struct RealArray_C mom_src_bh; ///< MEKE source from the biharmonic part of the lateral friction in the
                                   ///  momentum equations [R Z L2 T-3 ~> W m-2].
    struct RealArray_C GME_snk;    ///< MEKE sink from GME backscatter in the momentum equations
                                   ///  [R Z L2 T-3 ~> W m-2].
    double backscatter_Ro_pow;     ///< Power in Rossby number function for backscatter [nondim].
    double backscatter_Ro_c;       ///< Coefficient in Rossby number function for backscatter [nondim].
};
/// @brief Alias letting C callers spell the MEKE mirror without the `struct` keyword.
typedef struct hor_visc_MEKE_C hor_visc_MEKE_C;

/// @brief bind(C) mirror of hor_visc_VarMix_type, field-for-field, same order.
struct hor_visc_VarMix_C {
    struct RealArray_C Res_fn_h;   ///< Non-dimensional function of the ratio the first baroclinic
                                   ///  deformation radius to the grid spacing at h points [nondim].
    struct RealArray_C Res_fn_q;   ///< Non-dimensional function of the ratio the first baroclinic
                                   ///  deformation radius to the grid spacing at q points [nondim].
    struct RealArray_C BS_struct;  ///< Vertical structure function used in backscatter [nondim]
    bool use_variable_mixing;      ///< If true, use the variable mixing.
    bool Resoln_scaled_Kh;         ///< If true, scale away the Laplacian viscosity
                                   ///  when the deformation radius is well resolved.
};
/// @brief Alias letting C callers spell the variable mixing mirror without the `struct` keyword.
typedef struct hor_visc_VarMix_C hor_visc_VarMix_C;

/// @brief bind(C) mirror of hor_visc_STOCH_type, field-for-field, same order.
struct hor_visc_STOCH_C {
    struct RealArray_C skeb_diss;  ///< Dissipation rate used to set amplitude of SKEBS [L2 T-3 ~> m2 s-3]
                                   ///  Index into this at h points.
    bool skeb_use_frict;           ///< If true, adds viscous dissipation rate to the amplitude of SKEBS
    double skeb_frict_coef;        ///< If skeb_use_frict is true, then skeb_gm_coef * GM_work is added to the
                                   ///  dissipation rate used to set the amplitude of SKEBS [nondim]
};
/// @brief Alias letting C callers spell the STOCH mirror without the `struct` keyword.
typedef struct hor_visc_STOCH_C hor_visc_STOCH_C;

/// @brief bind(C) mirror of hv_metrics_CS, field-for-field, same order.
struct hv_metrics_CS_C {
    struct RealArray_C Kh_bg_xx;         ///< The background Laplacian viscosity at h points [L2 T-1 ~> m2 s-1]. The
                                         ///  actual viscosity may be the larger of this viscosity and the Smagorinsky
                                         ///  and Leith viscosities.
    struct RealArray_C Ah_bg_xx;         ///< The background biharmonic viscosity at h points [L4 T-1 ~> m4 s-1]. The
                                         ///  actual viscosity may be the larger of this viscosity and the Smagorinsky
                                         ///  and Leith viscosities.
    struct RealArray_C reduction_xx;     ///< The amount by which stresses through h points are reduced due to partial
                                         ///  barriers [nondim].
    struct RealArray_C Kh_Max_xx;        ///< The maximum permitted Laplacian viscosity [L2 T-1 ~> m2 s-1].
    struct RealArray_C Ah_Max_xx;        ///< The maximum permitted biharmonic viscosity [L4 T-1 ~> m4 s-1].
    struct RealArray_C Ah_Max_xx_KS;     ///< The maximum permitted biharmonic viscosity for the kill switch
                                         ///  [L4 T-1 ~> m4 s-1].
    struct RealArray_C n1n2_h;           ///< Factor n1*n2 in the anisotropic direction tensor at h-points [nondim]
    struct RealArray_C n1n1_m_n2n2_h;    ///< Factor n1**2-n2**2 in the anisotropic direction tensor at h-points
                                         ///  [nondim]
    struct RealArray_C grid_sp_h2;       ///< Harmonic mean of the squares of the grid [L2 ~> m2]
    struct RealArray_C grid_sp_h3;       ///< Harmonic mean of the squares of the grid^(3/2) [L3 ~> m3]
    struct RealArray_C Kh_bg_xy;         ///< The background Laplacian viscosity at q points [L2 T-1 ~> m2 s-1]. The
                                         ///  actual viscosity may be the larger of this viscosity and the Smagorinsky
                                         ///  and Leith viscosities.
    struct RealArray_C Ah_bg_xy;         ///< The background biharmonic viscosity at q points [L4 T-1 ~> m4 s-1]. The
                                         ///  actual viscosity may be the larger of this viscosity and the Smagorinsky
                                         ///  and Leith viscosities.
    struct RealArray_C reduction_xy;     ///< The amount by which stresses through q points are reduced due to partial
                                         ///  barriers [nondim].
    struct RealArray_C Kh_Max_xy;        ///< The maximum permitted Laplacian viscosity [L2 T-1 ~> m2 s-1].
    struct RealArray_C Ah_Max_xy;        ///< The maximum permitted biharmonic viscosity [L4 T-1 ~> m4 s-1].
    struct RealArray_C Ah_Max_xy_KS;     ///< The maximum permitted biharmonic viscosity for the  kill switch
                                         ///  [L4 T-1 ~> m4 s-1].
    struct RealArray_C n1n2_q;           ///< Factor n1*n2 in the anisotropic direction tensor at q-points [nondim]
    struct RealArray_C n1n1_m_n2n2_q;    ///< Factor n1**2-n2**2 in the anisotropic direction tensor at q-points
                                         ///  [nondim]
    struct RealArray_C dx2h;             ///< Pre-calculated dx^2 at h points [L2 ~> m2]
    struct RealArray_C dy2h;             ///< Pre-calculated dy^2 at h points [L2 ~> m2]
    struct RealArray_C dx_dyT;           ///< Pre-calculated dx/dy at h points [nondim]
    struct RealArray_C dy_dxT;           ///< Pre-calculated dy/dx at h points [nondim]
    struct RealArray_C m_const_leithy;   ///< Pre-calculated .5*sqrt(c_K)*max{dx,dy} [L ~> m]
    struct RealArray_C m_leithy_max;     ///< Pre-calculated 4./max(dx,dy)^2 at h points [L-2 ~> m-2]
    struct RealArray_C dx2q;             ///< Pre-calculated dx^2 at q points [L2 ~> m2]
    struct RealArray_C dy2q;             ///< Pre-calculated dy^2 at q points [L2 ~> m2]
    struct RealArray_C dx_dyBu;          ///< Pre-calculated dx/dy at q points [nondim]
    struct RealArray_C dy_dxBu;          ///< Pre-calculated dy/dx at q points [nondim]
    struct RealArray_C Idx2dyCu;         ///< 1/(dx^2 dy) at u points [L-3 ~> m-3]
    struct RealArray_C Idxdy2u;          ///< 1/(dx dy^2) at u points [L-3 ~> m-3]
    struct RealArray_C Idx2dyCv;         ///< 1/(dx^2 dy) at v points [L-3 ~> m-3]
    struct RealArray_C Idxdy2v;          ///< 1/(dx dy^2) at v points [L-3 ~> m-3]
    struct RealArray_C Laplac2_const_xx; ///< Laplacian metric-dependent constants [L2 ~> m2]
    struct RealArray_C Biharm6_const_xx; ///< Biharmonic metric-dependent constants [L6 ~> m6]
    struct RealArray_C Laplac3_const_xx; ///< Laplacian metric-dependent constants [L3 ~> m3]
    struct RealArray_C Biharm_const_xx;  ///< Biharmonic metric-dependent constants [L4 ~> m4]
    struct RealArray_C Biharm_const2_xx; ///< Biharmonic metric-dependent constants [T L4 ~> s m4]
    struct RealArray_C Re_Ah_const_xx;   ///< Biharmonic metric-dependent constants [L3 ~> m3]
    struct RealArray_C Laplac2_const_xy; ///< Laplacian metric-dependent constants [L2 ~> m2]
    struct RealArray_C Biharm6_const_xy; ///< Biharmonic metric-dependent constants [L6 ~> m6]
    struct RealArray_C Laplac3_const_xy; ///< Laplacian metric-dependent constants [L3 ~> m3]
    struct RealArray_C Biharm_const_xy;  ///< Biharmonic metric-dependent constants [L4 ~> m4]
    struct RealArray_C Biharm_const2_xy; ///< Biharmonic metric-dependent constants [T L4 ~> s m4]
    struct RealArray_C Re_Ah_const_xy;   ///< Biharmonic metric-dependent constants [L3 ~> m3]
};
/// @brief Alias letting C callers spell the metrics mirror without the `struct` keyword.
typedef struct hv_metrics_CS_C hv_metrics_CS_C;

#ifdef __cplusplus
extern "C" {
#endif

/// @brief Bridge for the MOM6 smooth_x9_h subroutine.
///
/// @param Gcore_HOST   Persistent copies of the core ocean grid fields
/// @param bxH0_HOST    The h-point box of the whole column, [isc:iec, jsc:jec, 1:nz]
/// @param field_h_HOST h-point field to be smoothed [arbitrary]
/// @param zero_land    If false, return the average of the surrounding ocean points when
///                     smoothing, otherwise use a value of 0 for land points and include them
///                     in the averages.
void turbotmp_smooth_x9_h_bridge(const grid_core_C* Gcore_HOST,
                                 const Box_C* bxH0_HOST,
                                 RealArray_C* field_h_HOST,
                                 const bool zero_land);

/// @brief Bridge for the MOM6 smooth_x9_uv subroutine.
///
/// @param Gcore_HOST   Persistent copies of the core ocean grid fields
/// @param bxU0_HOST    The u-point box of the whole column, [IscB:IecB, jsc:jec, 1:nz]
/// @param bxV0_HOST    The v-point box of the whole column, [isc:iec, JscB:JecB, 1:nz]
/// @param field_u_HOST u-point field to be smoothed [arbitrary]
/// @param field_v_HOST v-point field to be smoothed [arbitrary]
/// @param zero_land    If false, return the average of the surrounding ocean points when
///                     smoothing, otherwise use a value of 0 for land points and include them
///                     in the averages.
void turbotmp_smooth_x9_uv_bridge(const grid_core_C* Gcore_HOST,
                                  const Box_C* bxU0_HOST,
                                  const Box_C* bxV0_HOST,
                                  RealArray_C* field_u_HOST,
                                  RealArray_C* field_v_HOST,
                                  const bool zero_land);

/// @brief Bridge for the MOM6 smooth_GME_h_pass subroutine.
///
/// @param Gcore_HOST      Persistent copies of the core ocean grid fields
/// @param bxPass_HOST     The h-point range of this smoothing pass, the computational
///                        domain extended by the valid halo points
/// @param GME_flux_h_HOST GME diffusive flux at h points [L2 T-2 ~> m2 s-2]
void turbotmp_smooth_gme_h_pass_bridge(const grid_core_C* Gcore_HOST,
                                       const Box_C* bxPass_HOST,
                                       RealArray_C* GME_flux_h_HOST);

/// @brief Bridge for the MOM6 smooth_GME_q_pass subroutine.
///
/// @param Gcore_HOST      Persistent copies of the core ocean grid fields
/// @param bxPass_HOST     The q-point range of this smoothing pass, the computational
///                        domain extended by the valid halo points
/// @param GME_flux_q_HOST GME diffusive flux at q points [L2 T-2 ~> m2 s-2]
void turbotmp_smooth_gme_q_pass_bridge(const grid_core_C* Gcore_HOST,
                                       const Box_C* bxPass_HOST,
                                       RealArray_C* GME_flux_q_HOST);

/// @brief Bridge for the MOM6 hor_visc_backscatter_h subroutine.
///
/// @param bxQ_HOST               The B-grid-index iteration box of this tile,
///                               [IscB:IecB, JscB:JecB, 1:kmax]
/// @param CS_HOST                Options used by the horizontal viscosity kernels
/// @param MEKE_HOST              MEKE fields related to Mesoscale Eddy Kinetic Energy; each field
///                               may be unallocated (data == nullptr)
/// @param VarMix_HOST            Variable mixing control structure; each field may be unallocated
///                               (data == nullptr)
/// @param use_kh_struct          If true, shape the backscatter coefficient with VarMix%BS_struct
/// @param nkblock                The k-block size used to size the following arrays [nondim]
/// @param kstart                 The first absolute k-layer of the current k-block
/// @param visc_limit_h_flag_HOST determines whether backscatter is shut off [nondim]
/// @param sh_xx_HOST             horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
/// @param str_xx_HOST            The diagonal term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2]
/// @param BS_coeff_h_HOST        A diagnostic array of the backscatter coefficient [L2 T-1 ~> m2 s-1]
void turbotmp_hor_visc_backscatter_h_bridge(const Box_C* bxQ_HOST,
                                            const hv_opts_CS_C* CS_HOST,
                                            const hor_visc_MEKE_C* MEKE_HOST,
                                            const hor_visc_VarMix_C* VarMix_HOST,
                                            const bool use_kh_struct,
                                            const int nkblock,
                                            const int kstart,
                                            const RealArray_C* visc_limit_h_flag_HOST,
                                            const RealArray_C* sh_xx_HOST,
                                            RealArray_C* str_xx_HOST,
                                            RealArray_C* BS_coeff_h_HOST);

/// @brief Bridge for the MOM6 hor_visc_backscatter_q subroutine.
///
/// @param bxH_HOST               The h-point iteration box of this tile, [isc:iec, jsc:jec, 1:kmax]
/// @param CS_HOST                Options used by the horizontal viscosity kernels
/// @param MEKE_HOST              MEKE fields related to Mesoscale Eddy Kinetic Energy; each field
///                               may be unallocated (data == nullptr)
/// @param VarMix_HOST            Variable mixing control structure; each field may be unallocated
///                               (data == nullptr)
/// @param use_kh_struct          If true, shape the backscatter coefficient with VarMix%BS_struct
/// @param nkblock                The k-block size used to size the following arrays [nondim]
/// @param kstart                 The first absolute k-layer of the current k-block
/// @param visc_limit_q_flag_HOST determines whether backscatter is shut off [nondim]
/// @param sh_xy_HOST             horizontal shearing strain (du/dy + dv/dx) including metric terms
///                               [T-1 ~> s-1]
/// @param str_xy_HOST            The cross term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2]
/// @param BS_coeff_q_HOST        A diagnostic array of the backscatter coefficient [L2 T-1 ~> m2 s-1]
void turbotmp_hor_visc_backscatter_q_bridge(const Box_C* bxH_HOST,
                                            const hv_opts_CS_C* CS_HOST,
                                            const hor_visc_MEKE_C* MEKE_HOST,
                                            const hor_visc_VarMix_C* VarMix_HOST,
                                            const bool use_kh_struct,
                                            const int nkblock,
                                            const int kstart,
                                            const RealArray_C* visc_limit_q_flag_HOST,
                                            const RealArray_C* sh_xy_HOST,
                                            RealArray_C* str_xy_HOST,
                                            RealArray_C* BS_coeff_q_HOST);

/// @brief Bridge for the MOM6 hor_visc_GME_bt subroutine.
///
/// @param Gcore_HOST       Persistent copies of the core ocean grid fields
/// @param bxH0_HOST        The h-point box of the whole column, [isc:iec, jsc:jec, 1:nz]
/// @param bxQ0_HOST        The B-grid-index box of the whole column, [IscB:IecB, JscB:JecB, 1:nz]
/// @param CS_HOST          Options used by the horizontal viscosity kernels
/// @param metrics_HOST     Persistent copies of the horizontal viscosity metric and coefficient
///                        arrays; each field may be unallocated (data == nullptr)
/// @param h_HOST           Layer thicknesses [H ~> m or kg m-2]
/// @param ubtav_HOST       zonal barotropic velocity averaged over a baroclinic time-step [L T-1 ~> m s-1]
/// @param vbtav_HOST       meridional barotropic velocity averaged over a baroclinic time-step [L T-1 ~> m s-1]
/// @param dudx_bt_HOST     x-component in the barotropic horizontal tension [T-1 ~> s-1]
/// @param dvdy_bt_HOST     y-component in the barotropic horizontal tension [T-1 ~> s-1]
/// @param dvdx_bt_HOST     x-component in the barotropic shearing strain [T-1 ~> s-1]
/// @param dudy_bt_HOST     y-component in the barotropic shearing strain [T-1 ~> s-1]
/// @param sh_xx_bt_HOST    Barotropic horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
/// @param sh_xy_bt_HOST    Barotropic horizontal shearing strain (du/dy + dv/dx) including metric terms [T-1 ~> s-1]
/// @param GME_effic_h_HOST The filtered efficiency of the GME terms at h points [nondim]
/// @param GME_effic_q_HOST The filtered efficiency of the GME terms at q points [nondim]
void turbotmp_hor_visc_gme_bt_bridge(const grid_core_C* Gcore_HOST,
                                     const Box_C* bxH0_HOST,
                                     const Box_C* bxQ0_HOST,
                                     const hv_opts_CS_C* CS_HOST,
                                     const hv_metrics_CS_C* metrics_HOST,
                                     RealArray_C* h_HOST,
                                     const RealArray_C* ubtav_HOST,
                                     const RealArray_C* vbtav_HOST,
                                     RealArray_C* dudx_bt_HOST,
                                     RealArray_C* dvdy_bt_HOST,
                                     RealArray_C* dvdx_bt_HOST,
                                     RealArray_C* dudy_bt_HOST,
                                     RealArray_C* sh_xx_bt_HOST,
                                     RealArray_C* sh_xy_bt_HOST,
                                     RealArray_C* GME_effic_h_HOST,
                                     RealArray_C* GME_effic_q_HOST);

/// @brief Bridge for the MOM6 hor_visc_Leith_grad_a subroutine.
///
/// @param bxH_HOST                The h-point iteration box of this tile, [isc:iec, jsc:jec, 1:kmax]
/// @param bxKh_HOST               The thickness-point viscosity box of this tile,
///                                [is_Kh:ie_Kh, js_Kh:je_Kh, 1:kmax]
/// @param bxKh_i_HOST             The thickness-point viscosity range in i only, [is_Kh:ie_Kh, jsc:jec, 1:kmax]
/// @param bxKh_j_HOST             The thickness-point viscosity range in j only, [isc:iec, js_Kh:je_Kh, 1:kmax]
/// @param Gcore_HOST              Persistent copies of the core ocean grid fields
/// @param CS_HOST                 Options used by the horizontal viscosity kernels
/// @param nkblock                 The k-block size used to size the following arrays [nondim]
/// @param dudx_HOST               x-derivative of the horizontal tension [T-1 ~> s-1]
/// @param dvdy_HOST               y-derivative of the horizontal tension [T-1 ~> s-1]
/// @param vort_xy_HOST            Vertical vorticity (dv/dx - du/dy) including metric terms [T-1 ~> s-1]
/// @param vort_xy_smooth_HOST     Vertical vorticity including metric terms, smoothed [T-1 ~> s-1]
/// @param grad_vort_mag_h_2d_HOST Magnitude of 2d vorticity gradient at h-points [L-1 T-1 ~> m-1 s-1]
/// @param grad_div_mag_h_HOST     Magnitude of divergence gradient at h-points [L-1 T-1 ~> m-1 s-1]
/// @param grad_vort_mag_q_2d_HOST Magnitude of 2d vorticity gradient at q-points [L-1 T-1 ~> m-1 s-1]
/// @param grad_div_mag_q_HOST     Magnitude of divergence gradient at q-points [L-1 T-1 ~> m-1 s-1]
/// @param Del2vort_q_HOST         Laplacian of vorticity at q-points [L-2 T-1 ~> m-2 s-1]
/// @param vort_xy_dx_HOST         x-derivative of vertical vorticity (d/dx(dv/dx - du/dy)) [L-1 T-1 ~> m-1 s-1]
/// @param vort_xy_dx_smooth_HOST  x-derivative of smoothed vertical vorticity [L-1 T-1 ~> m-1 s-1]
/// @param div_xx_dy_HOST          y-derivative of horizontal divergence (d/dy(du/dx + dv/dy)) [L-1 T-1 ~> m-1 s-1]
/// @param vort_xy_dy_HOST         y-derivative of vertical vorticity (d/dy(dv/dx - du/dy)) [L-1 T-1 ~> m-1 s-1]
/// @param vort_xy_dy_smooth_HOST  y-derivative of smoothed vertical vorticity [L-1 T-1 ~> m-1 s-1]
/// @param div_xx_dx_HOST          x-derivative of horizontal divergence (d/dx(du/dx + dv/dy)) [L-1 T-1 ~> m-1 s-1]
void turbotmp_hor_visc_leith_grad_a_bridge(const Box_C* bxH_HOST,
                                           const Box_C* bxKh_HOST,
                                           const Box_C* bxKh_i_HOST,
                                           const Box_C* bxKh_j_HOST,
                                           const grid_core_C* Gcore_HOST,
                                           const hv_opts_CS_C* CS_HOST,
                                           const int nkblock,
                                           const RealArray_C* dudx_HOST,
                                           const RealArray_C* dvdy_HOST,
                                           const RealArray_C* vort_xy_HOST,
                                           const RealArray_C* vort_xy_smooth_HOST,
                                           RealArray_C* grad_vort_mag_h_2d_HOST,
                                           RealArray_C* grad_div_mag_h_HOST,
                                           RealArray_C* grad_vort_mag_q_2d_HOST,
                                           RealArray_C* grad_div_mag_q_HOST,
                                           RealArray_C* Del2vort_q_HOST,
                                           RealArray_C* vort_xy_dx_HOST,
                                           RealArray_C* vort_xy_dx_smooth_HOST,
                                           RealArray_C* div_xx_dy_HOST,
                                           RealArray_C* vort_xy_dy_HOST,
                                           RealArray_C* vort_xy_dy_smooth_HOST,
                                           RealArray_C* div_xx_dx_HOST);

/// @brief Bridge for the MOM6 hor_visc_Leith_grad_b subroutine.
///
/// @param bxH_HOST                  The h-point iteration box of this tile, [isc:iec, jsc:jec, 1:kmax]
/// @param bxKh_HOST                 The thickness-point viscosity box of this tile,
///                                  [is_Kh:ie_Kh, js_Kh:je_Kh, 1:kmax]
/// @param CS_HOST                   Options used by the horizontal viscosity kernels
/// @param nkblock                   The k-block size used to size the following arrays [nondim]
/// @param grad_vort_mag_h_HOST      Magnitude of vorticity gradient at h-points [L-1 T-1 ~> m-1 s-1]
/// @param grad_vort_mag_q_HOST      Magnitude of vorticity gradient at q-points [L-1 T-1 ~> m-1 s-1]
/// @param vert_vort_mag_smooth_HOST Magnitude of gradient of smoothed vertical vorticity (h or q) [L-1 T-1 ~> m-1 s-1]
/// @param vort_xy_dx_HOST           x-derivative of vertical vorticity (d/dx(dv/dx - du/dy)) [L-1 T-1 ~> m-1 s-1]
/// @param vort_xy_dx_smooth_HOST    x-derivative of smoothed vertical vorticity [L-1 T-1 ~> m-1 s-1]
/// @param vort_xy_dy_HOST           y-derivative of vertical vorticity (d/dy(dv/dx - du/dy)) [L-1 T-1 ~> m-1 s-1]
/// @param vort_xy_dy_smooth_HOST    y-derivative of smoothed vertical vorticity [L-1 T-1 ~> m-1 s-1]
void turbotmp_hor_visc_leith_grad_b_bridge(const Box_C* bxH_HOST,
                                           const Box_C* bxKh_HOST,
                                           const hv_opts_CS_C* CS_HOST,
                                           const int nkblock,
                                           RealArray_C* grad_vort_mag_h_HOST,
                                           RealArray_C* grad_vort_mag_q_HOST,
                                           RealArray_C* vert_vort_mag_smooth_HOST,
                                           RealArray_C* vort_xy_dx_HOST,
                                           RealArray_C* vort_xy_dx_smooth_HOST,
                                           RealArray_C* vort_xy_dy_HOST,
                                           RealArray_C* vort_xy_dy_smooth_HOST);

/// @brief Bridge for the MOM6 hor_visc_Leithy_m subroutine.
///
/// @param bxKh_HOST           The thickness-point viscosity box of this tile, [is_Kh:ie_Kh, js_Kh:je_Kh, 1:kmax]
/// @param Gcore_HOST          Persistent copies of the core ocean grid fields
/// @param CS_HOST             Options used by the horizontal viscosity kernels
/// @param metrics_HOST        Persistent copies of the horizontal viscosity metric and coefficient
///                            arrays; each field may be unallocated (data == nullptr)
/// @param nkblock             The k-block size used to size the following arrays [nondim]
/// @param inv_PI6             The inverse of pi to the sixth power [nondim]
/// @param Del2vort_q_HOST     Laplacian of vorticity at q-points [L-2 T-1 ~> m-2 s-1]
/// @param vert_vort_mag_HOST  Magnitude of the vertical vorticity gradient (h or q) [L-1 T-1 ~> m-1 s-1]
/// @param vort_xy_smooth_HOST Vertical vorticity including metric terms, smoothed [T-1 ~> s-1]
/// @param m_leithy_HOST       Kh=m_leithy*Ah in Leith+E parameterization [L-2 ~> m-2]
void turbotmp_hor_visc_leithy_m_bridge(const Box_C* bxKh_HOST,
                                       const grid_core_C* Gcore_HOST,
                                       const hv_opts_CS_C* CS_HOST,
                                       const hv_metrics_CS_C* metrics_HOST,
                                       const int nkblock,
                                       const double inv_PI6,
                                       const RealArray_C* Del2vort_q_HOST,
                                       const RealArray_C* vert_vort_mag_HOST,
                                       const RealArray_C* vort_xy_smooth_HOST,
                                       RealArray_C* m_leithy_HOST);

/// @brief Bridge for the MOM6 hor_visc_Leithy_Ah subroutine.
///
/// @param bxKh_HOST                 The thickness-point viscosity box of this tile,
///                                  [is_Kh:ie_Kh, js_Kh:je_Kh, 1:kmax]
/// @param CS_HOST                   Options used by the horizontal viscosity kernels
/// @param metrics_HOST              Persistent copies of the horizontal viscosity metric and coefficient
///                                  arrays; each field may be unallocated (data == nullptr)
/// @param nkblock                   The k-block size used to size the following arrays [nondim]
/// @param inv_PI6                   The inverse of pi to the sixth power [nondim]
/// @param Del2vort_q_HOST           Laplacian of vorticity at q-points [L-2 T-1 ~> m-2 s-1]
/// @param vert_vort_mag_smooth_HOST Magnitude of gradient of smoothed vertical vorticity (h or q) [L-1 T-1 ~> m-1 s-1]
/// @param Ah_HOST                   biharmonic viscosity (h or q) [L4 T-1 ~> m4 s-1]
/// @param m_leithy_HOST             Kh=m_leithy*Ah in Leith+E parameterization [L-2 ~> m-2]
/// @param Ah_sq_HOST                The square of the biharmonic viscosity [L8 T-2 ~> m8 s-2]
void turbotmp_hor_visc_leithy_ah_bridge(const Box_C* bxKh_HOST,
                                        const hv_opts_CS_C* CS_HOST,
                                        const hv_metrics_CS_C* metrics_HOST,
                                        const int nkblock,
                                        const double inv_PI6,
                                        const RealArray_C* Del2vort_q_HOST,
                                        const RealArray_C* vert_vort_mag_smooth_HOST,
                                        RealArray_C* Ah_HOST,
                                        const RealArray_C* m_leithy_HOST,
                                        RealArray_C* Ah_sq_HOST);

/// @brief Bridge for the MOM6 hor_visc_Leithy_Ah_h subroutine.
///
/// @param bxKh_HOST    The thickness-point viscosity box of this tile,
///                     [is_Kh:ie_Kh, js_Kh:je_Kh, 1:kmax]
/// @param CS_HOST      Options used by the horizontal viscosity kernels
/// @param metrics_HOST Persistent copies of the horizontal viscosity metric and coefficient
///                     arrays; each field may be unallocated (data == nullptr)
/// @param nkblock      The k-block size used to size the following arrays [nondim]
/// @param kstart       The first absolute k-layer of the current k-block
/// @param Ah_HOST      biharmonic viscosity (h or q) [L4 T-1 ~> m4 s-1]
/// @param Ah_h_HOST    biharmonic viscosity at thickness points [L4 T-1 ~> m4 s-1]
/// @param Ah_sq_HOST   The square of the biharmonic viscosity [L8 T-2 ~> m8 s-2]
void turbotmp_hor_visc_leithy_ah_h_bridge(const Box_C* bxKh_HOST,
                                          const hv_opts_CS_C* CS_HOST,
                                          const hv_metrics_CS_C* metrics_HOST,
                                          const int nkblock,
                                          const int kstart,
                                          RealArray_C* Ah_HOST,
                                          RealArray_C* Ah_h_HOST,
                                          const RealArray_C* Ah_sq_HOST);

/// @brief Bridge for the MOM6 hor_visc_strain subroutine.
///
/// @param Gcore_HOST                 Persistent copies of the core ocean grid fields
/// @param CS_HOST                    Options used by the horizontal viscosity kernels
/// @param metrics_HOST               Persistent copies of the horizontal viscosity metric and coefficient arrays; each
///                                   field may be unallocated (data == nullptr)
/// @param OBC                        Pointer to an open boundary condition type
/// @param bxH_HOST                   The h-point iteration box of this tile, [isc:iec, jsc:jec, 1:kmax]
/// @param bxQ_HOST                   The B-grid-index iteration box of this tile, [IscB:IecB, JscB:JecB, 1:kmax]
/// @param bxU_HOST                   The u-point iteration box of this tile, [IscB:IecB, jsc:jec, 1:kmax]
/// @param bxV_HOST                   The v-point iteration box of this tile, [isc:iec, JscB:JecB, 1:kmax]
/// @param bxKh_HOST                  The thickness-point viscosity box of this tile, [is_Kh:ie_Kh, js_Kh:je_Kh, 1:kmax]
/// @param bxVort_HOST                The vorticity box of this tile, [is_vort:ie_vort, js_vort:je_vort, 1:kmax]
/// @param nkblock                    The effective k-block size [nondim]
/// @param kstart                     The first absolute k-layer of the current k-block
/// @param use_vort_xy                True if vort_xy must be computed
/// @param apply_OBC                  If true, open boundary conditions are applied on this PE
/// @param apply_OBC_strain           If true, the specified open boundary strain configuration is applied
/// @param use_cont_huv               If true, use the continuity solver thicknesses at velocity points
/// @param is_vort                    Loop ranges for vorticity terms
/// @param ie_vort                    Loop ranges for vorticity terms
/// @param js_vort                    Loop ranges for vorticity terms
/// @param je_vort                    Loop ranges for vorticity terms
/// @param u_HOST                     The zonal velocity [L T-1 ~> m s-1].
/// @param v_HOST                     The meridional velocity [L T-1 ~> m s-1].
/// @param h_HOST                     Layer thicknesses [H ~> m or kg m-2].
/// @param hu_cont_HOST               Layer thickness at u-points [H ~> m or kg m-2]; may be unallocated (data ==
///                                   nullptr)
/// @param hv_cont_HOST               Layer thickness at v-points [H ~> m or kg m-2]; may be unallocated (data ==
///                                   nullptr)
/// @param u_smooth_HOST              Zonal velocity, smoothed with a spatial low-pass filter [L T-1 ~> m s-1]
/// @param v_smooth_HOST              Meridional velocity, smoothed with a spatial low-pass filter [L T-1 ~> m s-1]
/// @param NoSt_HOST                  A diagnostic array of normal stress [T-1 ~> s-1].
/// @param ShSt_HOST                  A diagnostic array of shear stress [T-1 ~> s-1].
/// @param Del2u_HOST                 The u-component of the Laplacian of velocity [L-1 T-1 ~> m-1 s-1]
/// @param h_u_HOST                   Thickness interpolated to u points [H ~> m or kg m-2].
/// @param Del2v_HOST                 The v-component of the Laplacian of velocity [L-1 T-1 ~> m-1 s-1]
/// @param h_v_HOST                   Thickness interpolated to v points [H ~> m or kg m-2].
/// @param sh_xx_HOST                 horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
/// @param sh_xx_smooth_HOST          horizontal tension from smoothed velocity including metric terms [T-1 ~> s-1]
/// @param dudx_HOST                  components in the horizontal tension [T-1 ~> s-1]
/// @param dvdy_HOST                  components in the horizontal tension [T-1 ~> s-1]
/// @param dvdx_HOST                  components in the shearing strain [T-1 ~> s-1]
/// @param dudy_HOST                  components in the shearing strain [T-1 ~> s-1]
/// @param dvdx_smooth_HOST           components in the shearing strain from smoothed velocity [T-1 ~> s-1]
/// @param dudy_smooth_HOST           components in the shearing strain from smoothed velocity [T-1 ~> s-1]
/// @param sh_xy_HOST                 horizontal shearing strain (du/dy + dv/dx) including metric terms [T-1 ~> s-1]
/// @param sh_xy_smooth_HOST          horizontal shearing strain from smoothed velocity including metric terms
///                                   [T-1 ~> s-1]
/// @param vort_xy_HOST               Vertical vorticity (dv/dx - du/dy) including metric terms [T-1 ~> s-1]
/// @param vort_xy_smooth_HOST        Vertical vorticity including metric terms, smoothed [T-1 ~> s-1]
void turbotmp_hor_visc_strain_bridge(const grid_core_C* Gcore_HOST,
                                     const hv_opts_CS_C* CS_HOST,
                                     const hv_metrics_CS_C* metrics_HOST,
                                     OceanOBC* OBC,
                                     const Box_C* bxH_HOST,
                                     const Box_C* bxQ_HOST,
                                     const Box_C* bxU_HOST,
                                     const Box_C* bxV_HOST,
                                     const Box_C* bxKh_HOST,
                                     const Box_C* bxVort_HOST,
                                     const int nkblock,
                                     const int kstart,
                                     const bool use_vort_xy,
                                     const bool apply_OBC,
                                     const bool apply_OBC_strain,
                                     const bool use_cont_huv,
                                     const int is_vort,
                                     const int ie_vort,
                                     const int js_vort,
                                     const int je_vort,
                                     const RealArray_C* u_HOST,
                                     const RealArray_C* v_HOST,
                                     const RealArray_C* h_HOST,
                                     RealArray_C* hu_cont_HOST,
                                     RealArray_C* hv_cont_HOST,
                                     const RealArray_C* u_smooth_HOST,
                                     const RealArray_C* v_smooth_HOST,
                                     RealArray_C* NoSt_HOST,
                                     RealArray_C* ShSt_HOST,
                                     RealArray_C* Del2u_HOST,
                                     RealArray_C* h_u_HOST,
                                     RealArray_C* Del2v_HOST,
                                     RealArray_C* h_v_HOST,
                                     RealArray_C* sh_xx_HOST,
                                     RealArray_C* sh_xx_smooth_HOST,
                                     RealArray_C* dudx_HOST,
                                     RealArray_C* dvdy_HOST,
                                     RealArray_C* dvdx_HOST,
                                     RealArray_C* dudy_HOST,
                                     RealArray_C* dvdx_smooth_HOST,
                                     RealArray_C* dudy_smooth_HOST,
                                     RealArray_C* sh_xy_HOST,
                                     RealArray_C* sh_xy_smooth_HOST,
                                     RealArray_C* vort_xy_HOST,
                                     RealArray_C* vort_xy_smooth_HOST);

/// @brief Bridge for the MOM6 hor_visc_h_coeffs_a subroutine.
///
/// @param MEKE_HOST               MEKE fields related to Mesoscale Eddy Kinetic Energy; each field may be unallocated
///                                (data == nullptr)
/// @param VarMix_HOST             Variable mixing control structure; each field may be unallocated (data == nullptr)
/// @param CS_HOST                 Options used by the horizontal viscosity kernels
/// @param metrics_HOST            Persistent copies of the horizontal viscosity metric and coefficient arrays; each
///                                field may be unallocated (data == nullptr)
/// @param bxH_HOST                The h-point iteration box of this tile, [isc:iec, jsc:jec, 1:kmax]
/// @param bxQ_HOST                The B-grid-index iteration box of this tile, [IscB:IecB, JscB:JecB, 1:kmax]
/// @param bxKh_HOST               The thickness-point viscosity box of this tile, [is_Kh:ie_Kh, js_Kh:je_Kh, 1:kmax]
/// @param nkblock                 The effective k-block size [nondim]
/// @param kstart                  The first absolute k-layer of the current k-block
/// @param h_neglect               thickness so small it can be lost in roundoff and so neglected [H ~> m or kg m-2]
/// @param inv_PI3                 Powers of the inverse of pi [nondim]
/// @param inv_PI6                 Powers of the inverse of pi [nondim]
/// @param use_kh_struct           If true, the MEKE viscosity contributions are shaped with VarMix%BS_struct
/// @param use_MEKE_Ku             If true, add the Laplacian viscosity derived from MEKE
/// @param rescale_Kh              If true, the Laplacian viscosity is scaled by the resolution function
/// @param u_HOST                  The zonal velocity [L T-1 ~> m s-1]
/// @param v_HOST                  The meridional velocity [L T-1 ~> m s-1]
/// @param h_HOST                  Layer thicknesses [H ~> m or kg m-2]
/// @param Kh_h_HOST               Laplacian viscosity at thickness points [L2 T-1 ~> m2 s-1]
/// @param grid_Re_Kh_HOST         Grid Reynolds number for Laplacian horizontal viscosity at h points [nondim]
/// @param div_xx_h_HOST           horizontal divergence [T-1 ~> s-1]
/// @param sh_xx_h_HOST            horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
/// @param h_u_HOST                Thickness interpolated to u points [H ~> m or kg m-2]
/// @param h_v_HOST                Thickness interpolated to v points [H ~> m or kg m-2]
/// @param sh_xx_HOST              horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
/// @param str_xx_HOST             str_xx is the diagonal term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2], but
///                                at some points in the code it is not yet layer integrated, so is in
///                                [L2 T-2 ~> m2 s-2]
/// @param grad_vort_mag_h_HOST    Magnitude of vorticity gradient at h-points [L-1 T-1 ~> m-1 s-1]
/// @param grad_vort_mag_h_2d_HOST Magnitude of 2d vorticity gradient at h-points [L-1 T-1 ~> m-1 s-1]
/// @param grad_div_mag_h_HOST     Magnitude of divergence gradient at h-points [L-1 T-1 ~> m-1 s-1]
/// @param dudx_HOST               components in the horizontal tension [T-1 ~> s-1]
/// @param dvdy_HOST               components in the horizontal tension [T-1 ~> s-1]
/// @param sh_xy_HOST              horizontal shearing strain (du/dy + dv/dx) including metric terms [T-1 ~> s-1]
/// @param Del2vort_q_HOST         Laplacian of vorticity at q-points [L-2 T-1 ~> m-2 s-1]
/// @param Ah_HOST                 biharmonic viscosity (h or q) [L4 T-1 ~> m4 s-1]
/// @param Kh_HOST                 Laplacian viscosity (h or q) [L2 T-1 ~> m2 s-1]
/// @param Shear_mag_HOST          magnitude of the shear (h or q) [T-1 ~> s-1]
/// @param vert_vort_mag_HOST      magnitude of the vertical vorticity gradient (h or q) [L-1 T-1 ~> m-1 s-1]
/// @param hrat_min_HOST           h_min divided by the thickness at the stress point (h or q) [nondim]
/// @param visc_bound_rem_HOST     fraction of overall viscous bounds that remain to be applied (h or q) [nondim]
void turbotmp_hor_visc_h_coeffs_a_bridge(const hor_visc_MEKE_C* MEKE_HOST,
                                         const hor_visc_VarMix_C* VarMix_HOST,
                                         const hv_opts_CS_C* CS_HOST,
                                         const hv_metrics_CS_C* metrics_HOST,
                                         const Box_C* bxH_HOST,
                                         const Box_C* bxQ_HOST,
                                         const Box_C* bxKh_HOST,
                                         const int nkblock,
                                         const int kstart,
                                         const double h_neglect,
                                         const double inv_PI3,
                                         const double inv_PI6,
                                         const bool use_kh_struct,
                                         const bool use_MEKE_Ku,
                                         const bool rescale_Kh,
                                         const RealArray_C* u_HOST,
                                         const RealArray_C* v_HOST,
                                         const RealArray_C* h_HOST,
                                         RealArray_C* Kh_h_HOST,
                                         RealArray_C* grid_Re_Kh_HOST,
                                         RealArray_C* div_xx_h_HOST,
                                         RealArray_C* sh_xx_h_HOST,
                                         RealArray_C* h_u_HOST,
                                         RealArray_C* h_v_HOST,
                                         RealArray_C* sh_xx_HOST,
                                         RealArray_C* str_xx_HOST,
                                         RealArray_C* grad_vort_mag_h_HOST,
                                         RealArray_C* grad_vort_mag_h_2d_HOST,
                                         RealArray_C* grad_div_mag_h_HOST,
                                         RealArray_C* dudx_HOST,
                                         RealArray_C* dvdy_HOST,
                                         RealArray_C* sh_xy_HOST,
                                         RealArray_C* Del2vort_q_HOST,
                                         RealArray_C* Ah_HOST,
                                         RealArray_C* Kh_HOST,
                                         RealArray_C* Shear_mag_HOST,
                                         RealArray_C* vert_vort_mag_HOST,
                                         RealArray_C* hrat_min_HOST,
                                         RealArray_C* visc_bound_rem_HOST);

/// @brief Bridge for the MOM6 hor_visc_stress_div subroutine.
///
/// @param MEKE_HOST              MEKE fields related to Mesoscale Eddy Kinetic Energy; each field may be unallocated
///                               (data == nullptr)
/// @param Gcore_HOST             Persistent copies of the core ocean grid fields
/// @param CS_HOST                Options used by the horizontal viscosity kernels
/// @param metrics_HOST           Persistent copies of the horizontal viscosity metric and coefficient arrays; each
///                               field may be unallocated (data == nullptr)
/// @param OBC                    Pointer to an open boundary condition type
/// @param STOCH_HOST             Stochastic control structure; skeb_diss may be unallocated (data == nullptr)
/// @param bxH_HOST               The h-point iteration box of this tile, [isc:iec, jsc:jec, 1:kmax]
/// @param bxQ_HOST               The B-grid-index iteration box of this tile, [IscB:IecB, JscB:JecB, 1:kmax]
/// @param bxU_HOST               The u-point iteration box of this tile, [IscB:IecB, jsc:jec, 1:kmax]
/// @param bxV_HOST               The v-point iteration box of this tile, [isc:iec, JscB:JecB, 1:kmax]
/// @param nkblock                The effective k-block size [nondim]
/// @param kstart                 The first absolute k-layer of the current k-block
/// @param h_neglect              thickness so small it can be lost in roundoff and so neglected [H ~> m or kg m-2]
/// @param apply_OBC              If true, open boundary conditions are applied on this PE
/// @param find_FrictWork         If true, the work done by lateral friction is calculated
/// @param find_FrictWork_bh      If true, the work done by biharmonic lateral friction is calculated
/// @param backscat_subround      The ratio of f over Shear_mag that is so small that the backscatter calculation gives
///                               the same value as if f were 0 [nondim].
/// @param skeb_use_frict         If true, accumulate the frictional dissipation for stochastic backscatter
/// @param u_HOST                 The zonal velocity [L T-1 ~> m s-1]
/// @param v_HOST                 The meridional velocity [L T-1 ~> m s-1]
/// @param h_HOST                 Layer thicknesses [H ~> m or kg m-2]
/// @param uh_HOST                The zonal volume transport [H L2 T-1 ~> m3 s-1]
/// @param vh_HOST                The meridional volume transport [H L2 T-1 ~> m3 s-1]
/// @param diffu_HOST             Zonal acceleration due to convergence of along-coordinate stress tensor
///                               [L T-2 ~> m s-2]
/// @param diffv_HOST             Meridional acceleration due to convergence of along-coordinate stress tensor
///                               [L T-2 ~> m s-2]
/// @param str_xx_GME_HOST        smoothed diagonal term in the stress tensor from GME [L2 T-2 ~> m2 s-2]
/// @param str_xy_GME_HOST        smoothed cross term in the stress tensor from GME [L2 T-2 ~> m2 s-2]
/// @param FrictWork_HOST         work done by MKE dissipation mechanisms [R Z L2 T-3 ~> W m-2]
/// @param FrictWork_bh_HOST      work done by the biharmonic MKE dissipation mechanisms [R Z L2 T-3 ~> W m-2]
/// @param FrictWork_GME_HOST     work done by GME [R Z L2 T-3 ~> W m-2]
/// @param visc_limit_h_flag_HOST determines whether backscatter is shut off [nondim]
/// @param h_u_HOST               Thickness interpolated to u points [H ~> m or kg m-2]
/// @param h_v_HOST               Thickness interpolated to v points [H ~> m or kg m-2]
/// @param sh_xx_HOST             horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
/// @param str_xx_HOST            str_xx is the diagonal term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2], but
///                               at some points in the code it is not yet layer integrated, so is in [L2 T-2 ~> m2 s-2]
/// @param bhstr_xx_HOST          A copy of str_xx that only contains the biharmonic contribution
///                               [H L2 T-2 ~> m3 s-2 or kg s-2]
/// @param sh_xy_HOST             horizontal shearing strain (du/dy + dv/dx) including metric terms [T-1 ~> s-1]
/// @param str_xy_HOST            str_xy is the cross term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2], but at
///                               some points in the code it is not yet layer integrated, so is in [L2 T-2 ~> m2 s-2]
/// @param bhstr_xy_HOST          A copy of str_xy that only contains the biharmonic contribution
///                               [H L2 T-2 ~> m3 s-2 or kg s-2]
/// @param hq_HOST                harmonic mean of the harmonic means of the u- & v point thicknesses [H ~> m or kg m-2]
///                               This form guarantees that hq/hu < 4.
/// @param H_to_RZ                A constant that translates thickness units to the units of mass per unit area
///                               [R Z H-1 ~> kg m-3 or 1]
/// @param s_to_T                 A constant that translates time intervals in seconds to the units of time [T s-1 ~> 1]
void turbotmp_hor_visc_stress_div_bridge(hor_visc_MEKE_C* MEKE_HOST,
                                         const grid_core_C* Gcore_HOST,
                                         const hv_opts_CS_C* CS_HOST,
                                         const hv_metrics_CS_C* metrics_HOST,
                                         OceanOBC* OBC,
                                         hor_visc_STOCH_C* STOCH_HOST,
                                         const Box_C* bxH_HOST,
                                         const Box_C* bxQ_HOST,
                                         const Box_C* bxU_HOST,
                                         const Box_C* bxV_HOST,
                                         const int nkblock,
                                         const int kstart,
                                         const double h_neglect,
                                         const bool apply_OBC,
                                         const bool find_FrictWork,
                                         const bool find_FrictWork_bh,
                                         const double backscat_subround,
                                         const bool skeb_use_frict,
                                         const RealArray_C* u_HOST,
                                         const RealArray_C* v_HOST,
                                         const RealArray_C* h_HOST,
                                         const RealArray_C* uh_HOST,
                                         const RealArray_C* vh_HOST,
                                         RealArray_C* diffu_HOST,
                                         RealArray_C* diffv_HOST,
                                         RealArray_C* str_xx_GME_HOST,
                                         RealArray_C* str_xy_GME_HOST,
                                         RealArray_C* FrictWork_HOST,
                                         RealArray_C* FrictWork_bh_HOST,
                                         RealArray_C* FrictWork_GME_HOST,
                                         RealArray_C* visc_limit_h_flag_HOST,
                                         RealArray_C* h_u_HOST,
                                         RealArray_C* h_v_HOST,
                                         RealArray_C* sh_xx_HOST,
                                         RealArray_C* str_xx_HOST,
                                         RealArray_C* bhstr_xx_HOST,
                                         RealArray_C* sh_xy_HOST,
                                         RealArray_C* str_xy_HOST,
                                         RealArray_C* bhstr_xy_HOST,
                                         RealArray_C* hq_HOST,
                                         const double H_to_RZ,
                                         const double s_to_T);

/// @brief Bridge for the MOM6 hor_visc_h_coeffs_b subroutine.
///
/// @param MEKE_HOST              MEKE fields related to Mesoscale Eddy Kinetic Energy; each field may be unallocated
///                               (data == nullptr)
/// @param VarMix_HOST            Variable mixing control structure; each field may be unallocated (data == nullptr)
/// @param Gcore_HOST             Persistent copies of the core ocean grid fields
/// @param CS_HOST                Options used by the horizontal viscosity kernels
/// @param metrics_HOST           Persistent copies of the horizontal viscosity metric and coefficient arrays; each
///                               field may be unallocated (data == nullptr)
/// @param OBC                    Pointer to an open boundary condition type
/// @param bxH_HOST               The h-point iteration box of this tile, [isc:iec, jsc:jec, 1:kmax]
/// @param bxQ_HOST               The B-grid-index iteration box of this tile, [IscB:IecB, JscB:JecB, 1:kmax]
/// @param bxKh_HOST              The thickness-point viscosity box of this tile, [is_Kh:ie_Kh, js_Kh:je_Kh, 1:kmax]
/// @param nkblock                The effective k-block size [nondim]
/// @param kstart                 The first absolute k-layer of the current k-block
/// @param h_neglect              thickness so small it can be lost in roundoff and so neglected [H ~> m or kg m-2]
/// @param h_neglect3             h_neglect^3 [H3 ~> m3 or kg3 m-6]
/// @param apply_OBC              If true, open boundary conditions are applied on this PE
/// @param use_kh_struct          If true, the MEKE viscosity contributions are shaped with VarMix%BS_struct
/// @param use_MEKE_Au            If true, add the biharmonic viscosity derived from MEKE
/// @param u_HOST                 The zonal velocity [L T-1 ~> m s-1]
/// @param v_HOST                 The meridional velocity [L T-1 ~> m s-1]
/// @param h_HOST                 Layer thicknesses [H ~> m or kg m-2]
/// @param Ah_h_HOST              biharmonic viscosity at thickness points [L4 T-1 ~> m4 s-1]
/// @param Kh_h_HOST              Laplacian viscosity at thickness points [L2 T-1 ~> m2 s-1]
/// @param grid_Re_Ah_HOST        Grid Reynolds number for Biharmonic horizontal viscosity at h points [nondim]
/// @param visc_limit_h_HOST      Used to stabilize the EY24_EBT_BS backscatter [nondim]
/// @param visc_limit_h_flag_HOST determines whether backscatter is shut off [nondim]
/// @param visc_limit_h_frac_HOST determines how close backscatter is to shutting off [nondim]
/// @param BS_coeff_h_HOST        A diagnostic array of the backscatter coefficient [L2 T-1 ~> m2 s-1]
/// @param Del2u_HOST             The u-component of the Laplacian of velocity [L-1 T-1 ~> m-1 s-1]
/// @param h_u_HOST               Thickness interpolated to u points [H ~> m or kg m-2]
/// @param Del2v_HOST             The v-component of the Laplacian of velocity [L-1 T-1 ~> m-1 s-1]
/// @param h_v_HOST               Thickness interpolated to v points [H ~> m or kg m-2]
/// @param sh_xx_HOST             horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
/// @param sh_xx_smooth_HOST      horizontal tension from smoothed velocity including metric terms [T-1 ~> s-1]
/// @param str_xx_HOST            str_xx is the diagonal term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2], but
///                               at some points in the code it is not yet layer integrated, so is in [L2 T-2 ~> m2 s-2]
/// @param bhstr_xx_HOST          A copy of str_xx that only contains the biharmonic contribution
///                               [H L2 T-2 ~> m3 s-2 or kg s-2]
/// @param m_leithy_HOST          Kh=m_leithy*Ah in Leith+E parameterization [L-2 ~> m-2]
/// @param dDel2vdx_HOST          Components in the biharmonic equivalent of the shearing strain [L-2 T-1 ~> m-2 s-1]
/// @param dDel2udy_HOST          Components in the biharmonic equivalent of the shearing strain [L-2 T-1 ~> m-2 s-1]
/// @param sh_xy_HOST             horizontal shearing strain (du/dy + dv/dx) including metric terms [T-1 ~> s-1]
/// @param hq_HOST                harmonic mean of the harmonic means of the u- & v point thicknesses [H ~> m or kg m-2]
///                               This form guarantees that hq/hu < 4.
/// @param Ah_HOST                biharmonic viscosity (h or q) [L4 T-1 ~> m4 s-1]
/// @param Kh_HOST                Laplacian viscosity (h or q) [L2 T-1 ~> m2 s-1]
/// @param Shear_mag_HOST         magnitude of the shear (h or q) [T-1 ~> s-1]
/// @param hrat_min_HOST          h_min divided by the thickness at the stress point (h or q) [nondim]
/// @param visc_bound_rem_HOST    fraction of overall viscous bounds that remain to be applied (h or q) [nondim]
void turbotmp_hor_visc_h_coeffs_b_bridge(const hor_visc_MEKE_C* MEKE_HOST,
                                         const hor_visc_VarMix_C* VarMix_HOST,
                                         const grid_core_C* Gcore_HOST,
                                         const hv_opts_CS_C* CS_HOST,
                                         const hv_metrics_CS_C* metrics_HOST,
                                         OceanOBC* OBC,
                                         const Box_C* bxH_HOST,
                                         const Box_C* bxQ_HOST,
                                         const Box_C* bxKh_HOST,
                                         const int nkblock,
                                         const int kstart,
                                         const double h_neglect,
                                         const double h_neglect3,
                                         const bool apply_OBC,
                                         const bool use_kh_struct,
                                         const bool use_MEKE_Au,
                                         const RealArray_C* u_HOST,
                                         const RealArray_C* v_HOST,
                                         const RealArray_C* h_HOST,
                                         RealArray_C* Ah_h_HOST,
                                         RealArray_C* Kh_h_HOST,
                                         RealArray_C* grid_Re_Ah_HOST,
                                         RealArray_C* visc_limit_h_HOST,
                                         RealArray_C* visc_limit_h_flag_HOST,
                                         RealArray_C* visc_limit_h_frac_HOST,
                                         RealArray_C* BS_coeff_h_HOST,
                                         RealArray_C* Del2u_HOST,
                                         RealArray_C* h_u_HOST,
                                         RealArray_C* Del2v_HOST,
                                         RealArray_C* h_v_HOST,
                                         RealArray_C* sh_xx_HOST,
                                         RealArray_C* sh_xx_smooth_HOST,
                                         RealArray_C* str_xx_HOST,
                                         RealArray_C* bhstr_xx_HOST,
                                         RealArray_C* m_leithy_HOST,
                                         RealArray_C* dDel2vdx_HOST,
                                         RealArray_C* dDel2udy_HOST,
                                         RealArray_C* sh_xy_HOST,
                                         RealArray_C* hq_HOST,
                                         RealArray_C* Ah_HOST,
                                         RealArray_C* Kh_HOST,
                                         RealArray_C* Shear_mag_HOST,
                                         RealArray_C* hrat_min_HOST,
                                         RealArray_C* visc_bound_rem_HOST);

/// @brief Bridge for the MOM6 hor_visc_q_coeffs subroutine.
///
/// @param MEKE_HOST               MEKE fields related to Mesoscale Eddy Kinetic Energy; each field may be unallocated
///                                (data == nullptr)
/// @param VarMix_HOST             Variable mixing control structure; each field may be unallocated (data == nullptr)
/// @param Gcore_HOST              Persistent copies of the core ocean grid fields
/// @param CS_HOST                 Options used by the horizontal viscosity kernels
/// @param metrics_HOST            Persistent copies of the horizontal viscosity metric and coefficient arrays; each
///                                field may be unallocated (data == nullptr)
/// @param bxH_HOST                The h-point iteration box of this tile, [isc:iec, jsc:jec, 1:kmax]
/// @param bxQ_HOST                The B-grid-index iteration box of this tile, [IscB:IecB, JscB:JecB, 1:kmax]
/// @param nkblock                 The effective k-block size [nondim]
/// @param kstart                  The first absolute k-layer of the current k-block
/// @param inv_PI3                 Powers of the inverse of pi [nondim]
/// @param inv_PI6                 Powers of the inverse of pi [nondim]
/// @param find_FrictWork          If true, the work done by lateral friction is calculated
/// @param use_kh_struct           If true, the MEKE viscosity contributions are shaped with VarMix%BS_struct
/// @param use_MEKE_Ku             If true, add the Laplacian viscosity derived from MEKE
/// @param use_MEKE_Au             If true, add the biharmonic viscosity derived from MEKE
/// @param rescale_Kh              If true, the Laplacian viscosity is scaled by the resolution function
/// @param u_HOST                  The zonal velocity [L T-1 ~> m s-1]
/// @param v_HOST                  The meridional velocity [L T-1 ~> m s-1]
/// @param sh_xx_bt_HOST           barotropic horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
/// @param sh_xy_bt_HOST           barotropic horizontal shearing strain (du/dy + dv/dx) inc. metric terms [T-1 ~> s-1]
/// @param GME_effic_h_HOST        The filtered efficiency of the GME terms at h points [nondim]
/// @param GME_effic_q_HOST        The filtered efficiency of the GME terms at q points [nondim]
/// @param KH_u_GME_HOST           Isopycnal height diffusivities in u-columns [L2 T-1 ~> m2 s-1]
/// @param KH_v_GME_HOST           Isopycnal height diffusivities in v-columns [L2 T-1 ~> m2 s-1]
/// @param GME_coeff_h_HOST        GME coefficient at h-points [L2 T-1 ~> m2 s-1]
/// @param GME_coeff_q_HOST        GME coeff. at q-points [L2 T-1 ~> m2 s-1]
/// @param str_xx_GME_HOST         smoothed diagonal term in the stress tensor from GME [L2 T-2 ~> m2 s-2]
/// @param str_xy_GME_HOST         smoothed cross term in the stress tensor from GME [L2 T-2 ~> m2 s-2]
/// @param Ah_h_HOST               biharmonic viscosity at thickness points [L4 T-1 ~> m4 s-1]
/// @param Ah_q_HOST               biharmonic viscosity at corner points [L4 T-1 ~> m4 s-1]
/// @param Kh_h_HOST               Laplacian viscosity at thickness points [L2 T-1 ~> m2 s-1]
/// @param Kh_q_HOST               Laplacian viscosity at corner points [L2 T-1 ~> m2 s-1]
/// @param vort_xy_q_HOST          vertical vorticity at corner points [T-1 ~> s-1]
/// @param sh_xy_q_HOST            horizontal shearing strain at corner points [T-1 ~> s-1]
/// @param visc_limit_q_HOST       used to stabilize the EY24_EBT_BS backscatter [nondim]
/// @param visc_limit_q_flag_HOST  determines whether backscatter is shut off [nondim]
/// @param visc_limit_q_frac_HOST  determines how close backscatter is to shutting off [nondim]
/// @param BS_coeff_q_HOST         A diagnostic array of the backscatter coefficient [L2 T-1 ~> m2 s-1]
/// @param sh_xx_HOST              horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
/// @param str_xx_HOST             str_xx is the diagonal term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2], but
///                                at some points in the code it is not yet layer integrated, so is in
///                                [L2 T-2 ~> m2 s-2]
/// @param dDel2vdx_HOST           Components in the biharmonic equivalent of the shearing strain [L-2 T-1 ~> m-2 s-1]
/// @param dDel2udy_HOST           Components in the biharmonic equivalent of the shearing strain [L-2 T-1 ~> m-2 s-1]
/// @param sh_xy_HOST              horizontal shearing strain (du/dy + dv/dx) including metric terms [T-1 ~> s-1]
/// @param sh_xy_smooth_HOST       horizontal shearing strain from smoothed velocity including metric terms [T-1 ~> s-1]
/// @param str_xy_HOST             str_xy is the cross term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2], but at
///                                some points in the code it is not yet layer integrated, so is in [L2 T-2 ~> m2 s-2]
/// @param bhstr_xy_HOST           A copy of str_xy that only contains the biharmonic contribution
///                                [H L2 T-2 ~> m3 s-2 or kg s-2]
/// @param vort_xy_HOST            Vertical vorticity (dv/dx - du/dy) including metric terms [T-1 ~> s-1]
/// @param grad_vort_mag_q_HOST    Magnitude of vorticity gradient at q-points [L-1 T-1 ~> m-1 s-1]
/// @param grad_vort_mag_q_2d_HOST Magnitude of 2d vorticity gradient at q-points [L-1 T-1 ~> m-1 s-1]
/// @param Del2vort_q_HOST         Laplacian of vorticity at q-points [L-2 T-1 ~> m-2 s-1]
/// @param grad_div_mag_q_HOST     Magnitude of divergence gradient at q-points [L-1 T-1 ~> m-1 s-1]
/// @param hq_HOST                 harmonic mean of the harmonic means of the u- & v point thicknesses
///                                [H ~> m or kg m-2] This form guarantees that hq/hu < 4.
/// @param Ah_HOST                 biharmonic viscosity (h or q) [L4 T-1 ~> m4 s-1]
/// @param Kh_HOST                 Laplacian viscosity (h or q) [L2 T-1 ~> m2 s-1]
/// @param Shear_mag_HOST          magnitude of the shear (h or q) [T-1 ~> s-1]
/// @param vert_vort_mag_HOST      magnitude of the vertical vorticity gradient (h or q) [L-1 T-1 ~> m-1 s-1]
/// @param hrat_min_HOST           h_min divided by the thickness at the stress point (h or q) [nondim]
/// @param visc_bound_rem_HOST     fraction of overall viscous bounds that remain to be applied (h or q) [nondim]
void turbotmp_hor_visc_q_coeffs_bridge(const hor_visc_MEKE_C* MEKE_HOST,
                                       const hor_visc_VarMix_C* VarMix_HOST,
                                       const grid_core_C* Gcore_HOST,
                                       const hv_opts_CS_C* CS_HOST,
                                       const hv_metrics_CS_C* metrics_HOST,
                                       const Box_C* bxH_HOST,
                                       const Box_C* bxQ_HOST,
                                       const int nkblock,
                                       const int kstart,
                                       const double inv_PI3,
                                       const double inv_PI6,
                                       const bool find_FrictWork,
                                       const bool use_kh_struct,
                                       const bool use_MEKE_Ku,
                                       const bool use_MEKE_Au,
                                       const bool rescale_Kh,
                                       const RealArray_C* u_HOST,
                                       const RealArray_C* v_HOST,
                                       const RealArray_C* sh_xx_bt_HOST,
                                       const RealArray_C* sh_xy_bt_HOST,
                                       const RealArray_C* GME_effic_h_HOST,
                                       const RealArray_C* GME_effic_q_HOST,
                                       const RealArray_C* KH_u_GME_HOST,
                                       const RealArray_C* KH_v_GME_HOST,
                                       RealArray_C* GME_coeff_h_HOST,
                                       RealArray_C* GME_coeff_q_HOST,
                                       RealArray_C* str_xx_GME_HOST,
                                       RealArray_C* str_xy_GME_HOST,
                                       RealArray_C* Ah_h_HOST,
                                       RealArray_C* Ah_q_HOST,
                                       RealArray_C* Kh_h_HOST,
                                       RealArray_C* Kh_q_HOST,
                                       RealArray_C* vort_xy_q_HOST,
                                       RealArray_C* sh_xy_q_HOST,
                                       RealArray_C* visc_limit_q_HOST,
                                       RealArray_C* visc_limit_q_flag_HOST,
                                       RealArray_C* visc_limit_q_frac_HOST,
                                       RealArray_C* BS_coeff_q_HOST,
                                       RealArray_C* sh_xx_HOST,
                                       RealArray_C* str_xx_HOST,
                                       RealArray_C* dDel2vdx_HOST,
                                       RealArray_C* dDel2udy_HOST,
                                       RealArray_C* sh_xy_HOST,
                                       RealArray_C* sh_xy_smooth_HOST,
                                       RealArray_C* str_xy_HOST,
                                       RealArray_C* bhstr_xy_HOST,
                                       RealArray_C* vort_xy_HOST,
                                       RealArray_C* grad_vort_mag_q_HOST,
                                       RealArray_C* grad_vort_mag_q_2d_HOST,
                                       RealArray_C* Del2vort_q_HOST,
                                       RealArray_C* grad_div_mag_q_HOST,
                                       RealArray_C* hq_HOST,
                                       RealArray_C* Ah_HOST,
                                       RealArray_C* Kh_HOST,
                                       RealArray_C* Shear_mag_HOST,
                                       RealArray_C* vert_vort_mag_HOST,
                                       RealArray_C* hrat_min_HOST,
                                       RealArray_C* visc_bound_rem_HOST);

/// @brief Bridge for the MOM6 horizontal_viscosity_TR subroutine.
///
/// @param bxH0_HOST              The h-point box of the whole column, [isc:iec, jsc:jec, 1:nz]
/// @param bxQ0_HOST              The B-grid-index box of the whole column, [IscB:IecB, JscB:JecB, 1:nz]
/// @param u_HOST                 The zonal velocity [L T-1 ~> m s-1].
/// @param v_HOST                 The meridional velocity [L T-1 ~> m s-1].
/// @param h_HOST                 Layer thicknesses [H ~> m or kg m-2].
/// @param uh_HOST                The zonal volume transport [H L2 T-1 ~> m3 s-1].
/// @param vh_HOST                The meridional volume transport [H L2 T-1 ~> m3 s-1].
/// @param diffu_HOST             Zonal acceleration due to convergence of along-coordinate stress tensor
///                               [L T-2 ~> m s-2]
/// @param diffv_HOST             Meridional acceleration due to convergence of along-coordinate stress tensor
///                               [L T-2 ~> m s-2].
/// @param MEKE_HOST              MEKE fields related to Mesoscale Eddy Kinetic Energy.; each field may be unallocated
///                               (data == nullptr)
/// @param VarMix_HOST            Variable mixing control structure; each field may be unallocated (data == nullptr)
/// @param Gcore_HOST             Persistent copies of the core ocean grid fields
/// @param CS_HOST                Options used by the horizontal viscosity kernels
/// @param metrics_HOST           Persistent copies of the horizontal viscosity metric and coefficient arrays; each
///                               field may be unallocated (data == nullptr)
/// @param nkblock                The effective k-block size [nondim]
/// @param OBC                    Pointer to an open boundary condition type
/// @param hu_cont_HOST           Layer thickness at u-points [H ~> m or kg m-2].
/// @param hv_cont_HOST           Layer thickness at v-points [H ~> m or kg m-2].
/// @param STOCH_HOST             Stochastic control structure; each field may be unallocated (data == nullptr)
/// @param u_smooth_HOST          Zonal velocity, smoothed with a spatial low-pass filter [L T-1 ~> m s-1]
/// @param v_smooth_HOST          Meridional velocity, smoothed with a spatial low-pass filter [L T-1 ~> m s-1]
/// @param dz_HOST                Height change across layers [Z ~> m]
/// @param slope_x_HOST           Isopycnal slope in i-direction [Z L-1 ~> nondim]
/// @param slope_y_HOST           Isopycnal slope in j-direction [Z L-1 ~> nondim]
/// @param sh_xx_bt_HOST          barotropic horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
/// @param sh_xy_bt_HOST          barotropic horizontal shearing strain (du/dy + dv/dx) inc. metric terms [T-1 ~> s-1]
/// @param GME_effic_h_HOST       The filtered efficiency of the GME terms at h points [nondim]
/// @param GME_effic_q_HOST       The filtered efficiency of the GME terms at q points [nondim]
/// @param KH_u_GME_HOST          Isopycnal height diffusivities in u-columns [L2 T-1 ~> m2 s-1]
/// @param KH_v_GME_HOST          Isopycnal height diffusivities in v-columns [L2 T-1 ~> m2 s-1]
/// @param GME_coeff_h_HOST       GME coefficient at h-points [L2 T-1 ~> m2 s-1]
/// @param GME_coeff_q_HOST       GME coeff. at q-points [L2 T-1 ~> m2 s-1]
/// @param str_xx_GME_HOST        smoothed diagonal term in the stress tensor from GME [L2 T-2 ~> m2 s-2]
/// @param str_xy_GME_HOST        smoothed cross term in the stress tensor from GME [L2 T-2 ~> m2 s-2]
/// @param NoSt_HOST              A diagnostic array of normal stress [T-1 ~> s-1].
/// @param ShSt_HOST              A diagnostic array of shear stress [T-1 ~> s-1].
/// @param FrictWork_HOST         work done by MKE dissipation mechanisms [R Z L2 T-3 ~> W m-2]
/// @param FrictWork_bh_HOST      work done by the biharmonic MKE dissipation mechanisms [R Z L2 T-3 ~> W m-2]
/// @param FrictWork_GME_HOST     work done by GME [R Z L2 T-3 ~> W m-2]
/// @param Ah_h_HOST              biharmonic viscosity at thickness points [L4 T-1 ~> m4 s-1]
/// @param Ah_q_HOST              biharmonic viscosity at corner points [L4 T-1 ~> m4 s-1]
/// @param Kh_h_HOST              Laplacian viscosity at thickness points [L2 T-1 ~> m2 s-1]
/// @param Kh_q_HOST              Laplacian viscosity at corner points [L2 T-1 ~> m2 s-1]
/// @param grid_Re_Kh_HOST        Grid Reynolds number for Laplacian horizontal viscosity at h points [nondim]
/// @param grid_Re_Ah_HOST        Grid Reynolds number for Biharmonic horizontal viscosity at h points [nondim]
/// @param div_xx_h_HOST          horizontal divergence [T-1 ~> s-1]
/// @param sh_xx_h_HOST           horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
/// @param vort_xy_q_HOST         vertical vorticity at corner points [T-1 ~> s-1]
/// @param sh_xy_q_HOST           horizontal shearing strain at corner points [T-1 ~> s-1]
/// @param visc_limit_h_HOST      Used to stabilize the EY24_EBT_BS backscatter [nondim]
/// @param visc_limit_q_HOST      used to stabilize the EY24_EBT_BS backscatter [nondim]
/// @param visc_limit_h_flag_HOST determines whether backscatter is shut off [nondim]
/// @param visc_limit_q_flag_HOST determines whether backscatter is shut off [nondim]
/// @param visc_limit_h_frac_HOST determines how close backscatter is to shutting off [nondim]
/// @param visc_limit_q_frac_HOST determines how close backscatter is to shutting off [nondim]
/// @param BS_coeff_h_HOST        A diagnostic array of the backscatter coefficient [L2 T-1 ~> m2 s-1]
/// @param BS_coeff_q_HOST        A diagnostic array of the backscatter coefficient [L2 T-1 ~> m2 s-1]
/// @param H_subroundoff          A thickness that is so small that it can be added to a thickness of Angstrom or larger
///                               without changing it at the bit level [H ~> m or kg m-2]
/// @param H_to_RZ                A constant that translates thickness units to the units of mass per unit area
///                               [R Z H-1 ~> kg m-3 or 1]
/// @param s_to_T                 A constant that translates time intervals in seconds to the units of time [T s-1 ~> 1]
/// @param host                   Opaque pointer to the hv_host_ctx that holds the Fortran objects used by the host
///                               helpers
void turbotmp_horizontal_viscosity_tr_bridge(const Box_C* bxH0_HOST,
                                             const Box_C* bxQ0_HOST,
                                             const RealArray_C* u_HOST,
                                             const RealArray_C* v_HOST,
                                             const RealArray_C* h_HOST,
                                             const RealArray_C* uh_HOST,
                                             const RealArray_C* vh_HOST,
                                             RealArray_C* diffu_HOST,
                                             RealArray_C* diffv_HOST,
                                             hor_visc_MEKE_C* MEKE_HOST,
                                             const hor_visc_VarMix_C* VarMix_HOST,
                                             const grid_core_C* Gcore_HOST,
                                             const hv_opts_CS_C* CS_HOST,
                                             const hv_metrics_CS_C* metrics_HOST,
                                             const int nkblock,
                                             OceanOBC* OBC,
                                             RealArray_C* hu_cont_HOST,
                                             RealArray_C* hv_cont_HOST,
                                             hor_visc_STOCH_C* STOCH_HOST,
                                             const RealArray_C* u_smooth_HOST,
                                             const RealArray_C* v_smooth_HOST,
                                             const RealArray_C* dz_HOST,
                                             RealArray_C* slope_x_HOST,
                                             RealArray_C* slope_y_HOST,
                                             const RealArray_C* sh_xx_bt_HOST,
                                             const RealArray_C* sh_xy_bt_HOST,
                                             const RealArray_C* GME_effic_h_HOST,
                                             const RealArray_C* GME_effic_q_HOST,
                                             const RealArray_C* KH_u_GME_HOST,
                                             const RealArray_C* KH_v_GME_HOST,
                                             RealArray_C* GME_coeff_h_HOST,
                                             RealArray_C* GME_coeff_q_HOST,
                                             RealArray_C* str_xx_GME_HOST,
                                             RealArray_C* str_xy_GME_HOST,
                                             RealArray_C* NoSt_HOST,
                                             RealArray_C* ShSt_HOST,
                                             RealArray_C* FrictWork_HOST,
                                             RealArray_C* FrictWork_bh_HOST,
                                             RealArray_C* FrictWork_GME_HOST,
                                             RealArray_C* Ah_h_HOST,
                                             RealArray_C* Ah_q_HOST,
                                             RealArray_C* Kh_h_HOST,
                                             RealArray_C* Kh_q_HOST,
                                             RealArray_C* grid_Re_Kh_HOST,
                                             RealArray_C* grid_Re_Ah_HOST,
                                             RealArray_C* div_xx_h_HOST,
                                             RealArray_C* sh_xx_h_HOST,
                                             RealArray_C* vort_xy_q_HOST,
                                             RealArray_C* sh_xy_q_HOST,
                                             RealArray_C* visc_limit_h_HOST,
                                             RealArray_C* visc_limit_q_HOST,
                                             RealArray_C* visc_limit_h_flag_HOST,
                                             RealArray_C* visc_limit_q_flag_HOST,
                                             RealArray_C* visc_limit_h_frac_HOST,
                                             RealArray_C* visc_limit_q_frac_HOST,
                                             RealArray_C* BS_coeff_h_HOST,
                                             RealArray_C* BS_coeff_q_HOST,
                                             const double H_subroundoff,
                                             const double H_to_RZ,
                                             const double s_to_T,
                                             void* host);

#ifdef __cplusplus
}
#endif
