// SKILLS: 0.3.2
#pragma once
/// @file turbotmp_mom_coriolisadv_bridge.h
/// @brief extern "C" bridge declarations between the MOM6 Fortran shims in
///        MOM_CoriolisAdv and the AMReX kernels (temporary turbotmp layer).

#include <stdbool.h>

#include "turbotmp_bridge_c_types.h"

/// Opaque: the OBC structure is never dereferenced on this side.
struct OceanOBC;
/// @brief Alias letting C callers spell the opaque OBC handle without the `struct` keyword.
typedef struct OceanOBC OceanOBC;
/// Opaque: the Stokes drift control structure is never dereferenced on this side.
struct WaveParametersCS;
/// @brief Alias letting C callers spell the opaque Stokes drift handle without the `struct` keyword.
typedef struct WaveParametersCS WaveParametersCS;
/// Defined in mom_coriolisadv.hpp -- field-for-field mirror of the Fortran bind(C) type
/// of the same name.
struct CorAdv_opts_CS_C;
/// @brief Alias letting C callers spell the options struct without the `struct` keyword.
typedef struct CorAdv_opts_CS_C CorAdv_opts_CS_C;
/// Defined in mom_coriolisadv.hpp -- field-for-field mirror of the Fortran bind(C) type
/// of the same name.
struct Coriolis_scheme_CS_C;
/// @brief Alias letting C callers spell the Coriolis scheme options without the `struct` keyword.
typedef struct Coriolis_scheme_CS_C Coriolis_scheme_CS_C;

/// @brief bind(C) mirror of CorAdv_AD_type, field-for-field, same order.
struct CorAdv_AD_C {
    struct RealArray_C gradKEu; ///< gradKEu = - d/dx(u2) [L T-2 ~> m s-2]
    struct RealArray_C gradKEv; ///< gradKEv = - d/dy(u2) [L T-2 ~> m s-2]
    struct RealArray_C rv_x_v;  ///< rv_x_v = rv * v at u [L T-2 ~> m s-2]
    struct RealArray_C rv_x_u;  ///< rv_x_u = rv * u at v [L T-2 ~> m s-2]
};
/// @brief Alias letting C callers spell the diagnostics mirror without the `struct` keyword.
typedef struct CorAdv_AD_C CorAdv_AD_C;

#ifdef __cplusplus
extern "C" {
#endif

void turbotmp_gradke_bridge(const Box_C* bxQ_HOST,
                            const Box_C* bxU_HOST,
                            const Box_C* bxV_HOST,
                            const RealArray_C* u_HOST,
                            const RealArray_C* v_HOST,
                            RealArray_C* KE_HOST,
                            RealArray_C* KEx_HOST,
                            RealArray_C* KEy_HOST,
                            const grid_core_C* Gcore_HOST,
                            const grid_OBC_C* Gobc_HOST,
                            const int KE_Scheme,
                            const bool KE_use_limiter);

void turbotmp_coradv_setup_bridge(const Box_C* bxH_HOST,
                                  const Box_C* bxQ_HOST,
                                  const Box_C* bxU_HOST,
                                  const Box_C* bxV_HOST,
                                  const Box_C* bxQs_HOST,
                                  const RealArray_C* u_HOST,
                                  const RealArray_C* v_HOST,
                                  const RealArray_C* h_HOST,
                                  OceanOBC* OBC,
                                  const CorAdv_AD_C* AD_HOST,
                                  const RealArray_C* por_face_areaU_HOST,
                                  const RealArray_C* por_face_areaV_HOST,
                                  WaveParametersCS* Waves,
                                  const grid_core_C* Gcore_HOST,
                                  const RealArray_C* Area_h_HOST,
                                  const RealArray_C* Area_q_HOST,
                                  const bool Stokes_VF,
                                  const bool use_weno,
                                  const double vol_neglect,
                                  const double area_neglect,
                                  RealArray_C* q_HOST,
                                  RealArray_C* Ih_q_HOST,
                                  RealArray_C* abs_vort_HOST,
                                  RealArray_C* h_q_HOST,
                                  RealArray_C* qS_HOST,
                                  RealArray_C* q2_HOST,
                                  RealArray_C* uh_center_HOST,
                                  RealArray_C* vh_center_HOST,
                                  RealArray_C* RV_HOST,
                                  RealArray_C* PV_HOST,
                                  const CorAdv_opts_CS_C* CS_HOST);

void turbotmp_coradv_common_terms_bridge(const Box_C* bxU_HOST,
                                         const Box_C* bxV_HOST,
                                         const RealArray_C* u_HOST,
                                         const RealArray_C* v_HOST,
                                         const RealArray_C* uh_HOST,
                                         const RealArray_C* vh_HOST,
                                         const RealArray_C* abs_vort_HOST,
                                         const RealArray_C* qS_HOST,
                                         const RealArray_C* q2_HOST,
                                         const RealArray_C* KEx_HOST,
                                         const RealArray_C* KEy_HOST,
                                         const bool Stokes_VF,
                                         const grid_core_C* Gcore_HOST,
                                         RealArray_C* CAu_HOST,
                                         RealArray_C* CAv_HOST,
                                         RealArray_C* CAuS_HOST,
                                         RealArray_C* CAvS_HOST,
                                         CorAdv_AD_C* AD_HOST,
                                         const CorAdv_opts_CS_C* CS_HOST);

void turbotmp_coradv_arakawa_bridge(const Box_C* bxH_HOST,
                                    const Box_C* bxQ_HOST,
                                    const Box_C* bxU_HOST,
                                    const Box_C* bxV_HOST,
                                    const RealArray_C* uh_HOST,
                                    const RealArray_C* vh_HOST,
                                    const RealArray_C* q_HOST,
                                    const RealArray_C* Ih_q_HOST,
                                    const double Fe_m2,
                                    const double rat_lin,
                                    const grid_core_C* Gcore_HOST,
                                    RealArray_C* CAu_HOST,
                                    RealArray_C* CAv_HOST,
                                    const double wt_lin_blend,
                                    const Coriolis_scheme_CS_C* CS_HOST);

void turbotmp_coradv_robust_enstro_bridge(const Box_C* bxU_HOST,
                                          const Box_C* bxV_HOST,
                                          const RealArray_C* u_HOST,
                                          const RealArray_C* v_HOST,
                                          const RealArray_C* h_HOST,
                                          const RealArray_C* uh_HOST,
                                          const RealArray_C* vh_HOST,
                                          const RealArray_C* abs_vort_HOST,
                                          const double eps_vel,
                                          const double h_tiny,
                                          const grid_core_C* Gcore_HOST,
                                          RealArray_C* CAu_HOST,
                                          RealArray_C* CAv_HOST,
                                          const int PV_Adv_Scheme);

void turbotmp_coradv_sadourny_bridge(const Box_C* bxH_HOST,
                                     const Box_C* bxU_HOST,
                                     const Box_C* bxV_HOST,
                                     const RealArray_C* u_HOST,
                                     const RealArray_C* v_HOST,
                                     const RealArray_C* uh_HOST,
                                     const RealArray_C* vh_HOST,
                                     const RealArray_C* q_HOST,
                                     const RealArray_C* uh_center_HOST,
                                     const RealArray_C* vh_center_HOST,
                                     const grid_core_C* Gcore_HOST,
                                     RealArray_C* CAu_HOST,
                                     RealArray_C* CAv_HOST,
                                     const Coriolis_scheme_CS_C* CS_HOST);

void turbotmp_coradv_weno_bridge(const Box_C* bxU_HOST,
                                 const Box_C* bxV_HOST,
                                 const RealArray_C* u_HOST,
                                 const RealArray_C* v_HOST,
                                 const RealArray_C* uh_HOST,
                                 const RealArray_C* vh_HOST,
                                 const RealArray_C* q_HOST,
                                 const RealArray_C* abs_vort_HOST,
                                 const RealArray_C* h_q_HOST,
                                 const grid_core_C* Gcore_HOST,
                                 const double H_subroundoff,
                                 RealArray_C* CAu_HOST,
                                 RealArray_C* CAv_HOST,
                                 const bool weno_velocity_smooth,
                                 const Coriolis_scheme_CS_C* CS_HOST);

void turbotmp_coradcalc_tr_bridge(const Box_C* bxH0_HOST,
                                  const Box_C* bxQ0_HOST,
                                  const RealArray_C* u_HOST,
                                  const RealArray_C* v_HOST,
                                  const RealArray_C* h_HOST,
                                  const RealArray_C* uh_HOST,
                                  const RealArray_C* vh_HOST,
                                  RealArray_C* CAu_HOST,
                                  RealArray_C* CAv_HOST,
                                  RealArray_C* RV_HOST,
                                  RealArray_C* PV_HOST,
                                  RealArray_C* CAuS_HOST,
                                  RealArray_C* CAvS_HOST,
                                  OceanOBC* OBC,
                                  CorAdv_AD_C* AD_HOST,
                                  const grid_core_C* Gcore_HOST,
                                  const grid_OBC_C* Gobc_HOST,
                                  const int nz,
                                  const double H_subroundoff,
                                  const double Angstrom_H,
                                  const double m_to_L,
                                  const double m_s_to_L_T,
                                  const CorAdv_opts_CS_C* CS_HOST,
                                  const RealArray_C* por_face_areaU_HOST,
                                  const RealArray_C* por_face_areaV_HOST,
                                  WaveParametersCS* Waves,
                                  const bool Stokes_VF);

#ifdef __cplusplus
}
#endif
