// SKILLS: 0.3.2
/// @file turbotmp_mom_coriolisadv_bridge.cpp
/// @brief Marshalling bridges between the MOM6 Fortran shims in MOM_CoriolisAdv
///        (mode = AMREX) and the AMReX kernels in mom_coriolisadv.cpp.

#include "mom_coriolisadv.hpp"
#include "turbotmp_helper.hpp"
#include "turbotmp_mom_coriolisadv_bridge.h"

using namespace amrex;

/// @brief Bridge for the MOM6 gradKE subroutine.
///
/// @param bxQ_HOST       The B-grid-index iteration box of this tile
/// @param bxU_HOST       The u-point iteration box of this tile
/// @param bxV_HOST       The v-point iteration box of this tile
/// @param u_HOST         Zonal velocity [L T-1 ~> m s-1]
/// @param v_HOST         Meridional velocity [L T-1 ~> m s-1]
/// @param KE_HOST        Kinetic energy per unit mass [L2 T-2 ~> m2 s-2]
/// @param KEx_HOST       Zonal acceleration due to kinetic energy gradient [L T-2 ~> m s-2]
/// @param KEy_HOST       Meridional acceleration due to kinetic energy gradient [L T-2 ~> m s-2]
/// @param Gcore_HOST     Persistent copies of the core ocean grid fields
/// @param Gobc_HOST      Persistent copies of the ocean grid OBC mask fields
/// @param KE_Scheme      Selects the discretization for the kinetic energy: KE_ARAKAWA,
///                       KE_SIMPLE_GUDONOV, KE_GUDONOV or KE_UP3
/// @param KE_use_limiter If true, use the Koren limiter for the KE_UP3 scheme
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
                            const bool KE_use_limiter)
{
    /// Reconstruct the iteration Boxes (Fortran 1-based -> AMReX 0-based)
    Box bxQ = turbotmp::to_box(bxQ_HOST);
    Box bxU = turbotmp::to_box(bxU_HOST);
    Box bxV = turbotmp::to_box(bxV_HOST);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto u_DEV             = turbotmp::make_array4(u_HOST);
    auto v_DEV             = turbotmp::make_array4(v_HOST);
    auto KE_DEV            = turbotmp::make_array4(KE_HOST);
    auto KEx_DEV           = turbotmp::make_array4(KEx_HOST);
    auto KEy_DEV           = turbotmp::make_array4(KEy_HOST);
    auto areaCu_DEV        = turbotmp::make_array4(&Gcore_HOST->areaCu);
    auto areaCv_DEV        = turbotmp::make_array4(&Gcore_HOST->areaCv);
    auto IareaT_DEV        = turbotmp::make_array4(&Gcore_HOST->IareaT);
    auto mask2dCu_DEV      = turbotmp::make_array4(&Gcore_HOST->mask2dCu);
    auto mask2dCv_DEV      = turbotmp::make_array4(&Gcore_HOST->mask2dCv);
    auto IdxCu_OBCmask_DEV = turbotmp::make_array4(&Gobc_HOST->IdxCu_OBCmask);
    auto IdyCv_OBCmask_DEV = turbotmp::make_array4(&Gobc_HOST->IdyCv_OBCmask);

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(u_HOST->data,                    u_DEV);
    turbotmp::copy_FortranHost_to_array4(v_HOST->data,                    v_DEV);
    turbotmp::copy_FortranHost_to_array4(KE_HOST->data,                   KE_DEV);
    turbotmp::copy_FortranHost_to_array4(KEx_HOST->data,                  KEx_DEV);
    turbotmp::copy_FortranHost_to_array4(KEy_HOST->data,                  KEy_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->areaCu.data,         areaCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->areaCv.data,         areaCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IareaT.data,         IareaT_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dCu.data,       mask2dCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dCv.data,       mask2dCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gobc_HOST->IdxCu_OBCmask.data,   IdxCu_OBCmask_DEV);
    turbotmp::copy_FortranHost_to_array4(Gobc_HOST->IdyCv_OBCmask.data,   IdyCv_OBCmask_DEV);

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::gradke(bxQ, bxU, bxV,
                u_DEV.arr, v_DEV.arr,
                KE_DEV.arr, KEx_DEV.arr, KEy_DEV.arr,
                areaCu_DEV.arr, areaCv_DEV.arr, IareaT_DEV.arr,
                mask2dCu_DEV.arr, mask2dCv_DEV.arr,
                IdxCu_OBCmask_DEV.arr, IdyCv_OBCmask_DEV.arr,
                KE_Scheme, KE_use_limiter);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(KE_DEV,  KE_HOST->data);
    turbotmp::copy_array4_to_FortranHost(KEx_DEV, KEx_HOST->data);
    turbotmp::copy_array4_to_FortranHost(KEy_DEV, KEy_HOST->data);

    /// Free everything we allocated
    turbotmp::free_array4(u_DEV);
    turbotmp::free_array4(v_DEV);
    turbotmp::free_array4(KE_DEV);
    turbotmp::free_array4(KEx_DEV);
    turbotmp::free_array4(KEy_DEV);
    turbotmp::free_array4(areaCu_DEV);
    turbotmp::free_array4(areaCv_DEV);
    turbotmp::free_array4(IareaT_DEV);
    turbotmp::free_array4(mask2dCu_DEV);
    turbotmp::free_array4(mask2dCv_DEV);
    turbotmp::free_array4(IdxCu_OBCmask_DEV);
    turbotmp::free_array4(IdyCv_OBCmask_DEV);
}

/// @brief Bridge for the MOM6 CorAdv_setup subroutine.
///
/// @param bxH_HOST            The h-point iteration box of this tile
/// @param bxQ_HOST            The B-grid-index iteration box of this tile
/// @param bxU_HOST            The u-point iteration box of this tile
/// @param bxV_HOST            The v-point iteration box of this tile
/// @param bxQs_HOST           The scheme-dependent vorticity-point box
/// @param u_HOST              Zonal velocity [L T-1 ~> m s-1]
/// @param v_HOST              Meridional velocity [L T-1 ~> m s-1]
/// @param h_HOST              Layer thickness [H ~> m or kg m-2]
/// @param OBC                 Open boundary control structure
/// @param AD_HOST             Copies of the acceleration diagnostics
/// @param por_face_areaU_HOST fractional open area of U-faces [nondim]
/// @param por_face_areaV_HOST fractional open area of V-faces [nondim]
/// @param Waves               An optional pointer to Stokes drift CS
/// @param Gcore_HOST          Persistent copies of the core ocean grid fields
/// @param Area_h_HOST         The ocean area at h points [L2 ~> m2]
/// @param Area_q_HOST         The sum of the ocean areas at the 4 adjacent thickness points [L2 ~> m2]
/// @param Stokes_VF           If true, include the Stokes drift
/// @param use_weno            True if using one of the WENO schemes
/// @param vol_neglect         A negligible volume [H L2 ~> m3 or kg]
/// @param area_neglect        A negligible area [L2 ~> m2]
/// @param q_HOST              Layer potential vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
/// @param Ih_q_HOST           Inverse of thickness interpolated to q points [H-1 ~> m-1 or m2 kg-1]
/// @param abs_vort_HOST       Absolute vorticity at q-points [T-1 ~> s-1]
/// @param h_q_HOST            Thickness interpolated to q points [H ~> m or kg m-2]
/// @param qS_HOST             Layer Stokes vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
/// @param q2_HOST             Relative vorticity over thickness [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
/// @param uh_center_HOST      Transport based on arithmetic mean h at u-points
///                            [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param vh_center_HOST      Transport based on arithmetic mean h at v-points
///                            [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param RV_HOST             Diagnostic relative vorticity [T-1 ~> s-1]
/// @param PV_HOST             Diagnostic potential vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
/// @param CS_HOST             Options used in the Coriolis and momentum advection calculations
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
                                  const CorAdv_opts_CS_C* CS_HOST)
{
    /// Reconstruct the iteration Boxes (Fortran 1-based -> AMReX 0-based)
    Box bxH  = turbotmp::to_box(bxH_HOST);
    Box bxQ  = turbotmp::to_box(bxQ_HOST);
    Box bxU  = turbotmp::to_box(bxU_HOST);
    Box bxV  = turbotmp::to_box(bxV_HOST);
    Box bxQs = turbotmp::to_box(bxQs_HOST);

    /// Only the presence of the rv_x_v / rv_x_u diagnostics is used (data == nullptr when absent)
    const bool rv_x_v_associated = (AD_HOST->rv_x_v.data != nullptr);
    const bool rv_x_u_associated = (AD_HOST->rv_x_u.data != nullptr);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto u_DEV              = turbotmp::make_array4(u_HOST);
    auto v_DEV              = turbotmp::make_array4(v_HOST);
    auto h_DEV              = turbotmp::make_array4(h_HOST);
    auto por_face_areaU_DEV = turbotmp::make_array4(por_face_areaU_HOST);
    auto por_face_areaV_DEV = turbotmp::make_array4(por_face_areaV_HOST);
    auto dyCv_DEV           = turbotmp::make_array4(&Gcore_HOST->dyCv);
    auto dxCu_DEV           = turbotmp::make_array4(&Gcore_HOST->dxCu);
    auto dy_Cu_DEV          = turbotmp::make_array4(&Gcore_HOST->dy_Cu);
    auto dx_Cv_DEV          = turbotmp::make_array4(&Gcore_HOST->dx_Cv);
    auto dyBu_DEV           = turbotmp::make_array4(&Gcore_HOST->dyBu);
    auto dxBu_DEV           = turbotmp::make_array4(&Gcore_HOST->dxBu);
    auto mask2dBu_DEV       = turbotmp::make_array4(&Gcore_HOST->mask2dBu);
    auto IareaBu_DEV        = turbotmp::make_array4(&Gcore_HOST->IareaBu);
    auto CoriolisBu_DEV     = turbotmp::make_array4(&Gcore_HOST->CoriolisBu);
    auto Area_h_DEV         = turbotmp::make_array4(Area_h_HOST);
    auto Area_q_DEV         = turbotmp::make_array4(Area_q_HOST);
    auto q_DEV              = turbotmp::make_array4(q_HOST);
    auto Ih_q_DEV           = turbotmp::make_array4(Ih_q_HOST);
    auto abs_vort_DEV       = turbotmp::make_array4(abs_vort_HOST);
    auto h_q_DEV            = turbotmp::make_array4(h_q_HOST);
    auto qS_DEV             = turbotmp::make_array4(qS_HOST);
    auto q2_DEV             = turbotmp::make_array4(q2_HOST);
    auto uh_center_DEV      = turbotmp::make_array4(uh_center_HOST);
    auto vh_center_DEV      = turbotmp::make_array4(vh_center_HOST);
    auto RV_DEV             = turbotmp::make_array4(RV_HOST);
    auto PV_DEV             = turbotmp::make_array4(PV_HOST);

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(u_HOST->data,                      u_DEV);
    turbotmp::copy_FortranHost_to_array4(v_HOST->data,                      v_DEV);
    turbotmp::copy_FortranHost_to_array4(h_HOST->data,                      h_DEV);
    turbotmp::copy_FortranHost_to_array4(por_face_areaU_HOST->data,         por_face_areaU_DEV);
    turbotmp::copy_FortranHost_to_array4(por_face_areaV_HOST->data,         por_face_areaV_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->dyCv.data,             dyCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->dxCu.data,             dxCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->dy_Cu.data,            dy_Cu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->dx_Cv.data,            dx_Cv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->dyBu.data,             dyBu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->dxBu.data,             dxBu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dBu.data,         mask2dBu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IareaBu.data,          IareaBu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->CoriolisBu.data,       CoriolisBu_DEV);
    turbotmp::copy_FortranHost_to_array4(Area_h_HOST->data,                 Area_h_DEV);
    turbotmp::copy_FortranHost_to_array4(Area_q_HOST->data,                 Area_q_DEV);
    turbotmp::copy_FortranHost_to_array4(q_HOST->data,                      q_DEV);
    turbotmp::copy_FortranHost_to_array4(Ih_q_HOST->data,                   Ih_q_DEV);
    turbotmp::copy_FortranHost_to_array4(abs_vort_HOST->data,               abs_vort_DEV);
    turbotmp::copy_FortranHost_to_array4(h_q_HOST->data,                    h_q_DEV);
    turbotmp::copy_FortranHost_to_array4(qS_HOST->data,                     qS_DEV);
    turbotmp::copy_FortranHost_to_array4(q2_HOST->data,                     q2_DEV);
    turbotmp::copy_FortranHost_to_array4(uh_center_HOST->data,              uh_center_DEV);
    turbotmp::copy_FortranHost_to_array4(vh_center_HOST->data,              vh_center_DEV);
    turbotmp::copy_FortranHost_to_array4(RV_HOST->data,                     RV_DEV);
    turbotmp::copy_FortranHost_to_array4(PV_HOST->data,                     PV_DEV);

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::coradv_setup(bxH, bxQ, bxU, bxV, bxQs,
                      u_DEV.arr, v_DEV.arr, h_DEV.arr,
                      OBC,
                      rv_x_v_associated, rv_x_u_associated,
                      por_face_areaU_DEV.arr, por_face_areaV_DEV.arr,
                      Waves,
                      dyCv_DEV.arr, dxCu_DEV.arr, dy_Cu_DEV.arr, dx_Cv_DEV.arr,
                      dyBu_DEV.arr, dxBu_DEV.arr, mask2dBu_DEV.arr, IareaBu_DEV.arr,
                      CoriolisBu_DEV.arr,
                      Area_h_DEV.arr, Area_q_DEV.arr,
                      Stokes_VF, use_weno, vol_neglect, area_neglect,
                      q_DEV.arr, Ih_q_DEV.arr, abs_vort_DEV.arr, h_q_DEV.arr,
                      qS_DEV.arr, q2_DEV.arr, uh_center_DEV.arr, vh_center_DEV.arr,
                      RV_DEV.arr, PV_DEV.arr,
                      *CS_HOST);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(q_DEV,              q_HOST->data);
    turbotmp::copy_array4_to_FortranHost(Ih_q_DEV,           Ih_q_HOST->data);
    turbotmp::copy_array4_to_FortranHost(abs_vort_DEV,       abs_vort_HOST->data);
    turbotmp::copy_array4_to_FortranHost(h_q_DEV,            h_q_HOST->data);
    turbotmp::copy_array4_to_FortranHost(qS_DEV,             qS_HOST->data);
    turbotmp::copy_array4_to_FortranHost(q2_DEV,             q2_HOST->data);
    turbotmp::copy_array4_to_FortranHost(uh_center_DEV,      uh_center_HOST->data);
    turbotmp::copy_array4_to_FortranHost(vh_center_DEV,      vh_center_HOST->data);
    turbotmp::copy_array4_to_FortranHost(RV_DEV,             RV_HOST->data);
    turbotmp::copy_array4_to_FortranHost(PV_DEV,             PV_HOST->data);

    /// Free everything we allocated
    turbotmp::free_array4(u_DEV);
    turbotmp::free_array4(v_DEV);
    turbotmp::free_array4(h_DEV);
    turbotmp::free_array4(por_face_areaU_DEV);
    turbotmp::free_array4(por_face_areaV_DEV);
    turbotmp::free_array4(dyCv_DEV);
    turbotmp::free_array4(dxCu_DEV);
    turbotmp::free_array4(dy_Cu_DEV);
    turbotmp::free_array4(dx_Cv_DEV);
    turbotmp::free_array4(dyBu_DEV);
    turbotmp::free_array4(dxBu_DEV);
    turbotmp::free_array4(mask2dBu_DEV);
    turbotmp::free_array4(IareaBu_DEV);
    turbotmp::free_array4(CoriolisBu_DEV);
    turbotmp::free_array4(Area_h_DEV);
    turbotmp::free_array4(Area_q_DEV);
    turbotmp::free_array4(q_DEV);
    turbotmp::free_array4(Ih_q_DEV);
    turbotmp::free_array4(abs_vort_DEV);
    turbotmp::free_array4(h_q_DEV);
    turbotmp::free_array4(qS_DEV);
    turbotmp::free_array4(q2_DEV);
    turbotmp::free_array4(uh_center_DEV);
    turbotmp::free_array4(vh_center_DEV);
    turbotmp::free_array4(RV_DEV);
    turbotmp::free_array4(PV_DEV);
}

/// @brief Bridge for the MOM6 CorAdv_common_terms subroutine.
///
/// @param bxU_HOST      The u-point iteration box of this tile
/// @param bxV_HOST      The v-point iteration box of this tile
/// @param u_HOST        Zonal velocity [L T-1 ~> m s-1]
/// @param v_HOST        Meridional velocity [L T-1 ~> m s-1]
/// @param uh_HOST       Zonal transport u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param vh_HOST       Meridional transport v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param abs_vort_HOST Absolute vorticity at q-points [T-1 ~> s-1]
/// @param qS_HOST       Layer Stokes vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
/// @param q2_HOST       Relative vorticity over thickness [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
/// @param KEx_HOST      Zonal gradient of kinetic energy per unit mass [L T-2 ~> m s-2]
/// @param KEy_HOST      Meridional gradient of kinetic energy per unit mass [L T-2 ~> m s-2]
/// @param Stokes_VF     If true, include the Stokes drift
/// @param Gcore_HOST    Persistent copies of the core ocean grid fields
/// @param CAu_HOST      Zonal acceleration due to Coriolis and momentum advection [L T-2 ~> m s-2].
/// @param CAv_HOST      Meridional acceleration due to Coriolis and momentum advection [L T-2 ~> m s-2].
/// @param CAuS_HOST     Stokes contribution to CAu [L T-2 ~> m s-2]
/// @param CAvS_HOST     Stokes contribution to CAv [L T-2 ~> m s-2]
/// @param AD_HOST       Copies of the acceleration diagnostics; each member may be absent
///                      (data == nullptr)
/// @param CS_HOST       Options used in the Coriolis and momentum advection calculations
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
                                         const CorAdv_opts_CS_C* CS_HOST)
{
    /// Reconstruct the iteration Boxes (Fortran 1-based -> AMReX 0-based)
    Box bxU = turbotmp::to_box(bxU_HOST);
    Box bxV = turbotmp::to_box(bxV_HOST);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto u_DEV        = turbotmp::make_array4(u_HOST);
    auto v_DEV        = turbotmp::make_array4(v_HOST);
    auto uh_DEV       = turbotmp::make_array4(uh_HOST);
    auto vh_DEV       = turbotmp::make_array4(vh_HOST);
    auto abs_vort_DEV = turbotmp::make_array4(abs_vort_HOST);
    auto qS_DEV       = turbotmp::make_array4(qS_HOST);
    auto q2_DEV       = turbotmp::make_array4(q2_HOST);
    auto KEx_DEV      = turbotmp::make_array4(KEx_HOST);
    auto KEy_DEV      = turbotmp::make_array4(KEy_HOST);
    auto IdxCu_DEV    = turbotmp::make_array4(&Gcore_HOST->IdxCu);
    auto IdyCv_DEV    = turbotmp::make_array4(&Gcore_HOST->IdyCv);
    auto CAu_DEV      = turbotmp::make_array4(CAu_HOST);
    auto CAv_DEV      = turbotmp::make_array4(CAv_HOST);
    auto CAuS_DEV     = turbotmp::make_array4(CAuS_HOST);
    auto CAvS_DEV     = turbotmp::make_array4(CAvS_HOST);

    /// The acceleration diagnostics may be absent (data == nullptr); only allocate/copy them when present.
    const bool has_gradKEu  = (AD_HOST->gradKEu.data != nullptr);
    const bool has_gradKEv  = (AD_HOST->gradKEv.data != nullptr);
    const bool has_rv_x_v   = (AD_HOST->rv_x_v.data != nullptr);
    const bool has_rv_x_u   = (AD_HOST->rv_x_u.data != nullptr);
    turbotmp::A4Box gradKEu_DEV {};
    turbotmp::A4Box gradKEv_DEV {};
    turbotmp::A4Box rv_x_v_DEV  {};
    turbotmp::A4Box rv_x_u_DEV  {};
    if (has_gradKEu) {
        gradKEu_DEV = turbotmp::make_array4(&AD_HOST->gradKEu);
    }
    if (has_gradKEv) {
        gradKEv_DEV = turbotmp::make_array4(&AD_HOST->gradKEv);
    }
    if (has_rv_x_v) {
        rv_x_v_DEV = turbotmp::make_array4(&AD_HOST->rv_x_v);
    }
    if (has_rv_x_u) {
        rv_x_u_DEV = turbotmp::make_array4(&AD_HOST->rv_x_u);
    }

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(u_HOST->data,                  u_DEV);
    turbotmp::copy_FortranHost_to_array4(v_HOST->data,                  v_DEV);
    turbotmp::copy_FortranHost_to_array4(uh_HOST->data,                 uh_DEV);
    turbotmp::copy_FortranHost_to_array4(vh_HOST->data,                 vh_DEV);
    turbotmp::copy_FortranHost_to_array4(abs_vort_HOST->data,           abs_vort_DEV);
    turbotmp::copy_FortranHost_to_array4(qS_HOST->data,                 qS_DEV);
    turbotmp::copy_FortranHost_to_array4(q2_HOST->data,                 q2_DEV);
    turbotmp::copy_FortranHost_to_array4(KEx_HOST->data,                KEx_DEV);
    turbotmp::copy_FortranHost_to_array4(KEy_HOST->data,                KEy_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdxCu.data,        IdxCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdyCv.data,        IdyCv_DEV);
    turbotmp::copy_FortranHost_to_array4(CAu_HOST->data,                CAu_DEV);
    turbotmp::copy_FortranHost_to_array4(CAv_HOST->data,                CAv_DEV);
    turbotmp::copy_FortranHost_to_array4(CAuS_HOST->data,               CAuS_DEV);
    turbotmp::copy_FortranHost_to_array4(CAvS_HOST->data,               CAvS_DEV);
    if (has_gradKEu) {
        turbotmp::copy_FortranHost_to_array4(AD_HOST->gradKEu.data, gradKEu_DEV);
    }
    if (has_gradKEv) {
        turbotmp::copy_FortranHost_to_array4(AD_HOST->gradKEv.data, gradKEv_DEV);
    }
    if (has_rv_x_v) {
        turbotmp::copy_FortranHost_to_array4(AD_HOST->rv_x_v.data, rv_x_v_DEV);
    }
    if (has_rv_x_u) {
        turbotmp::copy_FortranHost_to_array4(AD_HOST->rv_x_u.data, rv_x_u_DEV);
    }

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::coradv_common_terms(bxU, bxV,
                             u_DEV.arr, v_DEV.arr, uh_DEV.arr, vh_DEV.arr,
                             abs_vort_DEV.arr, qS_DEV.arr, q2_DEV.arr,
                             KEx_DEV.arr, KEy_DEV.arr,
                             Stokes_VF,
                             IdxCu_DEV.arr, IdyCv_DEV.arr,
                             CAu_DEV.arr, CAv_DEV.arr, CAuS_DEV.arr, CAvS_DEV.arr,
                             gradKEu_DEV.arr, gradKEv_DEV.arr, rv_x_v_DEV.arr, rv_x_u_DEV.arr,
                             *CS_HOST);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(CAu_DEV,      CAu_HOST->data);
    turbotmp::copy_array4_to_FortranHost(CAv_DEV,      CAv_HOST->data);
    turbotmp::copy_array4_to_FortranHost(CAuS_DEV,     CAuS_HOST->data);
    turbotmp::copy_array4_to_FortranHost(CAvS_DEV,     CAvS_HOST->data);
    if (has_gradKEu) {
        turbotmp::copy_array4_to_FortranHost(gradKEu_DEV, AD_HOST->gradKEu.data);
    }
    if (has_gradKEv) {
        turbotmp::copy_array4_to_FortranHost(gradKEv_DEV, AD_HOST->gradKEv.data);
    }
    if (has_rv_x_v) {
        turbotmp::copy_array4_to_FortranHost(rv_x_v_DEV, AD_HOST->rv_x_v.data);
    }
    if (has_rv_x_u) {
        turbotmp::copy_array4_to_FortranHost(rv_x_u_DEV, AD_HOST->rv_x_u.data);
    }

    /// Free everything we allocated (free_array4 on a never-allocated A4Box is a safe no-op)
    turbotmp::free_array4(u_DEV);
    turbotmp::free_array4(v_DEV);
    turbotmp::free_array4(uh_DEV);
    turbotmp::free_array4(vh_DEV);
    turbotmp::free_array4(abs_vort_DEV);
    turbotmp::free_array4(qS_DEV);
    turbotmp::free_array4(q2_DEV);
    turbotmp::free_array4(KEx_DEV);
    turbotmp::free_array4(KEy_DEV);
    turbotmp::free_array4(IdxCu_DEV);
    turbotmp::free_array4(IdyCv_DEV);
    turbotmp::free_array4(CAu_DEV);
    turbotmp::free_array4(CAv_DEV);
    turbotmp::free_array4(CAuS_DEV);
    turbotmp::free_array4(CAvS_DEV);
    turbotmp::free_array4(gradKEu_DEV);
    turbotmp::free_array4(gradKEv_DEV);
    turbotmp::free_array4(rv_x_v_DEV);
    turbotmp::free_array4(rv_x_u_DEV);
}

/// @brief Bridge for the MOM6 CorAdv_arakawa subroutine.
///
/// @param bxH_HOST     The h-point iteration box of this tile
/// @param bxQ_HOST     The B-grid-index iteration box of this tile
/// @param bxU_HOST     The u-point iteration box of this tile
/// @param bxV_HOST     The v-point iteration box of this tile
/// @param uh_HOST      Zonal transport u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param vh_HOST      Meridional transport v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param q_HOST       Layer potential vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
/// @param Ih_q_HOST    Inverse of thickness interpolated to q points [H-1 ~> m-1 or m2 kg-1]
/// @param Fe_m2        ARAKAWA_LAMB_BLEND blending parameter [nondim]
/// @param rat_lin      ARAKAWA_LAMB_BLEND blending parameter [nondim]
/// @param Gcore_HOST   Persistent copies of the core ocean grid fields
/// @param CAu_HOST     Zonal acceleration due to Coriolis and momentum advection [L T-2 ~> m s-2].
/// @param CAv_HOST     Meridional acceleration due to Coriolis and momentum advection [L T-2 ~> m s-2].
/// @param wt_lin_blend A weighting value beyond which the blending between Sadourny and
///                     Arakawa & Hsu goes linearly to 0 [nondim].
/// @param CS_HOST      Options selecting the Coriolis discretization
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
                                    const Coriolis_scheme_CS_C* CS_HOST)
{
    /// Reconstruct the iteration Boxes (Fortran 1-based -> AMReX 0-based)
    Box bxH = turbotmp::to_box(bxH_HOST);
    Box bxQ = turbotmp::to_box(bxQ_HOST);
    Box bxU = turbotmp::to_box(bxU_HOST);
    Box bxV = turbotmp::to_box(bxV_HOST);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto uh_DEV    = turbotmp::make_array4(uh_HOST);
    auto vh_DEV    = turbotmp::make_array4(vh_HOST);
    auto q_DEV     = turbotmp::make_array4(q_HOST);
    auto Ih_q_DEV  = turbotmp::make_array4(Ih_q_HOST);
    auto IdxCu_DEV = turbotmp::make_array4(&Gcore_HOST->IdxCu);
    auto IdyCv_DEV = turbotmp::make_array4(&Gcore_HOST->IdyCv);
    auto CAu_DEV   = turbotmp::make_array4(CAu_HOST);
    auto CAv_DEV   = turbotmp::make_array4(CAv_HOST);

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(uh_HOST->data,              uh_DEV);
    turbotmp::copy_FortranHost_to_array4(vh_HOST->data,              vh_DEV);
    turbotmp::copy_FortranHost_to_array4(q_HOST->data,               q_DEV);
    turbotmp::copy_FortranHost_to_array4(Ih_q_HOST->data,            Ih_q_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdxCu.data,     IdxCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdyCv.data,     IdyCv_DEV);
    turbotmp::copy_FortranHost_to_array4(CAu_HOST->data,             CAu_DEV);
    turbotmp::copy_FortranHost_to_array4(CAv_HOST->data,             CAv_DEV);

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::coradv_arakawa(bxH, bxQ, bxU, bxV,
                        uh_DEV.arr, vh_DEV.arr, q_DEV.arr, Ih_q_DEV.arr,
                        Fe_m2, rat_lin,
                        IdxCu_DEV.arr, IdyCv_DEV.arr,
                        CAu_DEV.arr, CAv_DEV.arr,
                        wt_lin_blend,
                        *CS_HOST);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(CAu_DEV,   CAu_HOST->data);
    turbotmp::copy_array4_to_FortranHost(CAv_DEV,   CAv_HOST->data);

    /// Free everything we allocated
    turbotmp::free_array4(uh_DEV);
    turbotmp::free_array4(vh_DEV);
    turbotmp::free_array4(q_DEV);
    turbotmp::free_array4(Ih_q_DEV);
    turbotmp::free_array4(IdxCu_DEV);
    turbotmp::free_array4(IdyCv_DEV);
    turbotmp::free_array4(CAu_DEV);
    turbotmp::free_array4(CAv_DEV);
}

/// @brief Bridge for the MOM6 CorAdv_robust_enstro subroutine.
///
/// @param bxU_HOST      The u-point iteration box of this tile
/// @param bxV_HOST      The v-point iteration box of this tile
/// @param u_HOST        Zonal velocity [L T-1 ~> m s-1]
/// @param v_HOST        Meridional velocity [L T-1 ~> m s-1]
/// @param h_HOST        Layer thickness [H ~> m or kg m-2]
/// @param uh_HOST       Zonal transport u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param vh_HOST       Meridional transport v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param abs_vort_HOST Absolute vorticity at q-points [T-1 ~> s-1]
/// @param eps_vel       A tiny, positive velocity [L T-1 ~> m s-1]
/// @param h_tiny        A very small thickness [H ~> m or kg m-2]
/// @param Gcore_HOST    Persistent copies of the core ocean grid fields
/// @param CAu_HOST      Zonal acceleration due to Coriolis and momentum advection [L T-2 ~> m s-2].
/// @param CAv_HOST      Meridional acceleration due to Coriolis and momentum advection [L T-2 ~> m s-2].
/// @param PV_Adv_Scheme Selects the discretization for PV advection: PV_ADV_CENTERED or
///                      PV_ADV_UPWIND1
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
                                          const int PV_Adv_Scheme)
{
    /// Reconstruct the iteration Boxes (Fortran 1-based -> AMReX 0-based)
    Box bxU = turbotmp::to_box(bxU_HOST);
    Box bxV = turbotmp::to_box(bxV_HOST);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto u_DEV        = turbotmp::make_array4(u_HOST);
    auto v_DEV        = turbotmp::make_array4(v_HOST);
    auto h_DEV        = turbotmp::make_array4(h_HOST);
    auto uh_DEV       = turbotmp::make_array4(uh_HOST);
    auto vh_DEV       = turbotmp::make_array4(vh_HOST);
    auto abs_vort_DEV = turbotmp::make_array4(abs_vort_HOST);
    auto IdxCv_DEV    = turbotmp::make_array4(&Gcore_HOST->IdxCv);
    auto IdxCu_DEV    = turbotmp::make_array4(&Gcore_HOST->IdxCu);
    auto IdyCu_DEV    = turbotmp::make_array4(&Gcore_HOST->IdyCu);
    auto IdyCv_DEV    = turbotmp::make_array4(&Gcore_HOST->IdyCv);
    auto CAu_DEV      = turbotmp::make_array4(CAu_HOST);
    auto CAv_DEV      = turbotmp::make_array4(CAv_HOST);

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(u_HOST->data,                  u_DEV);
    turbotmp::copy_FortranHost_to_array4(v_HOST->data,                  v_DEV);
    turbotmp::copy_FortranHost_to_array4(h_HOST->data,                  h_DEV);
    turbotmp::copy_FortranHost_to_array4(uh_HOST->data,                 uh_DEV);
    turbotmp::copy_FortranHost_to_array4(vh_HOST->data,                 vh_DEV);
    turbotmp::copy_FortranHost_to_array4(abs_vort_HOST->data,           abs_vort_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdxCv.data,        IdxCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdxCu.data,        IdxCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdyCu.data,        IdyCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdyCv.data,        IdyCv_DEV);
    turbotmp::copy_FortranHost_to_array4(CAu_HOST->data,                CAu_DEV);
    turbotmp::copy_FortranHost_to_array4(CAv_HOST->data,                CAv_DEV);

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::coradv_robust_enstro(bxU, bxV,
                              u_DEV.arr, v_DEV.arr, h_DEV.arr,
                              uh_DEV.arr, vh_DEV.arr, abs_vort_DEV.arr,
                              eps_vel, h_tiny,
                              IdxCv_DEV.arr, IdxCu_DEV.arr, IdyCu_DEV.arr, IdyCv_DEV.arr,
                              CAu_DEV.arr, CAv_DEV.arr,
                              PV_Adv_Scheme);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(CAu_DEV,      CAu_HOST->data);
    turbotmp::copy_array4_to_FortranHost(CAv_DEV,      CAv_HOST->data);

    /// Free everything we allocated
    turbotmp::free_array4(u_DEV);
    turbotmp::free_array4(v_DEV);
    turbotmp::free_array4(h_DEV);
    turbotmp::free_array4(uh_DEV);
    turbotmp::free_array4(vh_DEV);
    turbotmp::free_array4(abs_vort_DEV);
    turbotmp::free_array4(IdxCv_DEV);
    turbotmp::free_array4(IdxCu_DEV);
    turbotmp::free_array4(IdyCu_DEV);
    turbotmp::free_array4(IdyCv_DEV);
    turbotmp::free_array4(CAu_DEV);
    turbotmp::free_array4(CAv_DEV);
}

/// @brief Bridge for the MOM6 CorAdv_sadourny subroutine.
///
/// @param bxH_HOST       The h-point iteration box of this tile
/// @param bxU_HOST       The u-point iteration box of this tile
/// @param bxV_HOST       The v-point iteration box of this tile
/// @param u_HOST         Zonal velocity [L T-1 ~> m s-1]
/// @param v_HOST         Meridional velocity [L T-1 ~> m s-1]
/// @param uh_HOST        Zonal transport u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param vh_HOST        Meridional transport v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param q_HOST         Layer potential vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
/// @param uh_center_HOST Transport based on arithmetic mean h at u-points
///                       [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param vh_center_HOST Transport based on arithmetic mean h at v-points
///                       [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param Gcore_HOST     Persistent copies of the core ocean grid fields
/// @param CAu_HOST       Zonal acceleration due to Coriolis and momentum advection [L T-2 ~> m s-2].
/// @param CAv_HOST       Meridional acceleration due to Coriolis and momentum advection [L T-2 ~> m s-2].
/// @param CS_HOST        Options selecting the Coriolis discretization
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
                                     const Coriolis_scheme_CS_C* CS_HOST)
{
    /// Reconstruct the iteration Boxes (Fortran 1-based -> AMReX 0-based)
    Box bxH = turbotmp::to_box(bxH_HOST);
    Box bxU = turbotmp::to_box(bxU_HOST);
    Box bxV = turbotmp::to_box(bxV_HOST);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto u_DEV         = turbotmp::make_array4(u_HOST);
    auto v_DEV         = turbotmp::make_array4(v_HOST);
    auto uh_DEV        = turbotmp::make_array4(uh_HOST);
    auto vh_DEV        = turbotmp::make_array4(vh_HOST);
    auto q_DEV         = turbotmp::make_array4(q_HOST);
    auto uh_center_DEV = turbotmp::make_array4(uh_center_HOST);
    auto vh_center_DEV = turbotmp::make_array4(vh_center_HOST);
    auto dy_Cu_DEV     = turbotmp::make_array4(&Gcore_HOST->dy_Cu);
    auto dx_Cv_DEV     = turbotmp::make_array4(&Gcore_HOST->dx_Cv);
    auto IdxCu_DEV     = turbotmp::make_array4(&Gcore_HOST->IdxCu);
    auto IdyCv_DEV     = turbotmp::make_array4(&Gcore_HOST->IdyCv);
    auto CAu_DEV       = turbotmp::make_array4(CAu_HOST);
    auto CAv_DEV       = turbotmp::make_array4(CAv_HOST);

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(u_HOST->data,                   u_DEV);
    turbotmp::copy_FortranHost_to_array4(v_HOST->data,                   v_DEV);
    turbotmp::copy_FortranHost_to_array4(uh_HOST->data,                  uh_DEV);
    turbotmp::copy_FortranHost_to_array4(vh_HOST->data,                  vh_DEV);
    turbotmp::copy_FortranHost_to_array4(q_HOST->data,                   q_DEV);
    turbotmp::copy_FortranHost_to_array4(uh_center_HOST->data,           uh_center_DEV);
    turbotmp::copy_FortranHost_to_array4(vh_center_HOST->data,           vh_center_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->dy_Cu.data,         dy_Cu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->dx_Cv.data,         dx_Cv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdxCu.data,         IdxCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdyCv.data,         IdyCv_DEV);
    turbotmp::copy_FortranHost_to_array4(CAu_HOST->data,                 CAu_DEV);
    turbotmp::copy_FortranHost_to_array4(CAv_HOST->data,                 CAv_DEV);

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::coradv_sadourny(bxH, bxU, bxV,
                         u_DEV.arr, v_DEV.arr, uh_DEV.arr, vh_DEV.arr, q_DEV.arr,
                         uh_center_DEV.arr, vh_center_DEV.arr,
                         dy_Cu_DEV.arr, dx_Cv_DEV.arr, IdxCu_DEV.arr, IdyCv_DEV.arr,
                         CAu_DEV.arr, CAv_DEV.arr,
                         *CS_HOST);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(CAu_DEV,       CAu_HOST->data);
    turbotmp::copy_array4_to_FortranHost(CAv_DEV,       CAv_HOST->data);

    /// Free everything we allocated
    turbotmp::free_array4(u_DEV);
    turbotmp::free_array4(v_DEV);
    turbotmp::free_array4(uh_DEV);
    turbotmp::free_array4(vh_DEV);
    turbotmp::free_array4(q_DEV);
    turbotmp::free_array4(uh_center_DEV);
    turbotmp::free_array4(vh_center_DEV);
    turbotmp::free_array4(dy_Cu_DEV);
    turbotmp::free_array4(dx_Cv_DEV);
    turbotmp::free_array4(IdxCu_DEV);
    turbotmp::free_array4(IdyCv_DEV);
    turbotmp::free_array4(CAu_DEV);
    turbotmp::free_array4(CAv_DEV);
}

/// @brief Bridge for the MOM6 CorAdv_weno subroutine.
///
/// @param bxU_HOST             The u-point iteration box of this tile
/// @param bxV_HOST             The v-point iteration box of this tile
/// @param u_HOST               Zonal velocity [L T-1 ~> m s-1]
/// @param v_HOST               Meridional velocity [L T-1 ~> m s-1]
/// @param uh_HOST              Zonal transport u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param vh_HOST              Meridional transport v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param q_HOST               Layer potential vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
/// @param abs_vort_HOST        Absolute vorticity at q-points [T-1 ~> s-1]
/// @param h_q_HOST             Thickness interpolated to q points [H ~> m or kg m-2]
/// @param Gcore_HOST           Persistent copies of the core ocean grid fields
/// @param H_subroundoff        A thickness that is so small that it can be added to a thickness of
///                             Angstrom or larger without changing it at the bit level
///                             [H ~> m or kg m-2].
/// @param CAu_HOST             Zonal acceleration due to Coriolis and momentum advection [L T-2 ~> m s-2].
/// @param CAv_HOST             Meridional acceleration due to Coriolis and momentum advection
///                             [L T-2 ~> m s-2].
/// @param weno_velocity_smooth If true, use velocity to compute the smoothness indicator for WENO
/// @param CS_HOST              Options selecting the Coriolis discretization
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
                                 const Coriolis_scheme_CS_C* CS_HOST)
{
    /// Reconstruct the iteration Boxes (Fortran 1-based -> AMReX 0-based)
    Box bxU = turbotmp::to_box(bxU_HOST);
    Box bxV = turbotmp::to_box(bxV_HOST);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto u_DEV        = turbotmp::make_array4(u_HOST);
    auto v_DEV        = turbotmp::make_array4(v_HOST);
    auto uh_DEV       = turbotmp::make_array4(uh_HOST);
    auto vh_DEV       = turbotmp::make_array4(vh_HOST);
    auto q_DEV        = turbotmp::make_array4(q_HOST);
    auto abs_vort_DEV = turbotmp::make_array4(abs_vort_HOST);
    auto h_q_DEV      = turbotmp::make_array4(h_q_HOST);
    auto IdxCu_DEV    = turbotmp::make_array4(&Gcore_HOST->IdxCu);
    auto mask2dCu_DEV = turbotmp::make_array4(&Gcore_HOST->mask2dCu);
    auto IdyCv_DEV    = turbotmp::make_array4(&Gcore_HOST->IdyCv);
    auto mask2dCv_DEV = turbotmp::make_array4(&Gcore_HOST->mask2dCv);
    auto CAu_DEV      = turbotmp::make_array4(CAu_HOST);
    auto CAv_DEV      = turbotmp::make_array4(CAv_HOST);

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(u_HOST->data,                  u_DEV);
    turbotmp::copy_FortranHost_to_array4(v_HOST->data,                  v_DEV);
    turbotmp::copy_FortranHost_to_array4(uh_HOST->data,                 uh_DEV);
    turbotmp::copy_FortranHost_to_array4(vh_HOST->data,                 vh_DEV);
    turbotmp::copy_FortranHost_to_array4(q_HOST->data,                  q_DEV);
    turbotmp::copy_FortranHost_to_array4(abs_vort_HOST->data,           abs_vort_DEV);
    turbotmp::copy_FortranHost_to_array4(h_q_HOST->data,                h_q_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdxCu.data,        IdxCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dCu.data,     mask2dCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdyCv.data,        IdyCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dCv.data,     mask2dCv_DEV);
    turbotmp::copy_FortranHost_to_array4(CAu_HOST->data,                CAu_DEV);
    turbotmp::copy_FortranHost_to_array4(CAv_HOST->data,                CAv_DEV);

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::coradv_weno(bxU, bxV,
                     u_DEV.arr, v_DEV.arr, uh_DEV.arr, vh_DEV.arr,
                     q_DEV.arr, abs_vort_DEV.arr, h_q_DEV.arr,
                     IdxCu_DEV.arr, mask2dCu_DEV.arr, IdyCv_DEV.arr, mask2dCv_DEV.arr,
                     H_subroundoff,
                     CAu_DEV.arr, CAv_DEV.arr,
                     weno_velocity_smooth,
                     *CS_HOST);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(CAu_DEV,      CAu_HOST->data);
    turbotmp::copy_array4_to_FortranHost(CAv_DEV,      CAv_HOST->data);

    /// Free everything we allocated
    turbotmp::free_array4(u_DEV);
    turbotmp::free_array4(v_DEV);
    turbotmp::free_array4(uh_DEV);
    turbotmp::free_array4(vh_DEV);
    turbotmp::free_array4(q_DEV);
    turbotmp::free_array4(abs_vort_DEV);
    turbotmp::free_array4(h_q_DEV);
    turbotmp::free_array4(IdxCu_DEV);
    turbotmp::free_array4(mask2dCu_DEV);
    turbotmp::free_array4(IdyCv_DEV);
    turbotmp::free_array4(mask2dCv_DEV);
    turbotmp::free_array4(CAu_DEV);
    turbotmp::free_array4(CAv_DEV);
}

/// @brief Bridge for the MOM6 CorAdCalc_TR subroutine.
///
/// @param bxH0_HOST           The h-point box of the whole column
/// @param bxQ0_HOST           The B-grid-index box of the whole column
/// @param u_HOST              Zonal velocity [L T-1 ~> m s-1]
/// @param v_HOST              Meridional velocity [L T-1 ~> m s-1]
/// @param h_HOST              Layer thickness [H ~> m or kg m-2]
/// @param uh_HOST             Zonal transport u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param vh_HOST             Meridional transport v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param CAu_HOST            Zonal acceleration due to Coriolis and momentum advection [L T-2 ~> m s-2].
/// @param CAv_HOST            Meridional acceleration due to Coriolis and momentum advection
///                            [L T-2 ~> m s-2].
/// @param RV_HOST             Diagnostic relative vorticity [T-1 ~> s-1]
/// @param PV_HOST             Diagnostic potential vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
/// @param CAuS_HOST           Stokes contribution to CAu [L T-2 ~> m s-2]
/// @param CAvS_HOST           Stokes contribution to CAv [L T-2 ~> m s-2]
/// @param OBC                 Open boundary control structure
/// @param AD_HOST             Copies of the acceleration diagnostics; each member may be absent
///                            (data == nullptr)
/// @param Gcore_HOST          Persistent copies of the core ocean grid fields
/// @param Gobc_HOST           Persistent copies of the ocean grid OBC mask fields
/// @param nz                  The number of layers/levels in the vertical
/// @param H_subroundoff       A thickness that is so small that it can be added to a thickness of
///                            Angstrom or larger without changing it at the bit level
///                            [H ~> m or kg m-2].
/// @param Angstrom_H          A one-Angstrom thickness in the model thickness units [H ~> m or kg m-2].
/// @param m_to_L              A constant that translates lengths in meters to the units of
///                            horizontal lengths [L m-1 ~> 1]
/// @param m_s_to_L_T          Convert lateral velocities from m s-1 to L T-1 [L s T-1 m-1 ~> 1]
/// @param CS_HOST             Options used in the Coriolis and momentum advection calculations
/// @param por_face_areaU_HOST fractional open area of U-faces [nondim]
/// @param por_face_areaV_HOST fractional open area of V-faces [nondim]
/// @param Waves               An optional pointer to Stokes drift CS
/// @param Stokes_VF           If true, include the Stokes drift
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
                                  const bool Stokes_VF)
{
    /// Reconstruct the iteration Boxes (Fortran 1-based -> AMReX 0-based)
    Box bxH0 = turbotmp::to_box(bxH0_HOST);
    Box bxQ0 = turbotmp::to_box(bxQ0_HOST);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto u_DEV              = turbotmp::make_array4(u_HOST);
    auto v_DEV              = turbotmp::make_array4(v_HOST);
    auto h_DEV              = turbotmp::make_array4(h_HOST);
    auto uh_DEV             = turbotmp::make_array4(uh_HOST);
    auto vh_DEV             = turbotmp::make_array4(vh_HOST);
    auto CAu_DEV            = turbotmp::make_array4(CAu_HOST);
    auto CAv_DEV            = turbotmp::make_array4(CAv_HOST);
    auto RV_DEV             = turbotmp::make_array4(RV_HOST);
    auto PV_DEV             = turbotmp::make_array4(PV_HOST);
    auto CAuS_DEV           = turbotmp::make_array4(CAuS_HOST);
    auto CAvS_DEV           = turbotmp::make_array4(CAvS_HOST);
    auto mask2dT_DEV        = turbotmp::make_array4(&Gcore_HOST->mask2dT);
    auto areaT_DEV          = turbotmp::make_array4(&Gcore_HOST->areaT);
    auto areaCu_DEV         = turbotmp::make_array4(&Gcore_HOST->areaCu);
    auto areaCv_DEV         = turbotmp::make_array4(&Gcore_HOST->areaCv);
    auto IareaT_DEV         = turbotmp::make_array4(&Gcore_HOST->IareaT);
    auto mask2dCu_DEV       = turbotmp::make_array4(&Gcore_HOST->mask2dCu);
    auto mask2dCv_DEV       = turbotmp::make_array4(&Gcore_HOST->mask2dCv);
    auto dyCv_DEV           = turbotmp::make_array4(&Gcore_HOST->dyCv);
    auto dxCu_DEV           = turbotmp::make_array4(&Gcore_HOST->dxCu);
    auto dy_Cu_DEV          = turbotmp::make_array4(&Gcore_HOST->dy_Cu);
    auto dx_Cv_DEV          = turbotmp::make_array4(&Gcore_HOST->dx_Cv);
    auto dyBu_DEV           = turbotmp::make_array4(&Gcore_HOST->dyBu);
    auto dxBu_DEV           = turbotmp::make_array4(&Gcore_HOST->dxBu);
    auto mask2dBu_DEV       = turbotmp::make_array4(&Gcore_HOST->mask2dBu);
    auto IareaBu_DEV        = turbotmp::make_array4(&Gcore_HOST->IareaBu);
    auto CoriolisBu_DEV     = turbotmp::make_array4(&Gcore_HOST->CoriolisBu);
    auto IdxCu_DEV          = turbotmp::make_array4(&Gcore_HOST->IdxCu);
    auto IdyCv_DEV          = turbotmp::make_array4(&Gcore_HOST->IdyCv);
    auto IdxCv_DEV          = turbotmp::make_array4(&Gcore_HOST->IdxCv);
    auto IdyCu_DEV          = turbotmp::make_array4(&Gcore_HOST->IdyCu);
    auto IdxCu_OBCmask_DEV  = turbotmp::make_array4(&Gobc_HOST->IdxCu_OBCmask);
    auto IdyCv_OBCmask_DEV  = turbotmp::make_array4(&Gobc_HOST->IdyCv_OBCmask);
    auto por_face_areaU_DEV = turbotmp::make_array4(por_face_areaU_HOST);
    auto por_face_areaV_DEV = turbotmp::make_array4(por_face_areaV_HOST);

    /// The acceleration diagnostics may be absent (data == nullptr); only allocate/copy them when present.
    const bool has_gradKEu = (AD_HOST->gradKEu.data != nullptr);
    const bool has_gradKEv = (AD_HOST->gradKEv.data != nullptr);
    const bool has_rv_x_v  = (AD_HOST->rv_x_v.data != nullptr);
    const bool has_rv_x_u  = (AD_HOST->rv_x_u.data != nullptr);
    turbotmp::A4Box gradKEu_DEV{};
    turbotmp::A4Box gradKEv_DEV{};
    turbotmp::A4Box rv_x_v_DEV {};
    turbotmp::A4Box rv_x_u_DEV {};
    if (has_gradKEu) {
        gradKEu_DEV = turbotmp::make_array4(&AD_HOST->gradKEu);
    }
    if (has_gradKEv) {
        gradKEv_DEV = turbotmp::make_array4(&AD_HOST->gradKEv);
    }
    if (has_rv_x_v) {
        rv_x_v_DEV = turbotmp::make_array4(&AD_HOST->rv_x_v);
    }
    if (has_rv_x_u) {
        rv_x_u_DEV = turbotmp::make_array4(&AD_HOST->rv_x_u);
    }

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(u_HOST->data,                        u_DEV);
    turbotmp::copy_FortranHost_to_array4(v_HOST->data,                        v_DEV);
    turbotmp::copy_FortranHost_to_array4(h_HOST->data,                        h_DEV);
    turbotmp::copy_FortranHost_to_array4(uh_HOST->data,                       uh_DEV);
    turbotmp::copy_FortranHost_to_array4(vh_HOST->data,                       vh_DEV);
    turbotmp::copy_FortranHost_to_array4(CAu_HOST->data,                      CAu_DEV);
    turbotmp::copy_FortranHost_to_array4(CAv_HOST->data,                      CAv_DEV);
    turbotmp::copy_FortranHost_to_array4(RV_HOST->data,                       RV_DEV);
    turbotmp::copy_FortranHost_to_array4(PV_HOST->data,                       PV_DEV);
    turbotmp::copy_FortranHost_to_array4(CAuS_HOST->data,                     CAuS_DEV);
    turbotmp::copy_FortranHost_to_array4(CAvS_HOST->data,                     CAvS_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dT.data,            mask2dT_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->areaT.data,              areaT_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->areaCu.data,             areaCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->areaCv.data,             areaCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IareaT.data,             IareaT_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dCu.data,           mask2dCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dCv.data,           mask2dCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->dyCv.data,               dyCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->dxCu.data,               dxCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->dy_Cu.data,              dy_Cu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->dx_Cv.data,              dx_Cv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->dyBu.data,               dyBu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->dxBu.data,               dxBu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dBu.data,           mask2dBu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IareaBu.data,            IareaBu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->CoriolisBu.data,         CoriolisBu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdxCu.data,              IdxCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdyCv.data,              IdyCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdxCv.data,              IdxCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdyCu.data,              IdyCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gobc_HOST->IdxCu_OBCmask.data,       IdxCu_OBCmask_DEV);
    turbotmp::copy_FortranHost_to_array4(Gobc_HOST->IdyCv_OBCmask.data,       IdyCv_OBCmask_DEV);
    turbotmp::copy_FortranHost_to_array4(por_face_areaU_HOST->data,           por_face_areaU_DEV);
    turbotmp::copy_FortranHost_to_array4(por_face_areaV_HOST->data,           por_face_areaV_DEV);
    if (has_gradKEu) {
        turbotmp::copy_FortranHost_to_array4(AD_HOST->gradKEu.data, gradKEu_DEV);
    }
    if (has_gradKEv) {
        turbotmp::copy_FortranHost_to_array4(AD_HOST->gradKEv.data, gradKEv_DEV);
    }
    if (has_rv_x_v) {
        turbotmp::copy_FortranHost_to_array4(AD_HOST->rv_x_v.data, rv_x_v_DEV);
    }
    if (has_rv_x_u) {
        turbotmp::copy_FortranHost_to_array4(AD_HOST->rv_x_u.data, rv_x_u_DEV);
    }

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::coradcalc_tr(bxH0, bxQ0,
                      u_DEV.arr, v_DEV.arr, h_DEV.arr, uh_DEV.arr, vh_DEV.arr,
                      CAu_DEV.arr, CAv_DEV.arr, RV_DEV.arr, PV_DEV.arr, CAuS_DEV.arr, CAvS_DEV.arr,
                      OBC,
                      gradKEu_DEV.arr, gradKEv_DEV.arr, rv_x_v_DEV.arr, rv_x_u_DEV.arr,
                      mask2dT_DEV.arr, areaT_DEV.arr, areaCu_DEV.arr, areaCv_DEV.arr, IareaT_DEV.arr,
                      mask2dCu_DEV.arr, mask2dCv_DEV.arr, dyCv_DEV.arr, dxCu_DEV.arr,
                      dy_Cu_DEV.arr, dx_Cv_DEV.arr, dyBu_DEV.arr, dxBu_DEV.arr,
                      mask2dBu_DEV.arr, IareaBu_DEV.arr, CoriolisBu_DEV.arr,
                      IdxCu_DEV.arr, IdyCv_DEV.arr, IdxCv_DEV.arr, IdyCu_DEV.arr,
                      IdxCu_OBCmask_DEV.arr, IdyCv_OBCmask_DEV.arr,
                      nz, H_subroundoff, Angstrom_H, m_to_L, m_s_to_L_T,
                      *CS_HOST,
                      por_face_areaU_DEV.arr, por_face_areaV_DEV.arr,
                      Waves, Stokes_VF);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(CAu_DEV,            CAu_HOST->data);
    turbotmp::copy_array4_to_FortranHost(CAv_DEV,            CAv_HOST->data);
    turbotmp::copy_array4_to_FortranHost(RV_DEV,             RV_HOST->data);
    turbotmp::copy_array4_to_FortranHost(PV_DEV,             PV_HOST->data);
    turbotmp::copy_array4_to_FortranHost(CAuS_DEV,           CAuS_HOST->data);
    turbotmp::copy_array4_to_FortranHost(CAvS_DEV,           CAvS_HOST->data);
    if (has_gradKEu) {
        turbotmp::copy_array4_to_FortranHost(gradKEu_DEV, AD_HOST->gradKEu.data);
    }
    if (has_gradKEv) {
        turbotmp::copy_array4_to_FortranHost(gradKEv_DEV, AD_HOST->gradKEv.data);
    }
    if (has_rv_x_v) {
        turbotmp::copy_array4_to_FortranHost(rv_x_v_DEV, AD_HOST->rv_x_v.data);
    }
    if (has_rv_x_u) {
        turbotmp::copy_array4_to_FortranHost(rv_x_u_DEV, AD_HOST->rv_x_u.data);
    }

    /// Free everything we allocated (free_array4 on a never-allocated A4Box is a safe no-op)
    turbotmp::free_array4(u_DEV);
    turbotmp::free_array4(v_DEV);
    turbotmp::free_array4(h_DEV);
    turbotmp::free_array4(uh_DEV);
    turbotmp::free_array4(vh_DEV);
    turbotmp::free_array4(CAu_DEV);
    turbotmp::free_array4(CAv_DEV);
    turbotmp::free_array4(RV_DEV);
    turbotmp::free_array4(PV_DEV);
    turbotmp::free_array4(CAuS_DEV);
    turbotmp::free_array4(CAvS_DEV);
    turbotmp::free_array4(mask2dT_DEV);
    turbotmp::free_array4(areaT_DEV);
    turbotmp::free_array4(areaCu_DEV);
    turbotmp::free_array4(areaCv_DEV);
    turbotmp::free_array4(IareaT_DEV);
    turbotmp::free_array4(mask2dCu_DEV);
    turbotmp::free_array4(mask2dCv_DEV);
    turbotmp::free_array4(dyCv_DEV);
    turbotmp::free_array4(dxCu_DEV);
    turbotmp::free_array4(dy_Cu_DEV);
    turbotmp::free_array4(dx_Cv_DEV);
    turbotmp::free_array4(dyBu_DEV);
    turbotmp::free_array4(dxBu_DEV);
    turbotmp::free_array4(mask2dBu_DEV);
    turbotmp::free_array4(IareaBu_DEV);
    turbotmp::free_array4(CoriolisBu_DEV);
    turbotmp::free_array4(IdxCu_DEV);
    turbotmp::free_array4(IdyCv_DEV);
    turbotmp::free_array4(IdxCv_DEV);
    turbotmp::free_array4(IdyCu_DEV);
    turbotmp::free_array4(IdxCu_OBCmask_DEV);
    turbotmp::free_array4(IdyCv_OBCmask_DEV);
    turbotmp::free_array4(por_face_areaU_DEV);
    turbotmp::free_array4(por_face_areaV_DEV);
    turbotmp::free_array4(gradKEu_DEV);
    turbotmp::free_array4(gradKEv_DEV);
    turbotmp::free_array4(rv_x_v_DEV);
    turbotmp::free_array4(rv_x_u_DEV);
}
