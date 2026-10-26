// SKILLS: 0.3.3
/// @file turbotmp_mom_hor_visc_bridge.cpp
/// @brief Marshalling bridges between the MOM6 Fortran shims in MOM_hor_visc
///        (mode = AMREX) and the AMReX kernels in mom_hor_visc.cpp.

#include <AMReX_Arena.H>
#include <AMReX_Gpu.H>

#include "mom_hor_visc.hpp"
#include "turbotmp_helper.hpp"
#include "turbotmp_mom_hor_visc_bridge.h"

using namespace amrex;

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
                                 const bool zero_land)
{
    /// Reconstruct the iteration Box (Fortran 1-based -> AMReX 0-based)
    Box bxH0 = turbotmp::to_box(bxH0_HOST);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto mask2dT_DEV = turbotmp::make_array4(&Gcore_HOST->mask2dT);
    auto field_h_DEV = turbotmp::make_array4(field_h_HOST);

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dT.data, mask2dT_DEV);
    turbotmp::copy_FortranHost_to_array4(field_h_HOST->data,       field_h_DEV);

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::smooth_x9_h(bxH0, mask2dT_DEV.arr, field_h_DEV.arr, zero_land);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(field_h_DEV, field_h_HOST->data);

    /// Free everything we allocated
    turbotmp::free_array4(mask2dT_DEV);
    turbotmp::free_array4(field_h_DEV);
}

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
                                  const bool zero_land)
{
    /// Reconstruct the iteration Boxes (Fortran 1-based -> AMReX 0-based)
    Box bxU0 = turbotmp::to_box(bxU0_HOST);
    Box bxV0 = turbotmp::to_box(bxV0_HOST);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto mask2dCu_DEV = turbotmp::make_array4(&Gcore_HOST->mask2dCu);
    auto mask2dCv_DEV = turbotmp::make_array4(&Gcore_HOST->mask2dCv);
    auto field_u_DEV  = turbotmp::make_array4(field_u_HOST);
    auto field_v_DEV  = turbotmp::make_array4(field_v_HOST);

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dCu.data, mask2dCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dCv.data, mask2dCv_DEV);
    turbotmp::copy_FortranHost_to_array4(field_u_HOST->data,        field_u_DEV);
    turbotmp::copy_FortranHost_to_array4(field_v_HOST->data,        field_v_DEV);

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::smooth_x9_uv(bxU0, bxV0,
                      mask2dCu_DEV.arr, mask2dCv_DEV.arr,
                      field_u_DEV.arr, field_v_DEV.arr,
                      zero_land);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(field_u_DEV, field_u_HOST->data);
    turbotmp::copy_array4_to_FortranHost(field_v_DEV, field_v_HOST->data);

    /// Free everything we allocated
    turbotmp::free_array4(mask2dCu_DEV);
    turbotmp::free_array4(mask2dCv_DEV);
    turbotmp::free_array4(field_u_DEV);
    turbotmp::free_array4(field_v_DEV);
}

/// @brief Bridge for the MOM6 smooth_GME_h_pass subroutine.
///
/// @param Gcore_HOST      Persistent copies of the core ocean grid fields
/// @param bxPass_HOST     The h-point range of this smoothing pass, the computational
///                        domain extended by the valid halo points
/// @param GME_flux_h_HOST GME diffusive flux at h points [L2 T-2 ~> m2 s-2]
void turbotmp_smooth_gme_h_pass_bridge(const grid_core_C* Gcore_HOST,
                                       const Box_C* bxPass_HOST,
                                       RealArray_C* GME_flux_h_HOST)
{
    /// Reconstruct the iteration Box (Fortran 1-based -> AMReX 0-based)
    Box bxPass = turbotmp::to_box(bxPass_HOST);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto mask2dT_DEV    = turbotmp::make_array4(&Gcore_HOST->mask2dT);
    auto GME_flux_h_DEV = turbotmp::make_array4(GME_flux_h_HOST);

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dT.data, mask2dT_DEV);
    turbotmp::copy_FortranHost_to_array4(GME_flux_h_HOST->data,    GME_flux_h_DEV);

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::smooth_gme_h_pass(bxPass, mask2dT_DEV.arr, GME_flux_h_DEV.arr);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(GME_flux_h_DEV, GME_flux_h_HOST->data);

    /// Free everything we allocated
    turbotmp::free_array4(mask2dT_DEV);
    turbotmp::free_array4(GME_flux_h_DEV);
}

/// @brief Bridge for the MOM6 smooth_GME_q_pass subroutine.
///
/// @param Gcore_HOST      Persistent copies of the core ocean grid fields
/// @param bxPass_HOST     The q-point range of this smoothing pass, the computational
///                        domain extended by the valid halo points
/// @param GME_flux_q_HOST GME diffusive flux at q points [L2 T-2 ~> m2 s-2]
void turbotmp_smooth_gme_q_pass_bridge(const grid_core_C* Gcore_HOST,
                                       const Box_C* bxPass_HOST,
                                       RealArray_C* GME_flux_q_HOST)
{
    /// Reconstruct the iteration Box (Fortran 1-based -> AMReX 0-based)
    Box bxPass = turbotmp::to_box(bxPass_HOST);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto mask2dBu_DEV   = turbotmp::make_array4(&Gcore_HOST->mask2dBu);
    auto GME_flux_q_DEV = turbotmp::make_array4(GME_flux_q_HOST);

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dBu.data, mask2dBu_DEV);
    turbotmp::copy_FortranHost_to_array4(GME_flux_q_HOST->data,     GME_flux_q_DEV);

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::smooth_gme_q_pass(bxPass, mask2dBu_DEV.arr, GME_flux_q_DEV.arr);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(GME_flux_q_DEV, GME_flux_q_HOST->data);

    /// Free everything we allocated
    turbotmp::free_array4(mask2dBu_DEV);
    turbotmp::free_array4(GME_flux_q_DEV);
}

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
                                            RealArray_C* BS_coeff_h_HOST)
{
    (void)nkblock;

    /// Reconstruct the iteration Box (Fortran 1-based -> AMReX 0-based)
    Box bxQ = turbotmp::to_box(bxQ_HOST);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto visc_limit_h_flag_DEV = turbotmp::make_array4(visc_limit_h_flag_HOST);
    auto sh_xx_DEV             = turbotmp::make_array4(sh_xx_HOST);
    auto str_xx_DEV            = turbotmp::make_array4(str_xx_HOST);
    auto BS_coeff_h_DEV        = turbotmp::make_array4(BS_coeff_h_HOST);

    /// MEKE%Ku and VarMix%BS_struct may be unallocated (data == nullptr); only allocate/copy them when allocated.
    const bool MEKE_Ku_allocated   = (MEKE_HOST->Ku.data != nullptr);
    const bool BS_struct_allocated = (VarMix_HOST->BS_struct.data != nullptr);
    turbotmp::A4Box MEKE_Ku_DEV{};
    turbotmp::A4Box BS_struct_DEV{};
    if (MEKE_Ku_allocated)   MEKE_Ku_DEV   = turbotmp::make_array4(&MEKE_HOST->Ku);
    if (BS_struct_allocated) BS_struct_DEV = turbotmp::make_array4(&VarMix_HOST->BS_struct);

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(visc_limit_h_flag_HOST->data, visc_limit_h_flag_DEV);
    turbotmp::copy_FortranHost_to_array4(sh_xx_HOST->data,             sh_xx_DEV);
    turbotmp::copy_FortranHost_to_array4(str_xx_HOST->data,            str_xx_DEV);
    turbotmp::copy_FortranHost_to_array4(BS_coeff_h_HOST->data,        BS_coeff_h_DEV);
    if (MEKE_Ku_allocated)   turbotmp::copy_FortranHost_to_array4(MEKE_HOST->Ku.data,          MEKE_Ku_DEV);
    if (BS_struct_allocated) turbotmp::copy_FortranHost_to_array4(VarMix_HOST->BS_struct.data, BS_struct_DEV);

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::hor_visc_backscatter_h(bxQ, *CS_HOST,
                                MEKE_Ku_DEV.arr, BS_struct_DEV.arr,
                                use_kh_struct, kstart,
                                visc_limit_h_flag_DEV.arr, sh_xx_DEV.arr,
                                str_xx_DEV.arr, BS_coeff_h_DEV.arr);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(str_xx_DEV,     str_xx_HOST->data);
    turbotmp::copy_array4_to_FortranHost(BS_coeff_h_DEV, BS_coeff_h_HOST->data);

    /// Free everything we allocated (freeing a never-allocated container is a no-op)
    turbotmp::free_array4(visc_limit_h_flag_DEV);
    turbotmp::free_array4(sh_xx_DEV);
    turbotmp::free_array4(str_xx_DEV);
    turbotmp::free_array4(BS_coeff_h_DEV);
    turbotmp::free_array4(MEKE_Ku_DEV);
    turbotmp::free_array4(BS_struct_DEV);
}

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
                                            RealArray_C* BS_coeff_q_HOST)
{
    (void)nkblock;

    /// Reconstruct the iteration Box (Fortran 1-based -> AMReX 0-based)
    Box bxH = turbotmp::to_box(bxH_HOST);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto visc_limit_q_flag_DEV = turbotmp::make_array4(visc_limit_q_flag_HOST);
    auto sh_xy_DEV             = turbotmp::make_array4(sh_xy_HOST);
    auto str_xy_DEV            = turbotmp::make_array4(str_xy_HOST);
    auto BS_coeff_q_DEV        = turbotmp::make_array4(BS_coeff_q_HOST);

    /// MEKE%Ku and VarMix%BS_struct may be unallocated (data == nullptr); only allocate/copy them when allocated.
    const bool MEKE_Ku_allocated   = (MEKE_HOST->Ku.data != nullptr);
    const bool BS_struct_allocated = (VarMix_HOST->BS_struct.data != nullptr);
    turbotmp::A4Box MEKE_Ku_DEV{};
    turbotmp::A4Box BS_struct_DEV{};
    if (MEKE_Ku_allocated)   MEKE_Ku_DEV   = turbotmp::make_array4(&MEKE_HOST->Ku);
    if (BS_struct_allocated) BS_struct_DEV = turbotmp::make_array4(&VarMix_HOST->BS_struct);

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(visc_limit_q_flag_HOST->data, visc_limit_q_flag_DEV);
    turbotmp::copy_FortranHost_to_array4(sh_xy_HOST->data,             sh_xy_DEV);
    turbotmp::copy_FortranHost_to_array4(str_xy_HOST->data,            str_xy_DEV);
    turbotmp::copy_FortranHost_to_array4(BS_coeff_q_HOST->data,        BS_coeff_q_DEV);
    if (MEKE_Ku_allocated)   turbotmp::copy_FortranHost_to_array4(MEKE_HOST->Ku.data,          MEKE_Ku_DEV);
    if (BS_struct_allocated) turbotmp::copy_FortranHost_to_array4(VarMix_HOST->BS_struct.data, BS_struct_DEV);

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::hor_visc_backscatter_q(bxH, *CS_HOST,
                                MEKE_Ku_DEV.arr, BS_struct_DEV.arr,
                                use_kh_struct, kstart,
                                visc_limit_q_flag_DEV.arr, sh_xy_DEV.arr,
                                str_xy_DEV.arr, BS_coeff_q_DEV.arr);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(str_xy_DEV,     str_xy_HOST->data);
    turbotmp::copy_array4_to_FortranHost(BS_coeff_q_DEV, BS_coeff_q_HOST->data);

    /// Free everything we allocated (freeing a never-allocated container is a no-op)
    turbotmp::free_array4(visc_limit_q_flag_DEV);
    turbotmp::free_array4(sh_xy_DEV);
    turbotmp::free_array4(str_xy_DEV);
    turbotmp::free_array4(BS_coeff_q_DEV);
    turbotmp::free_array4(MEKE_Ku_DEV);
    turbotmp::free_array4(BS_struct_DEV);
}

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
                                     RealArray_C* GME_effic_q_HOST)
{
    /// Reconstruct the iteration Boxes (Fortran 1-based -> AMReX 0-based)
    Box bxH0 = turbotmp::to_box(bxH0_HOST);
    Box bxQ0 = turbotmp::to_box(bxQ0_HOST);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto IdyCu_DEV       = turbotmp::make_array4(&Gcore_HOST->IdyCu);
    auto IdxCv_DEV       = turbotmp::make_array4(&Gcore_HOST->IdxCv);
    auto IdyCv_DEV       = turbotmp::make_array4(&Gcore_HOST->IdyCv);
    auto IdxCu_DEV       = turbotmp::make_array4(&Gcore_HOST->IdxCu);
    auto mask2dBu_DEV    = turbotmp::make_array4(&Gcore_HOST->mask2dBu);
    auto mask2dCu_DEV    = turbotmp::make_array4(&Gcore_HOST->mask2dCu);
    auto mask2dCv_DEV    = turbotmp::make_array4(&Gcore_HOST->mask2dCv);
    auto mask2dT_DEV     = turbotmp::make_array4(&Gcore_HOST->mask2dT);
    auto h_DEV           = turbotmp::make_array4(h_HOST);
    auto ubtav_DEV       = turbotmp::make_array4(ubtav_HOST);
    auto vbtav_DEV       = turbotmp::make_array4(vbtav_HOST);
    auto dudx_bt_DEV     = turbotmp::make_array4(dudx_bt_HOST);
    auto dvdy_bt_DEV     = turbotmp::make_array4(dvdy_bt_HOST);
    auto dvdx_bt_DEV     = turbotmp::make_array4(dvdx_bt_HOST);
    auto dudy_bt_DEV     = turbotmp::make_array4(dudy_bt_HOST);
    auto sh_xx_bt_DEV    = turbotmp::make_array4(sh_xx_bt_HOST);
    auto sh_xy_bt_DEV    = turbotmp::make_array4(sh_xy_bt_HOST);
    auto GME_effic_h_DEV = turbotmp::make_array4(GME_effic_h_HOST);
    auto GME_effic_q_DEV = turbotmp::make_array4(GME_effic_q_HOST);

    /// The metric arrays may be unallocated (data == nullptr); only allocate/copy them when allocated.
    const bool dx_dyT_allocated  = (metrics_HOST->dx_dyT.data != nullptr);
    const bool dy_dxT_allocated  = (metrics_HOST->dy_dxT.data != nullptr);
    const bool dx_dyBu_allocated = (metrics_HOST->dx_dyBu.data != nullptr);
    const bool dy_dxBu_allocated = (metrics_HOST->dy_dxBu.data != nullptr);
    turbotmp::A4Box dx_dyT_DEV{};
    turbotmp::A4Box dy_dxT_DEV{};
    turbotmp::A4Box dx_dyBu_DEV{};
    turbotmp::A4Box dy_dxBu_DEV{};
    if (dx_dyT_allocated)  dx_dyT_DEV  = turbotmp::make_array4(&metrics_HOST->dx_dyT);
    if (dy_dxT_allocated)  dy_dxT_DEV  = turbotmp::make_array4(&metrics_HOST->dy_dxT);
    if (dx_dyBu_allocated) dx_dyBu_DEV = turbotmp::make_array4(&metrics_HOST->dx_dyBu);
    if (dy_dxBu_allocated) dy_dxBu_DEV = turbotmp::make_array4(&metrics_HOST->dy_dxBu);

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdyCu.data,    IdyCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdxCv.data,    IdxCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdyCv.data,    IdyCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdxCu.data,    IdxCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dBu.data, mask2dBu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dCu.data, mask2dCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dCv.data, mask2dCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dT.data,  mask2dT_DEV);
    turbotmp::copy_FortranHost_to_array4(h_HOST->data,              h_DEV);
    turbotmp::copy_FortranHost_to_array4(ubtav_HOST->data,          ubtav_DEV);
    turbotmp::copy_FortranHost_to_array4(vbtav_HOST->data,          vbtav_DEV);
    turbotmp::copy_FortranHost_to_array4(dudx_bt_HOST->data,        dudx_bt_DEV);
    turbotmp::copy_FortranHost_to_array4(dvdy_bt_HOST->data,        dvdy_bt_DEV);
    turbotmp::copy_FortranHost_to_array4(dvdx_bt_HOST->data,        dvdx_bt_DEV);
    turbotmp::copy_FortranHost_to_array4(dudy_bt_HOST->data,        dudy_bt_DEV);
    turbotmp::copy_FortranHost_to_array4(sh_xx_bt_HOST->data,       sh_xx_bt_DEV);
    turbotmp::copy_FortranHost_to_array4(sh_xy_bt_HOST->data,       sh_xy_bt_DEV);
    turbotmp::copy_FortranHost_to_array4(GME_effic_h_HOST->data,    GME_effic_h_DEV);
    turbotmp::copy_FortranHost_to_array4(GME_effic_q_HOST->data,    GME_effic_q_DEV);
    if (dx_dyT_allocated)  turbotmp::copy_FortranHost_to_array4(metrics_HOST->dx_dyT.data, dx_dyT_DEV);
    if (dy_dxT_allocated)  turbotmp::copy_FortranHost_to_array4(metrics_HOST->dy_dxT.data, dy_dxT_DEV);
    if (dx_dyBu_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->dx_dyBu.data, dx_dyBu_DEV);
    if (dy_dxBu_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->dy_dxBu.data, dy_dxBu_DEV);

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::hor_visc_gme_bt(bxH0, bxQ0, *CS_HOST,
                         IdyCu_DEV.arr, IdxCv_DEV.arr, IdyCv_DEV.arr, IdxCu_DEV.arr,
                         mask2dBu_DEV.arr, mask2dCu_DEV.arr, mask2dCv_DEV.arr, mask2dT_DEV.arr,
                         dx_dyT_DEV.arr, dy_dxT_DEV.arr, dx_dyBu_DEV.arr, dy_dxBu_DEV.arr,
                         h_DEV.arr, ubtav_DEV.arr, vbtav_DEV.arr,
                         dudx_bt_DEV.arr, dvdy_bt_DEV.arr, dvdx_bt_DEV.arr, dudy_bt_DEV.arr,
                         sh_xx_bt_DEV.arr, sh_xy_bt_DEV.arr, GME_effic_h_DEV.arr, GME_effic_q_DEV.arr);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(h_DEV,           h_HOST->data);
    turbotmp::copy_array4_to_FortranHost(dudx_bt_DEV,     dudx_bt_HOST->data);
    turbotmp::copy_array4_to_FortranHost(dvdy_bt_DEV,     dvdy_bt_HOST->data);
    turbotmp::copy_array4_to_FortranHost(dvdx_bt_DEV,     dvdx_bt_HOST->data);
    turbotmp::copy_array4_to_FortranHost(dudy_bt_DEV,     dudy_bt_HOST->data);
    turbotmp::copy_array4_to_FortranHost(sh_xx_bt_DEV,    sh_xx_bt_HOST->data);
    turbotmp::copy_array4_to_FortranHost(sh_xy_bt_DEV,    sh_xy_bt_HOST->data);
    turbotmp::copy_array4_to_FortranHost(GME_effic_h_DEV, GME_effic_h_HOST->data);
    turbotmp::copy_array4_to_FortranHost(GME_effic_q_DEV, GME_effic_q_HOST->data);

    /// Free everything we allocated (freeing a never-allocated container is a no-op)
    turbotmp::free_array4(IdyCu_DEV);
    turbotmp::free_array4(IdxCv_DEV);
    turbotmp::free_array4(IdyCv_DEV);
    turbotmp::free_array4(IdxCu_DEV);
    turbotmp::free_array4(mask2dBu_DEV);
    turbotmp::free_array4(mask2dCu_DEV);
    turbotmp::free_array4(mask2dCv_DEV);
    turbotmp::free_array4(mask2dT_DEV);
    turbotmp::free_array4(h_DEV);
    turbotmp::free_array4(ubtav_DEV);
    turbotmp::free_array4(vbtav_DEV);
    turbotmp::free_array4(dudx_bt_DEV);
    turbotmp::free_array4(dvdy_bt_DEV);
    turbotmp::free_array4(dvdx_bt_DEV);
    turbotmp::free_array4(dudy_bt_DEV);
    turbotmp::free_array4(sh_xx_bt_DEV);
    turbotmp::free_array4(sh_xy_bt_DEV);
    turbotmp::free_array4(GME_effic_h_DEV);
    turbotmp::free_array4(GME_effic_q_DEV);
    turbotmp::free_array4(dx_dyT_DEV);
    turbotmp::free_array4(dy_dxT_DEV);
    turbotmp::free_array4(dx_dyBu_DEV);
    turbotmp::free_array4(dy_dxBu_DEV);
}

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
                                           RealArray_C* div_xx_dx_HOST)
{
    (void)nkblock;

    /// Reconstruct the iteration Boxes (Fortran 1-based -> AMReX 0-based)
    Box bxH    = turbotmp::to_box(bxH_HOST);
    Box bxKh   = turbotmp::to_box(bxKh_HOST);
    Box bxKh_i = turbotmp::to_box(bxKh_i_HOST);
    Box bxKh_j = turbotmp::to_box(bxKh_j_HOST);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto dyBu_DEV               = turbotmp::make_array4(&Gcore_HOST->dyBu);
    auto IdxBu_DEV              = turbotmp::make_array4(&Gcore_HOST->IdxBu);
    auto IdyCu_DEV              = turbotmp::make_array4(&Gcore_HOST->IdyCu);
    auto dxBu_DEV               = turbotmp::make_array4(&Gcore_HOST->dxBu);
    auto IdyBu_DEV              = turbotmp::make_array4(&Gcore_HOST->IdyBu);
    auto IdxCv_DEV              = turbotmp::make_array4(&Gcore_HOST->IdxCv);
    auto IdyCv_DEV              = turbotmp::make_array4(&Gcore_HOST->IdyCv);
    auto IdxCu_DEV              = turbotmp::make_array4(&Gcore_HOST->IdxCu);
    auto dF_dx_DEV              = turbotmp::make_array4(&Gcore_HOST->df_dx);
    auto dF_dy_DEV              = turbotmp::make_array4(&Gcore_HOST->df_dy);
    auto dudx_DEV               = turbotmp::make_array4(dudx_HOST);
    auto dvdy_DEV               = turbotmp::make_array4(dvdy_HOST);
    auto vort_xy_DEV            = turbotmp::make_array4(vort_xy_HOST);
    auto vort_xy_smooth_DEV     = turbotmp::make_array4(vort_xy_smooth_HOST);
    auto grad_vort_mag_h_2d_DEV = turbotmp::make_array4(grad_vort_mag_h_2d_HOST);
    auto grad_div_mag_h_DEV     = turbotmp::make_array4(grad_div_mag_h_HOST);
    auto grad_vort_mag_q_2d_DEV = turbotmp::make_array4(grad_vort_mag_q_2d_HOST);
    auto grad_div_mag_q_DEV     = turbotmp::make_array4(grad_div_mag_q_HOST);
    auto Del2vort_q_DEV         = turbotmp::make_array4(Del2vort_q_HOST);
    auto vort_xy_dx_DEV         = turbotmp::make_array4(vort_xy_dx_HOST);
    auto vort_xy_dx_smooth_DEV  = turbotmp::make_array4(vort_xy_dx_smooth_HOST);
    auto div_xx_dy_DEV          = turbotmp::make_array4(div_xx_dy_HOST);
    auto vort_xy_dy_DEV         = turbotmp::make_array4(vort_xy_dy_HOST);
    auto vort_xy_dy_smooth_DEV  = turbotmp::make_array4(vort_xy_dy_smooth_HOST);
    auto div_xx_dx_DEV          = turbotmp::make_array4(div_xx_dx_HOST);

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->dyBu.data,         dyBu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdxBu.data,        IdxBu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdyCu.data,        IdyCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->dxBu.data,         dxBu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdyBu.data,        IdyBu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdxCv.data,        IdxCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdyCv.data,        IdyCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdxCu.data,        IdxCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->df_dx.data,        dF_dx_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->df_dy.data,        dF_dy_DEV);
    turbotmp::copy_FortranHost_to_array4(dudx_HOST->data,               dudx_DEV);
    turbotmp::copy_FortranHost_to_array4(dvdy_HOST->data,               dvdy_DEV);
    turbotmp::copy_FortranHost_to_array4(vort_xy_HOST->data,            vort_xy_DEV);
    turbotmp::copy_FortranHost_to_array4(vort_xy_smooth_HOST->data,     vort_xy_smooth_DEV);
    turbotmp::copy_FortranHost_to_array4(grad_vort_mag_h_2d_HOST->data, grad_vort_mag_h_2d_DEV);
    turbotmp::copy_FortranHost_to_array4(grad_div_mag_h_HOST->data,     grad_div_mag_h_DEV);
    turbotmp::copy_FortranHost_to_array4(grad_vort_mag_q_2d_HOST->data, grad_vort_mag_q_2d_DEV);
    turbotmp::copy_FortranHost_to_array4(grad_div_mag_q_HOST->data,     grad_div_mag_q_DEV);
    turbotmp::copy_FortranHost_to_array4(Del2vort_q_HOST->data,         Del2vort_q_DEV);
    turbotmp::copy_FortranHost_to_array4(vort_xy_dx_HOST->data,         vort_xy_dx_DEV);
    turbotmp::copy_FortranHost_to_array4(vort_xy_dx_smooth_HOST->data,  vort_xy_dx_smooth_DEV);
    turbotmp::copy_FortranHost_to_array4(div_xx_dy_HOST->data,          div_xx_dy_DEV);
    turbotmp::copy_FortranHost_to_array4(vort_xy_dy_HOST->data,         vort_xy_dy_DEV);
    turbotmp::copy_FortranHost_to_array4(vort_xy_dy_smooth_HOST->data,  vort_xy_dy_smooth_DEV);
    turbotmp::copy_FortranHost_to_array4(div_xx_dx_HOST->data,          div_xx_dx_DEV);

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::hor_visc_leith_grad_a(bxH, bxKh, bxKh_i, bxKh_j, *CS_HOST,
                               dyBu_DEV.arr, IdxBu_DEV.arr, IdyCu_DEV.arr, dxBu_DEV.arr, IdyBu_DEV.arr,
                               IdxCv_DEV.arr, IdyCv_DEV.arr, IdxCu_DEV.arr, dF_dx_DEV.arr, dF_dy_DEV.arr,
                               dudx_DEV.arr, dvdy_DEV.arr, vort_xy_DEV.arr, vort_xy_smooth_DEV.arr,
                               grad_vort_mag_h_2d_DEV.arr, grad_div_mag_h_DEV.arr, grad_vort_mag_q_2d_DEV.arr, grad_div_mag_q_DEV.arr,
                               Del2vort_q_DEV.arr, vort_xy_dx_DEV.arr, vort_xy_dx_smooth_DEV.arr, div_xx_dy_DEV.arr,
                               vort_xy_dy_DEV.arr, vort_xy_dy_smooth_DEV.arr, div_xx_dx_DEV.arr);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(grad_vort_mag_h_2d_DEV, grad_vort_mag_h_2d_HOST->data);
    turbotmp::copy_array4_to_FortranHost(grad_div_mag_h_DEV,     grad_div_mag_h_HOST->data);
    turbotmp::copy_array4_to_FortranHost(grad_vort_mag_q_2d_DEV, grad_vort_mag_q_2d_HOST->data);
    turbotmp::copy_array4_to_FortranHost(grad_div_mag_q_DEV,     grad_div_mag_q_HOST->data);
    turbotmp::copy_array4_to_FortranHost(Del2vort_q_DEV,         Del2vort_q_HOST->data);
    turbotmp::copy_array4_to_FortranHost(vort_xy_dx_DEV,         vort_xy_dx_HOST->data);
    turbotmp::copy_array4_to_FortranHost(vort_xy_dx_smooth_DEV,  vort_xy_dx_smooth_HOST->data);
    turbotmp::copy_array4_to_FortranHost(div_xx_dy_DEV,          div_xx_dy_HOST->data);
    turbotmp::copy_array4_to_FortranHost(vort_xy_dy_DEV,         vort_xy_dy_HOST->data);
    turbotmp::copy_array4_to_FortranHost(vort_xy_dy_smooth_DEV,  vort_xy_dy_smooth_HOST->data);
    turbotmp::copy_array4_to_FortranHost(div_xx_dx_DEV,          div_xx_dx_HOST->data);

    /// Free everything we allocated
    turbotmp::free_array4(dyBu_DEV);
    turbotmp::free_array4(IdxBu_DEV);
    turbotmp::free_array4(IdyCu_DEV);
    turbotmp::free_array4(dxBu_DEV);
    turbotmp::free_array4(IdyBu_DEV);
    turbotmp::free_array4(IdxCv_DEV);
    turbotmp::free_array4(IdyCv_DEV);
    turbotmp::free_array4(IdxCu_DEV);
    turbotmp::free_array4(dF_dx_DEV);
    turbotmp::free_array4(dF_dy_DEV);
    turbotmp::free_array4(dudx_DEV);
    turbotmp::free_array4(dvdy_DEV);
    turbotmp::free_array4(vort_xy_DEV);
    turbotmp::free_array4(vort_xy_smooth_DEV);
    turbotmp::free_array4(grad_vort_mag_h_2d_DEV);
    turbotmp::free_array4(grad_div_mag_h_DEV);
    turbotmp::free_array4(grad_vort_mag_q_2d_DEV);
    turbotmp::free_array4(grad_div_mag_q_DEV);
    turbotmp::free_array4(Del2vort_q_DEV);
    turbotmp::free_array4(vort_xy_dx_DEV);
    turbotmp::free_array4(vort_xy_dx_smooth_DEV);
    turbotmp::free_array4(div_xx_dy_DEV);
    turbotmp::free_array4(vort_xy_dy_DEV);
    turbotmp::free_array4(vort_xy_dy_smooth_DEV);
    turbotmp::free_array4(div_xx_dx_DEV);
}

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
                                           RealArray_C* vort_xy_dy_smooth_HOST)
{
    (void)nkblock;

    /// Reconstruct the iteration Boxes (Fortran 1-based -> AMReX 0-based)
    Box bxH  = turbotmp::to_box(bxH_HOST);
    Box bxKh = turbotmp::to_box(bxKh_HOST);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto grad_vort_mag_h_DEV      = turbotmp::make_array4(grad_vort_mag_h_HOST);
    auto grad_vort_mag_q_DEV      = turbotmp::make_array4(grad_vort_mag_q_HOST);
    auto vert_vort_mag_smooth_DEV = turbotmp::make_array4(vert_vort_mag_smooth_HOST);
    auto vort_xy_dx_DEV           = turbotmp::make_array4(vort_xy_dx_HOST);
    auto vort_xy_dx_smooth_DEV    = turbotmp::make_array4(vort_xy_dx_smooth_HOST);
    auto vort_xy_dy_DEV           = turbotmp::make_array4(vort_xy_dy_HOST);
    auto vort_xy_dy_smooth_DEV    = turbotmp::make_array4(vort_xy_dy_smooth_HOST);

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(grad_vort_mag_h_HOST->data,       grad_vort_mag_h_DEV);
    turbotmp::copy_FortranHost_to_array4(grad_vort_mag_q_HOST->data,       grad_vort_mag_q_DEV);
    turbotmp::copy_FortranHost_to_array4(vert_vort_mag_smooth_HOST->data,  vert_vort_mag_smooth_DEV);
    turbotmp::copy_FortranHost_to_array4(vort_xy_dx_HOST->data,            vort_xy_dx_DEV);
    turbotmp::copy_FortranHost_to_array4(vort_xy_dx_smooth_HOST->data,     vort_xy_dx_smooth_DEV);
    turbotmp::copy_FortranHost_to_array4(vort_xy_dy_HOST->data,            vort_xy_dy_DEV);
    turbotmp::copy_FortranHost_to_array4(vort_xy_dy_smooth_HOST->data,     vort_xy_dy_smooth_DEV);

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::hor_visc_leith_grad_b(bxH, bxKh, *CS_HOST,
                               grad_vort_mag_h_DEV.arr, grad_vort_mag_q_DEV.arr, vert_vort_mag_smooth_DEV.arr,
                               vort_xy_dx_DEV.arr, vort_xy_dx_smooth_DEV.arr, vort_xy_dy_DEV.arr, vort_xy_dy_smooth_DEV.arr);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(grad_vort_mag_h_DEV,      grad_vort_mag_h_HOST->data);
    turbotmp::copy_array4_to_FortranHost(grad_vort_mag_q_DEV,      grad_vort_mag_q_HOST->data);
    turbotmp::copy_array4_to_FortranHost(vert_vort_mag_smooth_DEV, vert_vort_mag_smooth_HOST->data);
    turbotmp::copy_array4_to_FortranHost(vort_xy_dx_DEV,           vort_xy_dx_HOST->data);
    turbotmp::copy_array4_to_FortranHost(vort_xy_dx_smooth_DEV,    vort_xy_dx_smooth_HOST->data);
    turbotmp::copy_array4_to_FortranHost(vort_xy_dy_DEV,           vort_xy_dy_HOST->data);
    turbotmp::copy_array4_to_FortranHost(vort_xy_dy_smooth_DEV,    vort_xy_dy_smooth_HOST->data);

    /// Free everything we allocated
    turbotmp::free_array4(grad_vort_mag_h_DEV);
    turbotmp::free_array4(grad_vort_mag_q_DEV);
    turbotmp::free_array4(vert_vort_mag_smooth_DEV);
    turbotmp::free_array4(vort_xy_dx_DEV);
    turbotmp::free_array4(vort_xy_dx_smooth_DEV);
    turbotmp::free_array4(vort_xy_dy_DEV);
    turbotmp::free_array4(vort_xy_dy_smooth_DEV);
}

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
                                       RealArray_C* m_leithy_HOST)
{
    (void)nkblock;

    /// Reconstruct the iteration Box (Fortran 1-based -> AMReX 0-based)
    Box bxKh = turbotmp::to_box(bxKh_HOST);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto mask2dBu_DEV         = turbotmp::make_array4(&Gcore_HOST->mask2dBu);
    auto Del2vort_q_DEV       = turbotmp::make_array4(Del2vort_q_HOST);
    auto vert_vort_mag_DEV    = turbotmp::make_array4(vert_vort_mag_HOST);
    auto vort_xy_smooth_DEV   = turbotmp::make_array4(vort_xy_smooth_HOST);
    auto m_leithy_DEV         = turbotmp::make_array4(m_leithy_HOST);

    /// The metric arrays may be unallocated (data == nullptr); only allocate/copy them when allocated.
    const bool Ah_bg_xx_allocated         = (metrics_HOST->Ah_bg_xx.data != nullptr);
    const bool m_const_leithy_allocated   = (metrics_HOST->m_const_leithy.data != nullptr);
    const bool m_leithy_max_allocated     = (metrics_HOST->m_leithy_max.data != nullptr);
    const bool Biharm6_const_xx_allocated = (metrics_HOST->Biharm6_const_xx.data != nullptr);
    turbotmp::A4Box Ah_bg_xx_DEV{};
    turbotmp::A4Box m_const_leithy_DEV{};
    turbotmp::A4Box m_leithy_max_DEV{};
    turbotmp::A4Box Biharm6_const_xx_DEV{};
    if (Ah_bg_xx_allocated)         Ah_bg_xx_DEV         = turbotmp::make_array4(&metrics_HOST->Ah_bg_xx);
    if (m_const_leithy_allocated)   m_const_leithy_DEV   = turbotmp::make_array4(&metrics_HOST->m_const_leithy);
    if (m_leithy_max_allocated)     m_leithy_max_DEV     = turbotmp::make_array4(&metrics_HOST->m_leithy_max);
    if (Biharm6_const_xx_allocated) Biharm6_const_xx_DEV = turbotmp::make_array4(&metrics_HOST->Biharm6_const_xx);

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dBu.data, mask2dBu_DEV);
    turbotmp::copy_FortranHost_to_array4(Del2vort_q_HOST->data,     Del2vort_q_DEV);
    turbotmp::copy_FortranHost_to_array4(vert_vort_mag_HOST->data,  vert_vort_mag_DEV);
    turbotmp::copy_FortranHost_to_array4(vort_xy_smooth_HOST->data, vort_xy_smooth_DEV);
    turbotmp::copy_FortranHost_to_array4(m_leithy_HOST->data,       m_leithy_DEV);
    if (Ah_bg_xx_allocated)         turbotmp::copy_FortranHost_to_array4(metrics_HOST->Ah_bg_xx.data, Ah_bg_xx_DEV);
    if (m_const_leithy_allocated)   turbotmp::copy_FortranHost_to_array4(metrics_HOST->m_const_leithy.data, m_const_leithy_DEV);
    if (m_leithy_max_allocated)     turbotmp::copy_FortranHost_to_array4(metrics_HOST->m_leithy_max.data, m_leithy_max_DEV);
    if (Biharm6_const_xx_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Biharm6_const_xx.data, Biharm6_const_xx_DEV);

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::hor_visc_leithy_m(bxKh, *CS_HOST, mask2dBu_DEV.arr,
                           Ah_bg_xx_DEV.arr, m_const_leithy_DEV.arr, m_leithy_max_DEV.arr, Biharm6_const_xx_DEV.arr,
                           inv_PI6,
                           Del2vort_q_DEV.arr, vert_vort_mag_DEV.arr, vort_xy_smooth_DEV.arr,
                           m_leithy_DEV.arr);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(m_leithy_DEV, m_leithy_HOST->data);

    /// Free everything we allocated (freeing a never-allocated container is a no-op)
    turbotmp::free_array4(mask2dBu_DEV);
    turbotmp::free_array4(Del2vort_q_DEV);
    turbotmp::free_array4(vert_vort_mag_DEV);
    turbotmp::free_array4(vort_xy_smooth_DEV);
    turbotmp::free_array4(m_leithy_DEV);
    turbotmp::free_array4(Ah_bg_xx_DEV);
    turbotmp::free_array4(m_const_leithy_DEV);
    turbotmp::free_array4(m_leithy_max_DEV);
    turbotmp::free_array4(Biharm6_const_xx_DEV);
}

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
                                        RealArray_C* Ah_sq_HOST)
{
    (void)nkblock;

    /// Reconstruct the iteration Box (Fortran 1-based -> AMReX 0-based)
    Box bxKh = turbotmp::to_box(bxKh_HOST);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto Del2vort_q_DEV           = turbotmp::make_array4(Del2vort_q_HOST);
    auto vert_vort_mag_smooth_DEV = turbotmp::make_array4(vert_vort_mag_smooth_HOST);
    auto Ah_DEV                   = turbotmp::make_array4(Ah_HOST);
    auto m_leithy_DEV             = turbotmp::make_array4(m_leithy_HOST);
    auto Ah_sq_DEV                = turbotmp::make_array4(Ah_sq_HOST);

    /// The metric arrays may be unallocated (data == nullptr); only allocate/copy them when allocated.
    const bool Ah_bg_xx_allocated         = (metrics_HOST->Ah_bg_xx.data != nullptr);
    const bool Biharm6_const_xx_allocated = (metrics_HOST->Biharm6_const_xx.data != nullptr);
    turbotmp::A4Box Ah_bg_xx_DEV{};
    turbotmp::A4Box Biharm6_const_xx_DEV{};
    if (Ah_bg_xx_allocated)         Ah_bg_xx_DEV         = turbotmp::make_array4(&metrics_HOST->Ah_bg_xx);
    if (Biharm6_const_xx_allocated) Biharm6_const_xx_DEV = turbotmp::make_array4(&metrics_HOST->Biharm6_const_xx);

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(Del2vort_q_HOST->data,            Del2vort_q_DEV);
    turbotmp::copy_FortranHost_to_array4(vert_vort_mag_smooth_HOST->data,  vert_vort_mag_smooth_DEV);
    turbotmp::copy_FortranHost_to_array4(Ah_HOST->data,                    Ah_DEV);
    turbotmp::copy_FortranHost_to_array4(m_leithy_HOST->data,              m_leithy_DEV);
    turbotmp::copy_FortranHost_to_array4(Ah_sq_HOST->data,                 Ah_sq_DEV);
    if (Ah_bg_xx_allocated)         turbotmp::copy_FortranHost_to_array4(metrics_HOST->Ah_bg_xx.data, Ah_bg_xx_DEV);
    if (Biharm6_const_xx_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Biharm6_const_xx.data, Biharm6_const_xx_DEV);

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::hor_visc_leithy_ah(bxKh, *CS_HOST,
                            Ah_bg_xx_DEV.arr, Biharm6_const_xx_DEV.arr,
                            inv_PI6,
                            Del2vort_q_DEV.arr, vert_vort_mag_smooth_DEV.arr,
                            Ah_DEV.arr, m_leithy_DEV.arr, Ah_sq_DEV.arr);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(Ah_DEV,    Ah_HOST->data);
    turbotmp::copy_array4_to_FortranHost(Ah_sq_DEV, Ah_sq_HOST->data);

    /// Free everything we allocated (freeing a never-allocated container is a no-op)
    turbotmp::free_array4(Del2vort_q_DEV);
    turbotmp::free_array4(vert_vort_mag_smooth_DEV);
    turbotmp::free_array4(Ah_DEV);
    turbotmp::free_array4(m_leithy_DEV);
    turbotmp::free_array4(Ah_sq_DEV);
    turbotmp::free_array4(Ah_bg_xx_DEV);
    turbotmp::free_array4(Biharm6_const_xx_DEV);
}

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
                                          const RealArray_C* Ah_sq_HOST)
{
    (void)nkblock;

    /// Reconstruct the iteration Box (Fortran 1-based -> AMReX 0-based)
    Box bxKh = turbotmp::to_box(bxKh_HOST);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto Ah_DEV    = turbotmp::make_array4(Ah_HOST);
    auto Ah_h_DEV  = turbotmp::make_array4(Ah_h_HOST);
    auto Ah_sq_DEV = turbotmp::make_array4(Ah_sq_HOST);

    /// The metric arrays may be unallocated (data == nullptr); only allocate/copy them when allocated.
    const bool Ah_bg_xx_allocated = (metrics_HOST->Ah_bg_xx.data != nullptr);
    turbotmp::A4Box Ah_bg_xx_DEV{};
    if (Ah_bg_xx_allocated) Ah_bg_xx_DEV = turbotmp::make_array4(&metrics_HOST->Ah_bg_xx);

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(Ah_HOST->data,    Ah_DEV);
    turbotmp::copy_FortranHost_to_array4(Ah_h_HOST->data,  Ah_h_DEV);
    turbotmp::copy_FortranHost_to_array4(Ah_sq_HOST->data, Ah_sq_DEV);
    if (Ah_bg_xx_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Ah_bg_xx.data, Ah_bg_xx_DEV);

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::hor_visc_leithy_ah_h(bxKh, *CS_HOST,
                              Ah_bg_xx_DEV.arr,
                              kstart,
                              Ah_DEV.arr, Ah_h_DEV.arr, Ah_sq_DEV.arr);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(Ah_DEV,   Ah_HOST->data);
    turbotmp::copy_array4_to_FortranHost(Ah_h_DEV, Ah_h_HOST->data);

    /// Free everything we allocated (freeing a never-allocated container is a no-op)
    turbotmp::free_array4(Ah_DEV);
    turbotmp::free_array4(Ah_h_DEV);
    turbotmp::free_array4(Ah_sq_DEV);
    turbotmp::free_array4(Ah_bg_xx_DEV);
}
