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
/// @param bxHstr_HOST            The h points of this tile where the normal stresses are computed,
///                               [IscB:IecB+1, JscB:JecB+1, 1:kmax]
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
void turbotmp_hor_visc_backscatter_h_bridge(const Box_C* bxHstr_HOST,
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
    Box bxHstr = turbotmp::to_box(bxHstr_HOST);

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
    MOM::hor_visc_backscatter_h(bxHstr, *CS_HOST,
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
/// @param bxQstr_HOST            The q points of this tile where the shear stresses are computed,
///                               [isc-1:iec, jsc-1:jec, 1:kmax]
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
void turbotmp_hor_visc_backscatter_q_bridge(const Box_C* bxQstr_HOST,
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
    Box bxQstr = turbotmp::to_box(bxQstr_HOST);

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
    MOM::hor_visc_backscatter_q(bxQstr, *CS_HOST,
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
/// @param bxKh_HOST               The thickness-point viscosity box of this tile, [is_Kh:ie_Kh, js_Kh:je_Kh, 1:kmax]
/// @param bxKh_i_HOST             The thickness-point viscosity range in i only, [is_Kh:ie_Kh, jsc:jec, 1:kmax]
/// @param bxKh_j_HOST             The thickness-point viscosity range in j only, [isc:iec, js_Kh:je_Kh, 1:kmax]
/// @param bxQstr_HOST             The q points of this tile where the shear stresses are computed,
///                                [isc-1:iec, jsc-1:jec, 1:kmax]
/// @param bxQvsmooth_HOST         The q points of this tile where the smoothed vorticity terms are computed,
///                                [is_Kh-1:ie_Kh, js_Kh-1:je_Kh, 1:kmax]
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
                                           const Box_C* bxQstr_HOST,
                                           const Box_C* bxQvsmooth_HOST,
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
    Box bxH        = turbotmp::to_box(bxH_HOST);
    Box bxKh       = turbotmp::to_box(bxKh_HOST);
    Box bxKh_i     = turbotmp::to_box(bxKh_i_HOST);
    Box bxKh_j     = turbotmp::to_box(bxKh_j_HOST);
    Box bxQstr     = turbotmp::to_box(bxQstr_HOST);
    Box bxQvsmooth = turbotmp::to_box(bxQvsmooth_HOST);

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
    MOM::hor_visc_leith_grad_a(bxH, bxKh, bxKh_i, bxKh_j, bxQstr, bxQvsmooth, *CS_HOST,
                               dyBu_DEV.arr, IdxBu_DEV.arr, IdyCu_DEV.arr, dxBu_DEV.arr, IdyBu_DEV.arr,
                               IdxCv_DEV.arr, IdyCv_DEV.arr, IdxCu_DEV.arr, dF_dx_DEV.arr, dF_dy_DEV.arr,
                               dudx_DEV.arr, dvdy_DEV.arr, vort_xy_DEV.arr, vort_xy_smooth_DEV.arr,
                               grad_vort_mag_h_2d_DEV.arr, grad_div_mag_h_DEV.arr, grad_vort_mag_q_2d_DEV.arr,
                               grad_div_mag_q_DEV.arr,
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
/// @param bxKh_HOST                 The thickness-point viscosity box of this tile, [is_Kh:ie_Kh, js_Kh:je_Kh, 1:kmax]
/// @param bxQstr_HOST               The q points of this tile where the shear stresses are computed,
///                                  [isc-1:iec, jsc-1:jec, 1:kmax]
/// @param CS_HOST                   Options used by the horizontal viscosity kernels
/// @param nkblock                   The k-block size used to size the following arrays [nondim]
/// @param grad_vort_mag_h_HOST      Magnitude of vorticity gradient at h-points [L-1 T-1 ~> m-1 s-1]
/// @param grad_vort_mag_q_HOST      Magnitude of vorticity gradient at q-points [L-1 T-1 ~> m-1 s-1]
/// @param vert_vort_mag_smooth_HOST Magnitude of gradient of smoothed vertical vorticity (h or q) [L-1 T-1 ~> m-1 s-1]
/// @param vort_xy_dx_HOST           x-derivative of vertical vorticity (d/dx(dv/dx - du/dy)) [L-1 T-1 ~> m-1 s-1]
/// @param vort_xy_dx_smooth_HOST    x-derivative of smoothed vertical vorticity [L-1 T-1 ~> m-1 s-1]
/// @param vort_xy_dy_HOST           y-derivative of vertical vorticity (d/dy(dv/dx - du/dy)) [L-1 T-1 ~> m-1 s-1]
/// @param vort_xy_dy_smooth_HOST    y-derivative of smoothed vertical vorticity [L-1 T-1 ~> m-1 s-1]
void turbotmp_hor_visc_leith_grad_b_bridge(const Box_C* bxKh_HOST,
                                           const Box_C* bxQstr_HOST,
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
    Box bxKh   = turbotmp::to_box(bxKh_HOST);
    Box bxQstr = turbotmp::to_box(bxQstr_HOST);

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
    MOM::hor_visc_leith_grad_b(bxKh, bxQstr, *CS_HOST,
                               grad_vort_mag_h_DEV.arr, grad_vort_mag_q_DEV.arr, vert_vort_mag_smooth_DEV.arr,
                               vort_xy_dx_DEV.arr, vort_xy_dx_smooth_DEV.arr, vort_xy_dy_DEV.arr,
                               vort_xy_dy_smooth_DEV.arr);

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
    if (m_const_leithy_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->m_const_leithy.data, m_const_leithy_DEV);
    if (m_leithy_max_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->m_leithy_max.data, m_leithy_max_DEV);
    if (Biharm6_const_xx_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->Biharm6_const_xx.data, Biharm6_const_xx_DEV);

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
    if (Biharm6_const_xx_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->Biharm6_const_xx.data, Biharm6_const_xx_DEV);

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
/// @param bxHstr_HOST                The h points of this tile where the normal stresses are computed,
///                                   [IscB:IecB+1, JscB:JecB+1, 1:kmax]
/// @param bxQstr_HOST                The q points of this tile where the shear stresses are computed,
///                                   [isc-1:iec, jsc-1:jec, 1:kmax]
/// @param bxHgrad_HOST               The h points of this tile where the velocity gradients are computed,
///                                   [IscB-1:IecB+2, JscB-1:JecB+2, 1:kmax]
/// @param bxQshear_HOST              The q points of this tile where the shearing strain is computed,
///                                   [isc-2:iec+1, jsc-2:jec+1, 1:kmax]
/// @param bxQvsmooth_HOST            The q points of this tile where the smoothed vorticity terms are computed,
///                                   [is_Kh-1:ie_Kh, js_Kh-1:je_Kh, 1:kmax]
/// @param nkblock                    The effective k-block size [nondim]
/// @param kstart                     The first absolute k-layer of the current k-block
/// @param use_vort_xy                True if vort_xy must be computed
/// @param apply_OBC                  If true, open boundary conditions are applied on this PE
/// @param apply_OBC_strain           If true, the specified open boundary strain configuration is applied
/// @param use_cont_huv               If true, use the continuity solver thicknesses at velocity points
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
                                     const Box_C* bxHstr_HOST,
                                     const Box_C* bxQstr_HOST,
                                     const Box_C* bxHgrad_HOST,
                                     const Box_C* bxQshear_HOST,
                                     const Box_C* bxQvsmooth_HOST,
                                     const int nkblock,
                                     const int kstart,
                                     const bool use_vort_xy,
                                     const bool apply_OBC,
                                     const bool apply_OBC_strain,
                                     const bool use_cont_huv,
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
                                     RealArray_C* vort_xy_smooth_HOST)
{
    (void)nkblock;

    /// Reconstruct the iteration Boxes (Fortran 1-based -> AMReX 0-based)
    Box bxH = turbotmp::to_box(bxH_HOST);
    Box bxQ = turbotmp::to_box(bxQ_HOST);
    Box bxU = turbotmp::to_box(bxU_HOST);
    Box bxV = turbotmp::to_box(bxV_HOST);
    Box bxKh = turbotmp::to_box(bxKh_HOST);
    Box bxHstr = turbotmp::to_box(bxHstr_HOST);
    Box bxQstr = turbotmp::to_box(bxQstr_HOST);
    Box bxHgrad = turbotmp::to_box(bxHgrad_HOST);
    Box bxQshear = turbotmp::to_box(bxQshear_HOST);
    Box bxQvsmooth = turbotmp::to_box(bxQvsmooth_HOST);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto IdyCu_DEV          = turbotmp::make_array4(&Gcore_HOST->IdyCu);
    auto IdxCv_DEV          = turbotmp::make_array4(&Gcore_HOST->IdxCv);
    auto IdyCv_DEV          = turbotmp::make_array4(&Gcore_HOST->IdyCv);
    auto IdxCu_DEV          = turbotmp::make_array4(&Gcore_HOST->IdxCu);
    auto mask2dT_DEV        = turbotmp::make_array4(&Gcore_HOST->mask2dT);
    auto dxBu_DEV           = turbotmp::make_array4(&Gcore_HOST->dxBu);
    auto mask2dBu_DEV       = turbotmp::make_array4(&Gcore_HOST->mask2dBu);
    auto IareaBu_DEV        = turbotmp::make_array4(&Gcore_HOST->IareaBu);
    auto dyCv_DEV           = turbotmp::make_array4(&Gcore_HOST->dyCv);
    auto dxCu_DEV           = turbotmp::make_array4(&Gcore_HOST->dxCu);
    auto u_DEV              = turbotmp::make_array4(u_HOST);
    auto v_DEV              = turbotmp::make_array4(v_HOST);
    auto h_DEV              = turbotmp::make_array4(h_HOST);
    auto u_smooth_DEV       = turbotmp::make_array4(u_smooth_HOST);
    auto v_smooth_DEV       = turbotmp::make_array4(v_smooth_HOST);
    auto NoSt_DEV           = turbotmp::make_array4(NoSt_HOST);
    auto ShSt_DEV           = turbotmp::make_array4(ShSt_HOST);
    auto Del2u_DEV          = turbotmp::make_array4(Del2u_HOST);
    auto h_u_DEV            = turbotmp::make_array4(h_u_HOST);
    auto Del2v_DEV          = turbotmp::make_array4(Del2v_HOST);
    auto h_v_DEV            = turbotmp::make_array4(h_v_HOST);
    auto sh_xx_DEV          = turbotmp::make_array4(sh_xx_HOST);
    auto sh_xx_smooth_DEV   = turbotmp::make_array4(sh_xx_smooth_HOST);
    auto dudx_DEV           = turbotmp::make_array4(dudx_HOST);
    auto dvdy_DEV           = turbotmp::make_array4(dvdy_HOST);
    auto dvdx_DEV           = turbotmp::make_array4(dvdx_HOST);
    auto dudy_DEV           = turbotmp::make_array4(dudy_HOST);
    auto dvdx_smooth_DEV    = turbotmp::make_array4(dvdx_smooth_HOST);
    auto dudy_smooth_DEV    = turbotmp::make_array4(dudy_smooth_HOST);
    auto sh_xy_DEV          = turbotmp::make_array4(sh_xy_HOST);
    auto sh_xy_smooth_DEV   = turbotmp::make_array4(sh_xy_smooth_HOST);
    auto vort_xy_DEV        = turbotmp::make_array4(vort_xy_HOST);
    auto vort_xy_smooth_DEV = turbotmp::make_array4(vort_xy_smooth_HOST);

    /// The metric arrays may be unallocated (data == nullptr); only allocate/copy them when allocated.
    const bool dx2h_allocated     = (metrics_HOST->dx2h.data != nullptr);
    const bool dy2h_allocated     = (metrics_HOST->dy2h.data != nullptr);
    const bool dx_dyT_allocated   = (metrics_HOST->dx_dyT.data != nullptr);
    const bool dy_dxT_allocated   = (metrics_HOST->dy_dxT.data != nullptr);
    const bool dx2q_allocated     = (metrics_HOST->dx2q.data != nullptr);
    const bool dy2q_allocated     = (metrics_HOST->dy2q.data != nullptr);
    const bool dx_dyBu_allocated  = (metrics_HOST->dx_dyBu.data != nullptr);
    const bool dy_dxBu_allocated  = (metrics_HOST->dy_dxBu.data != nullptr);
    const bool Idx2dyCu_allocated = (metrics_HOST->Idx2dyCu.data != nullptr);
    const bool Idxdy2u_allocated  = (metrics_HOST->Idxdy2u.data != nullptr);
    const bool Idx2dyCv_allocated = (metrics_HOST->Idx2dyCv.data != nullptr);
    const bool Idxdy2v_allocated  = (metrics_HOST->Idxdy2v.data != nullptr);
    turbotmp::A4Box dx2h_DEV{};
    turbotmp::A4Box dy2h_DEV{};
    turbotmp::A4Box dx_dyT_DEV{};
    turbotmp::A4Box dy_dxT_DEV{};
    turbotmp::A4Box dx2q_DEV{};
    turbotmp::A4Box dy2q_DEV{};
    turbotmp::A4Box dx_dyBu_DEV{};
    turbotmp::A4Box dy_dxBu_DEV{};
    turbotmp::A4Box Idx2dyCu_DEV{};
    turbotmp::A4Box Idxdy2u_DEV{};
    turbotmp::A4Box Idx2dyCv_DEV{};
    turbotmp::A4Box Idxdy2v_DEV{};
    if (dx2h_allocated) dx2h_DEV = turbotmp::make_array4(&metrics_HOST->dx2h);
    if (dy2h_allocated) dy2h_DEV = turbotmp::make_array4(&metrics_HOST->dy2h);
    if (dx_dyT_allocated) dx_dyT_DEV = turbotmp::make_array4(&metrics_HOST->dx_dyT);
    if (dy_dxT_allocated) dy_dxT_DEV = turbotmp::make_array4(&metrics_HOST->dy_dxT);
    if (dx2q_allocated) dx2q_DEV = turbotmp::make_array4(&metrics_HOST->dx2q);
    if (dy2q_allocated) dy2q_DEV = turbotmp::make_array4(&metrics_HOST->dy2q);
    if (dx_dyBu_allocated) dx_dyBu_DEV = turbotmp::make_array4(&metrics_HOST->dx_dyBu);
    if (dy_dxBu_allocated) dy_dxBu_DEV = turbotmp::make_array4(&metrics_HOST->dy_dxBu);
    if (Idx2dyCu_allocated) Idx2dyCu_DEV = turbotmp::make_array4(&metrics_HOST->Idx2dyCu);
    if (Idxdy2u_allocated) Idxdy2u_DEV = turbotmp::make_array4(&metrics_HOST->Idxdy2u);
    if (Idx2dyCv_allocated) Idx2dyCv_DEV = turbotmp::make_array4(&metrics_HOST->Idx2dyCv);
    if (Idxdy2v_allocated) Idxdy2v_DEV = turbotmp::make_array4(&metrics_HOST->Idxdy2v);

    /// The continuity solver thicknesses may be unallocated (data == nullptr); only allocate/copy them when allocated.
    const bool hu_cont_allocated = (hu_cont_HOST->data != nullptr);
    const bool hv_cont_allocated = (hv_cont_HOST->data != nullptr);
    turbotmp::A4Box hu_cont_DEV{};
    turbotmp::A4Box hv_cont_DEV{};
    if (hu_cont_allocated) hu_cont_DEV = turbotmp::make_array4(hu_cont_HOST);
    if (hv_cont_allocated) hv_cont_DEV = turbotmp::make_array4(hv_cont_HOST);

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdyCu.data,    IdyCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdxCv.data,    IdxCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdyCv.data,    IdyCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdxCu.data,    IdxCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dT.data,  mask2dT_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->dxBu.data,     dxBu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dBu.data, mask2dBu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IareaBu.data,  IareaBu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->dyCv.data,     dyCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->dxCu.data,     dxCu_DEV);
    turbotmp::copy_FortranHost_to_array4(u_HOST->data,              u_DEV);
    turbotmp::copy_FortranHost_to_array4(v_HOST->data,              v_DEV);
    turbotmp::copy_FortranHost_to_array4(h_HOST->data,              h_DEV);
    turbotmp::copy_FortranHost_to_array4(u_smooth_HOST->data,       u_smooth_DEV);
    turbotmp::copy_FortranHost_to_array4(v_smooth_HOST->data,       v_smooth_DEV);
    turbotmp::copy_FortranHost_to_array4(NoSt_HOST->data,           NoSt_DEV);
    turbotmp::copy_FortranHost_to_array4(ShSt_HOST->data,           ShSt_DEV);
    turbotmp::copy_FortranHost_to_array4(Del2u_HOST->data,          Del2u_DEV);
    turbotmp::copy_FortranHost_to_array4(h_u_HOST->data,            h_u_DEV);
    turbotmp::copy_FortranHost_to_array4(Del2v_HOST->data,          Del2v_DEV);
    turbotmp::copy_FortranHost_to_array4(h_v_HOST->data,            h_v_DEV);
    turbotmp::copy_FortranHost_to_array4(sh_xx_HOST->data,          sh_xx_DEV);
    turbotmp::copy_FortranHost_to_array4(sh_xx_smooth_HOST->data,   sh_xx_smooth_DEV);
    turbotmp::copy_FortranHost_to_array4(dudx_HOST->data,           dudx_DEV);
    turbotmp::copy_FortranHost_to_array4(dvdy_HOST->data,           dvdy_DEV);
    turbotmp::copy_FortranHost_to_array4(dvdx_HOST->data,           dvdx_DEV);
    turbotmp::copy_FortranHost_to_array4(dudy_HOST->data,           dudy_DEV);
    turbotmp::copy_FortranHost_to_array4(dvdx_smooth_HOST->data,    dvdx_smooth_DEV);
    turbotmp::copy_FortranHost_to_array4(dudy_smooth_HOST->data,    dudy_smooth_DEV);
    turbotmp::copy_FortranHost_to_array4(sh_xy_HOST->data,          sh_xy_DEV);
    turbotmp::copy_FortranHost_to_array4(sh_xy_smooth_HOST->data,   sh_xy_smooth_DEV);
    turbotmp::copy_FortranHost_to_array4(vort_xy_HOST->data,        vort_xy_DEV);
    turbotmp::copy_FortranHost_to_array4(vort_xy_smooth_HOST->data, vort_xy_smooth_DEV);
    if (dx2h_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->dx2h.data, dx2h_DEV);
    if (dy2h_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->dy2h.data, dy2h_DEV);
    if (dx_dyT_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->dx_dyT.data, dx_dyT_DEV);
    if (dy_dxT_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->dy_dxT.data, dy_dxT_DEV);
    if (dx2q_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->dx2q.data, dx2q_DEV);
    if (dy2q_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->dy2q.data, dy2q_DEV);
    if (dx_dyBu_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->dx_dyBu.data, dx_dyBu_DEV);
    if (dy_dxBu_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->dy_dxBu.data, dy_dxBu_DEV);
    if (Idx2dyCu_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Idx2dyCu.data, Idx2dyCu_DEV);
    if (Idxdy2u_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Idxdy2u.data, Idxdy2u_DEV);
    if (Idx2dyCv_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Idx2dyCv.data, Idx2dyCv_DEV);
    if (Idxdy2v_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Idxdy2v.data, Idxdy2v_DEV);
    if (hu_cont_allocated) turbotmp::copy_FortranHost_to_array4(hu_cont_HOST->data, hu_cont_DEV);
    if (hv_cont_allocated) turbotmp::copy_FortranHost_to_array4(hv_cont_HOST->data, hv_cont_DEV);

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::hor_visc_strain(bxH, bxQ, bxU, bxV, bxKh, bxHstr, bxQstr, bxHgrad, bxQshear, bxQvsmooth,
                         *CS_HOST,
                         IdyCu_DEV.arr, IdxCv_DEV.arr, IdyCv_DEV.arr, IdxCu_DEV.arr, mask2dT_DEV.arr,
                         dxBu_DEV.arr, mask2dBu_DEV.arr, IareaBu_DEV.arr, dyCv_DEV.arr, dxCu_DEV.arr,
                         dx2h_DEV.arr, dy2h_DEV.arr, dx_dyT_DEV.arr, dy_dxT_DEV.arr,
                         dx2q_DEV.arr, dy2q_DEV.arr, dx_dyBu_DEV.arr, dy_dxBu_DEV.arr,
                         Idx2dyCu_DEV.arr, Idxdy2u_DEV.arr, Idx2dyCv_DEV.arr, Idxdy2v_DEV.arr,
                         OBC,
                         kstart,
                         use_vort_xy, apply_OBC, apply_OBC_strain, use_cont_huv,
                         u_DEV.arr, v_DEV.arr, h_DEV.arr,
                         hu_cont_DEV.arr, hv_cont_DEV.arr,
                         u_smooth_DEV.arr, v_smooth_DEV.arr,
                         NoSt_DEV.arr, ShSt_DEV.arr, Del2u_DEV.arr, h_u_DEV.arr,
                         Del2v_DEV.arr, h_v_DEV.arr, sh_xx_DEV.arr, sh_xx_smooth_DEV.arr,
                         dudx_DEV.arr, dvdy_DEV.arr, dvdx_DEV.arr, dudy_DEV.arr,
                         dvdx_smooth_DEV.arr, dudy_smooth_DEV.arr, sh_xy_DEV.arr, sh_xy_smooth_DEV.arr,
                         vort_xy_DEV.arr, vort_xy_smooth_DEV.arr);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(NoSt_DEV,           NoSt_HOST->data);
    turbotmp::copy_array4_to_FortranHost(ShSt_DEV,           ShSt_HOST->data);
    turbotmp::copy_array4_to_FortranHost(Del2u_DEV,          Del2u_HOST->data);
    turbotmp::copy_array4_to_FortranHost(h_u_DEV,            h_u_HOST->data);
    turbotmp::copy_array4_to_FortranHost(Del2v_DEV,          Del2v_HOST->data);
    turbotmp::copy_array4_to_FortranHost(h_v_DEV,            h_v_HOST->data);
    turbotmp::copy_array4_to_FortranHost(sh_xx_DEV,          sh_xx_HOST->data);
    turbotmp::copy_array4_to_FortranHost(sh_xx_smooth_DEV,   sh_xx_smooth_HOST->data);
    turbotmp::copy_array4_to_FortranHost(dudx_DEV,           dudx_HOST->data);
    turbotmp::copy_array4_to_FortranHost(dvdy_DEV,           dvdy_HOST->data);
    turbotmp::copy_array4_to_FortranHost(dvdx_DEV,           dvdx_HOST->data);
    turbotmp::copy_array4_to_FortranHost(dudy_DEV,           dudy_HOST->data);
    turbotmp::copy_array4_to_FortranHost(dvdx_smooth_DEV,    dvdx_smooth_HOST->data);
    turbotmp::copy_array4_to_FortranHost(dudy_smooth_DEV,    dudy_smooth_HOST->data);
    turbotmp::copy_array4_to_FortranHost(sh_xy_DEV,          sh_xy_HOST->data);
    turbotmp::copy_array4_to_FortranHost(sh_xy_smooth_DEV,   sh_xy_smooth_HOST->data);
    turbotmp::copy_array4_to_FortranHost(vort_xy_DEV,        vort_xy_HOST->data);
    turbotmp::copy_array4_to_FortranHost(vort_xy_smooth_DEV, vort_xy_smooth_HOST->data);
    if (hu_cont_allocated) turbotmp::copy_array4_to_FortranHost(hu_cont_DEV, hu_cont_HOST->data);
    if (hv_cont_allocated) turbotmp::copy_array4_to_FortranHost(hv_cont_DEV, hv_cont_HOST->data);

    /// Free everything we allocated (freeing a never-allocated container is a no-op)
    turbotmp::free_array4(IdyCu_DEV);
    turbotmp::free_array4(IdxCv_DEV);
    turbotmp::free_array4(IdyCv_DEV);
    turbotmp::free_array4(IdxCu_DEV);
    turbotmp::free_array4(mask2dT_DEV);
    turbotmp::free_array4(dxBu_DEV);
    turbotmp::free_array4(mask2dBu_DEV);
    turbotmp::free_array4(IareaBu_DEV);
    turbotmp::free_array4(dyCv_DEV);
    turbotmp::free_array4(dxCu_DEV);
    turbotmp::free_array4(u_DEV);
    turbotmp::free_array4(v_DEV);
    turbotmp::free_array4(h_DEV);
    turbotmp::free_array4(hu_cont_DEV);
    turbotmp::free_array4(hv_cont_DEV);
    turbotmp::free_array4(u_smooth_DEV);
    turbotmp::free_array4(v_smooth_DEV);
    turbotmp::free_array4(NoSt_DEV);
    turbotmp::free_array4(ShSt_DEV);
    turbotmp::free_array4(Del2u_DEV);
    turbotmp::free_array4(h_u_DEV);
    turbotmp::free_array4(Del2v_DEV);
    turbotmp::free_array4(h_v_DEV);
    turbotmp::free_array4(sh_xx_DEV);
    turbotmp::free_array4(sh_xx_smooth_DEV);
    turbotmp::free_array4(dudx_DEV);
    turbotmp::free_array4(dvdy_DEV);
    turbotmp::free_array4(dvdx_DEV);
    turbotmp::free_array4(dudy_DEV);
    turbotmp::free_array4(dvdx_smooth_DEV);
    turbotmp::free_array4(dudy_smooth_DEV);
    turbotmp::free_array4(sh_xy_DEV);
    turbotmp::free_array4(sh_xy_smooth_DEV);
    turbotmp::free_array4(vort_xy_DEV);
    turbotmp::free_array4(vort_xy_smooth_DEV);
    turbotmp::free_array4(dx2h_DEV);
    turbotmp::free_array4(dy2h_DEV);
    turbotmp::free_array4(dx_dyT_DEV);
    turbotmp::free_array4(dy_dxT_DEV);
    turbotmp::free_array4(dx2q_DEV);
    turbotmp::free_array4(dy2q_DEV);
    turbotmp::free_array4(dx_dyBu_DEV);
    turbotmp::free_array4(dy_dxBu_DEV);
    turbotmp::free_array4(Idx2dyCu_DEV);
    turbotmp::free_array4(Idxdy2u_DEV);
    turbotmp::free_array4(Idx2dyCv_DEV);
    turbotmp::free_array4(Idxdy2v_DEV);
}

/// @brief Bridge for the MOM6 hor_visc_h_coeffs_a subroutine.
///
/// @param MEKE_HOST               MEKE fields related to Mesoscale Eddy Kinetic Energy; each field may be unallocated
///                                (data == nullptr)
/// @param VarMix_HOST             Variable mixing control structure; each field may be unallocated (data == nullptr)
/// @param CS_HOST                 Options used by the horizontal viscosity kernels
/// @param metrics_HOST            Persistent copies of the horizontal viscosity metric and coefficient arrays; each
///                                field may be unallocated (data == nullptr)
/// @param bxH_HOST                The h-point iteration box of this tile, [isc:iec, jsc:jec, 1:kmax]
/// @param bxKh_HOST               The thickness-point viscosity box of this tile, [is_Kh:ie_Kh, js_Kh:je_Kh, 1:kmax]
/// @param bxHstr_HOST             The h points of this tile where the normal stresses are computed,
///                                [IscB:IecB+1, JscB:JecB+1, 1:kmax]
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
                                         const Box_C* bxKh_HOST,
                                         const Box_C* bxHstr_HOST,
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
                                         RealArray_C* visc_bound_rem_HOST)
{
    (void)nkblock;

    /// Reconstruct the iteration Boxes (Fortran 1-based -> AMReX 0-based)
    Box bxH = turbotmp::to_box(bxH_HOST);
    Box bxKh = turbotmp::to_box(bxKh_HOST);
    Box bxHstr = turbotmp::to_box(bxHstr_HOST);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto u_DEV                  = turbotmp::make_array4(u_HOST);
    auto v_DEV                  = turbotmp::make_array4(v_HOST);
    auto h_DEV                  = turbotmp::make_array4(h_HOST);
    auto Kh_h_DEV               = turbotmp::make_array4(Kh_h_HOST);
    auto grid_Re_Kh_DEV         = turbotmp::make_array4(grid_Re_Kh_HOST);
    auto div_xx_h_DEV           = turbotmp::make_array4(div_xx_h_HOST);
    auto sh_xx_h_DEV            = turbotmp::make_array4(sh_xx_h_HOST);
    auto h_u_DEV                = turbotmp::make_array4(h_u_HOST);
    auto h_v_DEV                = turbotmp::make_array4(h_v_HOST);
    auto sh_xx_DEV              = turbotmp::make_array4(sh_xx_HOST);
    auto str_xx_DEV             = turbotmp::make_array4(str_xx_HOST);
    auto grad_vort_mag_h_DEV    = turbotmp::make_array4(grad_vort_mag_h_HOST);
    auto grad_vort_mag_h_2d_DEV = turbotmp::make_array4(grad_vort_mag_h_2d_HOST);
    auto grad_div_mag_h_DEV     = turbotmp::make_array4(grad_div_mag_h_HOST);
    auto dudx_DEV               = turbotmp::make_array4(dudx_HOST);
    auto dvdy_DEV               = turbotmp::make_array4(dvdy_HOST);
    auto sh_xy_DEV              = turbotmp::make_array4(sh_xy_HOST);
    auto Del2vort_q_DEV         = turbotmp::make_array4(Del2vort_q_HOST);
    auto Ah_DEV                 = turbotmp::make_array4(Ah_HOST);
    auto Kh_DEV                 = turbotmp::make_array4(Kh_HOST);
    auto Shear_mag_DEV          = turbotmp::make_array4(Shear_mag_HOST);
    auto vert_vort_mag_DEV      = turbotmp::make_array4(vert_vort_mag_HOST);
    auto hrat_min_DEV           = turbotmp::make_array4(hrat_min_HOST);
    auto visc_bound_rem_DEV     = turbotmp::make_array4(visc_bound_rem_HOST);

    /// The MEKE, VarMix and metric arrays may be unallocated (data == nullptr); only allocate/copy them when
    /// allocated.
    const bool MEKE_Ku_allocated          = (MEKE_HOST->Ku.data != nullptr);
    const bool Res_fn_h_allocated         = (VarMix_HOST->Res_fn_h.data != nullptr);
    const bool BS_struct_allocated        = (VarMix_HOST->BS_struct.data != nullptr);
    const bool Kh_bg_xx_allocated         = (metrics_HOST->Kh_bg_xx.data != nullptr);
    const bool Ah_bg_xx_allocated         = (metrics_HOST->Ah_bg_xx.data != nullptr);
    const bool Kh_Max_xx_allocated        = (metrics_HOST->Kh_Max_xx.data != nullptr);
    const bool n1n2_h_allocated           = (metrics_HOST->n1n2_h.data != nullptr);
    const bool n1n1_m_n2n2_h_allocated    = (metrics_HOST->n1n1_m_n2n2_h.data != nullptr);
    const bool grid_sp_h2_allocated       = (metrics_HOST->grid_sp_h2.data != nullptr);
    const bool Laplac2_const_xx_allocated = (metrics_HOST->Laplac2_const_xx.data != nullptr);
    const bool Biharm6_const_xx_allocated = (metrics_HOST->Biharm6_const_xx.data != nullptr);
    const bool Laplac3_const_xx_allocated = (metrics_HOST->Laplac3_const_xx.data != nullptr);
    const bool Biharm_const_xx_allocated  = (metrics_HOST->Biharm_const_xx.data != nullptr);
    const bool Biharm_const2_xx_allocated = (metrics_HOST->Biharm_const2_xx.data != nullptr);
    turbotmp::A4Box MEKE_Ku_DEV{};
    turbotmp::A4Box Res_fn_h_DEV{};
    turbotmp::A4Box BS_struct_DEV{};
    turbotmp::A4Box Kh_bg_xx_DEV{};
    turbotmp::A4Box Ah_bg_xx_DEV{};
    turbotmp::A4Box Kh_Max_xx_DEV{};
    turbotmp::A4Box n1n2_h_DEV{};
    turbotmp::A4Box n1n1_m_n2n2_h_DEV{};
    turbotmp::A4Box grid_sp_h2_DEV{};
    turbotmp::A4Box Laplac2_const_xx_DEV{};
    turbotmp::A4Box Biharm6_const_xx_DEV{};
    turbotmp::A4Box Laplac3_const_xx_DEV{};
    turbotmp::A4Box Biharm_const_xx_DEV{};
    turbotmp::A4Box Biharm_const2_xx_DEV{};
    if (MEKE_Ku_allocated) MEKE_Ku_DEV = turbotmp::make_array4(&MEKE_HOST->Ku);
    if (Res_fn_h_allocated) Res_fn_h_DEV = turbotmp::make_array4(&VarMix_HOST->Res_fn_h);
    if (BS_struct_allocated) BS_struct_DEV = turbotmp::make_array4(&VarMix_HOST->BS_struct);
    if (Kh_bg_xx_allocated) Kh_bg_xx_DEV = turbotmp::make_array4(&metrics_HOST->Kh_bg_xx);
    if (Ah_bg_xx_allocated) Ah_bg_xx_DEV = turbotmp::make_array4(&metrics_HOST->Ah_bg_xx);
    if (Kh_Max_xx_allocated) Kh_Max_xx_DEV = turbotmp::make_array4(&metrics_HOST->Kh_Max_xx);
    if (n1n2_h_allocated) n1n2_h_DEV = turbotmp::make_array4(&metrics_HOST->n1n2_h);
    if (n1n1_m_n2n2_h_allocated) n1n1_m_n2n2_h_DEV = turbotmp::make_array4(&metrics_HOST->n1n1_m_n2n2_h);
    if (grid_sp_h2_allocated) grid_sp_h2_DEV = turbotmp::make_array4(&metrics_HOST->grid_sp_h2);
    if (Laplac2_const_xx_allocated) Laplac2_const_xx_DEV = turbotmp::make_array4(&metrics_HOST->Laplac2_const_xx);
    if (Biharm6_const_xx_allocated) Biharm6_const_xx_DEV = turbotmp::make_array4(&metrics_HOST->Biharm6_const_xx);
    if (Laplac3_const_xx_allocated) Laplac3_const_xx_DEV = turbotmp::make_array4(&metrics_HOST->Laplac3_const_xx);
    if (Biharm_const_xx_allocated) Biharm_const_xx_DEV = turbotmp::make_array4(&metrics_HOST->Biharm_const_xx);
    if (Biharm_const2_xx_allocated) Biharm_const2_xx_DEV = turbotmp::make_array4(&metrics_HOST->Biharm_const2_xx);

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(u_HOST->data,                  u_DEV);
    turbotmp::copy_FortranHost_to_array4(v_HOST->data,                  v_DEV);
    turbotmp::copy_FortranHost_to_array4(h_HOST->data,                  h_DEV);
    turbotmp::copy_FortranHost_to_array4(Kh_h_HOST->data,               Kh_h_DEV);
    turbotmp::copy_FortranHost_to_array4(grid_Re_Kh_HOST->data,         grid_Re_Kh_DEV);
    turbotmp::copy_FortranHost_to_array4(div_xx_h_HOST->data,           div_xx_h_DEV);
    turbotmp::copy_FortranHost_to_array4(sh_xx_h_HOST->data,            sh_xx_h_DEV);
    turbotmp::copy_FortranHost_to_array4(h_u_HOST->data,                h_u_DEV);
    turbotmp::copy_FortranHost_to_array4(h_v_HOST->data,                h_v_DEV);
    turbotmp::copy_FortranHost_to_array4(sh_xx_HOST->data,              sh_xx_DEV);
    turbotmp::copy_FortranHost_to_array4(str_xx_HOST->data,             str_xx_DEV);
    turbotmp::copy_FortranHost_to_array4(grad_vort_mag_h_HOST->data,    grad_vort_mag_h_DEV);
    turbotmp::copy_FortranHost_to_array4(grad_vort_mag_h_2d_HOST->data, grad_vort_mag_h_2d_DEV);
    turbotmp::copy_FortranHost_to_array4(grad_div_mag_h_HOST->data,     grad_div_mag_h_DEV);
    turbotmp::copy_FortranHost_to_array4(dudx_HOST->data,               dudx_DEV);
    turbotmp::copy_FortranHost_to_array4(dvdy_HOST->data,               dvdy_DEV);
    turbotmp::copy_FortranHost_to_array4(sh_xy_HOST->data,              sh_xy_DEV);
    turbotmp::copy_FortranHost_to_array4(Del2vort_q_HOST->data,         Del2vort_q_DEV);
    turbotmp::copy_FortranHost_to_array4(Ah_HOST->data,                 Ah_DEV);
    turbotmp::copy_FortranHost_to_array4(Kh_HOST->data,                 Kh_DEV);
    turbotmp::copy_FortranHost_to_array4(Shear_mag_HOST->data,          Shear_mag_DEV);
    turbotmp::copy_FortranHost_to_array4(vert_vort_mag_HOST->data,      vert_vort_mag_DEV);
    turbotmp::copy_FortranHost_to_array4(hrat_min_HOST->data,           hrat_min_DEV);
    turbotmp::copy_FortranHost_to_array4(visc_bound_rem_HOST->data,     visc_bound_rem_DEV);
    if (MEKE_Ku_allocated) turbotmp::copy_FortranHost_to_array4(MEKE_HOST->Ku.data, MEKE_Ku_DEV);
    if (Res_fn_h_allocated) turbotmp::copy_FortranHost_to_array4(VarMix_HOST->Res_fn_h.data, Res_fn_h_DEV);
    if (BS_struct_allocated) turbotmp::copy_FortranHost_to_array4(VarMix_HOST->BS_struct.data, BS_struct_DEV);
    if (Kh_bg_xx_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Kh_bg_xx.data, Kh_bg_xx_DEV);
    if (Ah_bg_xx_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Ah_bg_xx.data, Ah_bg_xx_DEV);
    if (Kh_Max_xx_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Kh_Max_xx.data, Kh_Max_xx_DEV);
    if (n1n2_h_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->n1n2_h.data, n1n2_h_DEV);
    if (n1n1_m_n2n2_h_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->n1n1_m_n2n2_h.data, n1n1_m_n2n2_h_DEV);
    if (grid_sp_h2_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->grid_sp_h2.data, grid_sp_h2_DEV);
    if (Laplac2_const_xx_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->Laplac2_const_xx.data, Laplac2_const_xx_DEV);
    if (Biharm6_const_xx_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->Biharm6_const_xx.data, Biharm6_const_xx_DEV);
    if (Laplac3_const_xx_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->Laplac3_const_xx.data, Laplac3_const_xx_DEV);
    if (Biharm_const_xx_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->Biharm_const_xx.data, Biharm_const_xx_DEV);
    if (Biharm_const2_xx_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->Biharm_const2_xx.data, Biharm_const2_xx_DEV);

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::hor_visc_h_coeffs_a(bxH, bxKh, bxHstr, *CS_HOST,
                             MEKE_Ku_DEV.arr, Res_fn_h_DEV.arr, BS_struct_DEV.arr,
                             Kh_bg_xx_DEV.arr, Ah_bg_xx_DEV.arr, Kh_Max_xx_DEV.arr, n1n2_h_DEV.arr,
                             n1n1_m_n2n2_h_DEV.arr, grid_sp_h2_DEV.arr, Laplac2_const_xx_DEV.arr,
                             Biharm6_const_xx_DEV.arr, Laplac3_const_xx_DEV.arr, Biharm_const_xx_DEV.arr,
                             Biharm_const2_xx_DEV.arr,
                             kstart, h_neglect, inv_PI3, inv_PI6,
                             use_kh_struct, use_MEKE_Ku, rescale_Kh,
                             u_DEV.arr, v_DEV.arr, h_DEV.arr,
                             Kh_h_DEV.arr, grid_Re_Kh_DEV.arr, div_xx_h_DEV.arr, sh_xx_h_DEV.arr,
                             h_u_DEV.arr, h_v_DEV.arr, sh_xx_DEV.arr, str_xx_DEV.arr,
                             grad_vort_mag_h_DEV.arr, grad_vort_mag_h_2d_DEV.arr, grad_div_mag_h_DEV.arr, dudx_DEV.arr,
                             dvdy_DEV.arr, sh_xy_DEV.arr, Del2vort_q_DEV.arr, Ah_DEV.arr,
                             Kh_DEV.arr, Shear_mag_DEV.arr, vert_vort_mag_DEV.arr, hrat_min_DEV.arr,
                             visc_bound_rem_DEV.arr);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(Kh_h_DEV,               Kh_h_HOST->data);
    turbotmp::copy_array4_to_FortranHost(grid_Re_Kh_DEV,         grid_Re_Kh_HOST->data);
    turbotmp::copy_array4_to_FortranHost(div_xx_h_DEV,           div_xx_h_HOST->data);
    turbotmp::copy_array4_to_FortranHost(sh_xx_h_DEV,            sh_xx_h_HOST->data);
    turbotmp::copy_array4_to_FortranHost(h_u_DEV,                h_u_HOST->data);
    turbotmp::copy_array4_to_FortranHost(h_v_DEV,                h_v_HOST->data);
    turbotmp::copy_array4_to_FortranHost(sh_xx_DEV,              sh_xx_HOST->data);
    turbotmp::copy_array4_to_FortranHost(str_xx_DEV,             str_xx_HOST->data);
    turbotmp::copy_array4_to_FortranHost(grad_vort_mag_h_DEV,    grad_vort_mag_h_HOST->data);
    turbotmp::copy_array4_to_FortranHost(grad_vort_mag_h_2d_DEV, grad_vort_mag_h_2d_HOST->data);
    turbotmp::copy_array4_to_FortranHost(grad_div_mag_h_DEV,     grad_div_mag_h_HOST->data);
    turbotmp::copy_array4_to_FortranHost(dudx_DEV,               dudx_HOST->data);
    turbotmp::copy_array4_to_FortranHost(dvdy_DEV,               dvdy_HOST->data);
    turbotmp::copy_array4_to_FortranHost(sh_xy_DEV,              sh_xy_HOST->data);
    turbotmp::copy_array4_to_FortranHost(Del2vort_q_DEV,         Del2vort_q_HOST->data);
    turbotmp::copy_array4_to_FortranHost(Ah_DEV,                 Ah_HOST->data);
    turbotmp::copy_array4_to_FortranHost(Kh_DEV,                 Kh_HOST->data);
    turbotmp::copy_array4_to_FortranHost(Shear_mag_DEV,          Shear_mag_HOST->data);
    turbotmp::copy_array4_to_FortranHost(vert_vort_mag_DEV,      vert_vort_mag_HOST->data);
    turbotmp::copy_array4_to_FortranHost(hrat_min_DEV,           hrat_min_HOST->data);
    turbotmp::copy_array4_to_FortranHost(visc_bound_rem_DEV,     visc_bound_rem_HOST->data);

    /// Free everything we allocated (freeing a never-allocated container is a no-op)
    turbotmp::free_array4(u_DEV);
    turbotmp::free_array4(v_DEV);
    turbotmp::free_array4(h_DEV);
    turbotmp::free_array4(Kh_h_DEV);
    turbotmp::free_array4(grid_Re_Kh_DEV);
    turbotmp::free_array4(div_xx_h_DEV);
    turbotmp::free_array4(sh_xx_h_DEV);
    turbotmp::free_array4(h_u_DEV);
    turbotmp::free_array4(h_v_DEV);
    turbotmp::free_array4(sh_xx_DEV);
    turbotmp::free_array4(str_xx_DEV);
    turbotmp::free_array4(grad_vort_mag_h_DEV);
    turbotmp::free_array4(grad_vort_mag_h_2d_DEV);
    turbotmp::free_array4(grad_div_mag_h_DEV);
    turbotmp::free_array4(dudx_DEV);
    turbotmp::free_array4(dvdy_DEV);
    turbotmp::free_array4(sh_xy_DEV);
    turbotmp::free_array4(Del2vort_q_DEV);
    turbotmp::free_array4(Ah_DEV);
    turbotmp::free_array4(Kh_DEV);
    turbotmp::free_array4(Shear_mag_DEV);
    turbotmp::free_array4(vert_vort_mag_DEV);
    turbotmp::free_array4(hrat_min_DEV);
    turbotmp::free_array4(visc_bound_rem_DEV);
    turbotmp::free_array4(MEKE_Ku_DEV);
    turbotmp::free_array4(Res_fn_h_DEV);
    turbotmp::free_array4(BS_struct_DEV);
    turbotmp::free_array4(Kh_bg_xx_DEV);
    turbotmp::free_array4(Ah_bg_xx_DEV);
    turbotmp::free_array4(Kh_Max_xx_DEV);
    turbotmp::free_array4(n1n2_h_DEV);
    turbotmp::free_array4(n1n1_m_n2n2_h_DEV);
    turbotmp::free_array4(grid_sp_h2_DEV);
    turbotmp::free_array4(Laplac2_const_xx_DEV);
    turbotmp::free_array4(Biharm6_const_xx_DEV);
    turbotmp::free_array4(Laplac3_const_xx_DEV);
    turbotmp::free_array4(Biharm_const_xx_DEV);
    turbotmp::free_array4(Biharm_const2_xx_DEV);
}

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
/// @param bxU_HOST               The u-point iteration box of this tile, [IscB:IecB, jsc:jec, 1:kmax]
/// @param bxV_HOST               The v-point iteration box of this tile, [isc:iec, JscB:JecB, 1:kmax]
/// @param bxHstr_HOST            The h points of this tile where the normal stresses are computed,
///                               [IscB:IecB+1, JscB:JecB+1, 1:kmax]
/// @param bxQstr_HOST            The q points of this tile where the shear stresses are computed,
///                               [isc-1:iec, jsc-1:jec, 1:kmax]
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
                                         const Box_C* bxU_HOST,
                                         const Box_C* bxV_HOST,
                                         const Box_C* bxHstr_HOST,
                                         const Box_C* bxQstr_HOST,
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
                                         const double s_to_T)
{
    (void)nkblock;

    /// Reconstruct the iteration Boxes (Fortran 1-based -> AMReX 0-based)
    Box bxH = turbotmp::to_box(bxH_HOST);
    Box bxU = turbotmp::to_box(bxU_HOST);
    Box bxV = turbotmp::to_box(bxV_HOST);
    Box bxHstr = turbotmp::to_box(bxHstr_HOST);
    Box bxQstr = turbotmp::to_box(bxQstr_HOST);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto mask2dBu_DEV          = turbotmp::make_array4(&Gcore_HOST->mask2dBu);
    auto IdxCu_DEV             = turbotmp::make_array4(&Gcore_HOST->IdxCu);
    auto IdyCu_DEV             = turbotmp::make_array4(&Gcore_HOST->IdyCu);
    auto IareaCu_DEV           = turbotmp::make_array4(&Gcore_HOST->IareaCu);
    auto IdyCv_DEV             = turbotmp::make_array4(&Gcore_HOST->IdyCv);
    auto IdxCv_DEV             = turbotmp::make_array4(&Gcore_HOST->IdxCv);
    auto IareaCv_DEV           = turbotmp::make_array4(&Gcore_HOST->IareaCv);
    auto IdxT_DEV              = turbotmp::make_array4(&Gcore_HOST->IdxT);
    auto IdyT_DEV              = turbotmp::make_array4(&Gcore_HOST->IdyT);
    auto IdyBu_DEV             = turbotmp::make_array4(&Gcore_HOST->IdyBu);
    auto IdxBu_DEV             = turbotmp::make_array4(&Gcore_HOST->IdxBu);
    auto IareaT_DEV            = turbotmp::make_array4(&Gcore_HOST->IareaT);
    auto dxCu_DEV              = turbotmp::make_array4(&Gcore_HOST->dxCu);
    auto dyCv_DEV              = turbotmp::make_array4(&Gcore_HOST->dyCv);
    auto CoriolisBu_DEV        = turbotmp::make_array4(&Gcore_HOST->CoriolisBu);
    auto u_DEV                 = turbotmp::make_array4(u_HOST);
    auto v_DEV                 = turbotmp::make_array4(v_HOST);
    auto h_DEV                 = turbotmp::make_array4(h_HOST);
    auto uh_DEV                = turbotmp::make_array4(uh_HOST);
    auto vh_DEV                = turbotmp::make_array4(vh_HOST);
    auto diffu_DEV             = turbotmp::make_array4(diffu_HOST);
    auto diffv_DEV             = turbotmp::make_array4(diffv_HOST);
    auto str_xx_GME_DEV        = turbotmp::make_array4(str_xx_GME_HOST);
    auto str_xy_GME_DEV        = turbotmp::make_array4(str_xy_GME_HOST);
    auto FrictWork_DEV         = turbotmp::make_array4(FrictWork_HOST);
    auto FrictWork_bh_DEV      = turbotmp::make_array4(FrictWork_bh_HOST);
    auto FrictWork_GME_DEV     = turbotmp::make_array4(FrictWork_GME_HOST);
    auto visc_limit_h_flag_DEV = turbotmp::make_array4(visc_limit_h_flag_HOST);
    auto h_u_DEV               = turbotmp::make_array4(h_u_HOST);
    auto h_v_DEV               = turbotmp::make_array4(h_v_HOST);
    auto sh_xx_DEV             = turbotmp::make_array4(sh_xx_HOST);
    auto str_xx_DEV            = turbotmp::make_array4(str_xx_HOST);
    auto bhstr_xx_DEV          = turbotmp::make_array4(bhstr_xx_HOST);
    auto sh_xy_DEV             = turbotmp::make_array4(sh_xy_HOST);
    auto str_xy_DEV            = turbotmp::make_array4(str_xy_HOST);
    auto bhstr_xy_DEV          = turbotmp::make_array4(bhstr_xy_HOST);
    auto hq_DEV                = turbotmp::make_array4(hq_HOST);

    /// The MEKE, STOCH and metric arrays may be unallocated (data == nullptr); only allocate/copy them
    /// when allocated.
    const bool MEKE_mom_src_allocated = (MEKE_HOST->mom_src.data != nullptr);
    const bool MEKE_mom_src_bh_allocated = (MEKE_HOST->mom_src_bh.data != nullptr);
    const bool MEKE_GME_snk_allocated = (MEKE_HOST->GME_snk.data != nullptr);
    const bool skeb_diss_allocated = (STOCH_HOST->skeb_diss.data != nullptr);
    const bool reduction_xx_allocated = (metrics_HOST->reduction_xx.data != nullptr);
    const bool reduction_xy_allocated = (metrics_HOST->reduction_xy.data != nullptr);
    const bool dx2h_allocated = (metrics_HOST->dx2h.data != nullptr);
    const bool dy2h_allocated = (metrics_HOST->dy2h.data != nullptr);
    const bool dx2q_allocated = (metrics_HOST->dx2q.data != nullptr);
    const bool dy2q_allocated = (metrics_HOST->dy2q.data != nullptr);
    turbotmp::A4Box MEKE_mom_src_DEV{};
    turbotmp::A4Box MEKE_mom_src_bh_DEV{};
    turbotmp::A4Box MEKE_GME_snk_DEV{};
    turbotmp::A4Box skeb_diss_DEV{};
    turbotmp::A4Box reduction_xx_DEV{};
    turbotmp::A4Box reduction_xy_DEV{};
    turbotmp::A4Box dx2h_DEV{};
    turbotmp::A4Box dy2h_DEV{};
    turbotmp::A4Box dx2q_DEV{};
    turbotmp::A4Box dy2q_DEV{};
    if (MEKE_mom_src_allocated) MEKE_mom_src_DEV = turbotmp::make_array4(&MEKE_HOST->mom_src);
    if (MEKE_mom_src_bh_allocated) MEKE_mom_src_bh_DEV = turbotmp::make_array4(&MEKE_HOST->mom_src_bh);
    if (MEKE_GME_snk_allocated) MEKE_GME_snk_DEV = turbotmp::make_array4(&MEKE_HOST->GME_snk);
    if (skeb_diss_allocated) skeb_diss_DEV = turbotmp::make_array4(&STOCH_HOST->skeb_diss);
    if (reduction_xx_allocated) reduction_xx_DEV = turbotmp::make_array4(&metrics_HOST->reduction_xx);
    if (reduction_xy_allocated) reduction_xy_DEV = turbotmp::make_array4(&metrics_HOST->reduction_xy);
    if (dx2h_allocated) dx2h_DEV = turbotmp::make_array4(&metrics_HOST->dx2h);
    if (dy2h_allocated) dy2h_DEV = turbotmp::make_array4(&metrics_HOST->dy2h);
    if (dx2q_allocated) dx2q_DEV = turbotmp::make_array4(&metrics_HOST->dx2q);
    if (dy2q_allocated) dy2q_DEV = turbotmp::make_array4(&metrics_HOST->dy2q);

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dBu.data,    mask2dBu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdxCu.data,       IdxCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdyCu.data,       IdyCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IareaCu.data,     IareaCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdyCv.data,       IdyCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdxCv.data,       IdxCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IareaCv.data,     IareaCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdxT.data,        IdxT_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdyT.data,        IdyT_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdyBu.data,       IdyBu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdxBu.data,       IdxBu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IareaT.data,      IareaT_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->dxCu.data,        dxCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->dyCv.data,        dyCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->CoriolisBu.data,  CoriolisBu_DEV);
    turbotmp::copy_FortranHost_to_array4(u_HOST->data,                 u_DEV);
    turbotmp::copy_FortranHost_to_array4(v_HOST->data,                 v_DEV);
    turbotmp::copy_FortranHost_to_array4(h_HOST->data,                 h_DEV);
    turbotmp::copy_FortranHost_to_array4(uh_HOST->data,                uh_DEV);
    turbotmp::copy_FortranHost_to_array4(vh_HOST->data,                vh_DEV);
    turbotmp::copy_FortranHost_to_array4(diffu_HOST->data,             diffu_DEV);
    turbotmp::copy_FortranHost_to_array4(diffv_HOST->data,             diffv_DEV);
    turbotmp::copy_FortranHost_to_array4(str_xx_GME_HOST->data,        str_xx_GME_DEV);
    turbotmp::copy_FortranHost_to_array4(str_xy_GME_HOST->data,        str_xy_GME_DEV);
    turbotmp::copy_FortranHost_to_array4(FrictWork_HOST->data,         FrictWork_DEV);
    turbotmp::copy_FortranHost_to_array4(FrictWork_bh_HOST->data,      FrictWork_bh_DEV);
    turbotmp::copy_FortranHost_to_array4(FrictWork_GME_HOST->data,     FrictWork_GME_DEV);
    turbotmp::copy_FortranHost_to_array4(visc_limit_h_flag_HOST->data, visc_limit_h_flag_DEV);
    turbotmp::copy_FortranHost_to_array4(h_u_HOST->data,               h_u_DEV);
    turbotmp::copy_FortranHost_to_array4(h_v_HOST->data,               h_v_DEV);
    turbotmp::copy_FortranHost_to_array4(sh_xx_HOST->data,             sh_xx_DEV);
    turbotmp::copy_FortranHost_to_array4(str_xx_HOST->data,            str_xx_DEV);
    turbotmp::copy_FortranHost_to_array4(bhstr_xx_HOST->data,          bhstr_xx_DEV);
    turbotmp::copy_FortranHost_to_array4(sh_xy_HOST->data,             sh_xy_DEV);
    turbotmp::copy_FortranHost_to_array4(str_xy_HOST->data,            str_xy_DEV);
    turbotmp::copy_FortranHost_to_array4(bhstr_xy_HOST->data,          bhstr_xy_DEV);
    turbotmp::copy_FortranHost_to_array4(hq_HOST->data,                hq_DEV);
    if (MEKE_mom_src_allocated) turbotmp::copy_FortranHost_to_array4(MEKE_HOST->mom_src.data, MEKE_mom_src_DEV);
    if (MEKE_mom_src_bh_allocated)
        turbotmp::copy_FortranHost_to_array4(MEKE_HOST->mom_src_bh.data, MEKE_mom_src_bh_DEV);
    if (MEKE_GME_snk_allocated) turbotmp::copy_FortranHost_to_array4(MEKE_HOST->GME_snk.data, MEKE_GME_snk_DEV);
    if (skeb_diss_allocated) turbotmp::copy_FortranHost_to_array4(STOCH_HOST->skeb_diss.data, skeb_diss_DEV);
    if (reduction_xx_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->reduction_xx.data, reduction_xx_DEV);
    if (reduction_xy_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->reduction_xy.data, reduction_xy_DEV);
    if (dx2h_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->dx2h.data, dx2h_DEV);
    if (dy2h_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->dy2h.data, dy2h_DEV);
    if (dx2q_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->dx2q.data, dx2q_DEV);
    if (dy2q_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->dy2q.data, dy2q_DEV);

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::hor_visc_stress_div(bxH, bxU, bxV, bxHstr, bxQstr, *CS_HOST,
                             MEKE_mom_src_DEV.arr, MEKE_mom_src_bh_DEV.arr, MEKE_GME_snk_DEV.arr,
                             MEKE_HOST->backscatter_Ro_pow, MEKE_HOST->backscatter_Ro_c,
                             skeb_diss_DEV.arr, STOCH_HOST->skeb_frict_coef,
                             mask2dBu_DEV.arr, IdxCu_DEV.arr, IdyCu_DEV.arr, IareaCu_DEV.arr, IdyCv_DEV.arr,
                             IdxCv_DEV.arr, IareaCv_DEV.arr, IdxT_DEV.arr, IdyT_DEV.arr, IdyBu_DEV.arr,
                             IdxBu_DEV.arr, IareaT_DEV.arr, dxCu_DEV.arr, dyCv_DEV.arr, CoriolisBu_DEV.arr,
                             reduction_xx_DEV.arr, reduction_xy_DEV.arr, dx2h_DEV.arr,
                             dy2h_DEV.arr, dx2q_DEV.arr, dy2q_DEV.arr,
                             OBC,
                             kstart, h_neglect, apply_OBC, find_FrictWork, find_FrictWork_bh, backscat_subround,
                             skeb_use_frict,
                             u_DEV.arr, v_DEV.arr, h_DEV.arr, uh_DEV.arr, vh_DEV.arr,
                             diffu_DEV.arr, diffv_DEV.arr, str_xx_GME_DEV.arr, str_xy_GME_DEV.arr,
                             FrictWork_DEV.arr, FrictWork_bh_DEV.arr, FrictWork_GME_DEV.arr, visc_limit_h_flag_DEV.arr,
                             h_u_DEV.arr, h_v_DEV.arr, sh_xx_DEV.arr, str_xx_DEV.arr,
                             bhstr_xx_DEV.arr, sh_xy_DEV.arr, str_xy_DEV.arr, bhstr_xy_DEV.arr,
                             hq_DEV.arr,
                             H_to_RZ, s_to_T);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(diffu_DEV,             diffu_HOST->data);
    turbotmp::copy_array4_to_FortranHost(diffv_DEV,             diffv_HOST->data);
    turbotmp::copy_array4_to_FortranHost(str_xx_GME_DEV,        str_xx_GME_HOST->data);
    turbotmp::copy_array4_to_FortranHost(str_xy_GME_DEV,        str_xy_GME_HOST->data);
    turbotmp::copy_array4_to_FortranHost(FrictWork_DEV,         FrictWork_HOST->data);
    turbotmp::copy_array4_to_FortranHost(FrictWork_bh_DEV,      FrictWork_bh_HOST->data);
    turbotmp::copy_array4_to_FortranHost(FrictWork_GME_DEV,     FrictWork_GME_HOST->data);
    turbotmp::copy_array4_to_FortranHost(visc_limit_h_flag_DEV, visc_limit_h_flag_HOST->data);
    turbotmp::copy_array4_to_FortranHost(h_u_DEV,               h_u_HOST->data);
    turbotmp::copy_array4_to_FortranHost(h_v_DEV,               h_v_HOST->data);
    turbotmp::copy_array4_to_FortranHost(sh_xx_DEV,             sh_xx_HOST->data);
    turbotmp::copy_array4_to_FortranHost(str_xx_DEV,            str_xx_HOST->data);
    turbotmp::copy_array4_to_FortranHost(bhstr_xx_DEV,          bhstr_xx_HOST->data);
    turbotmp::copy_array4_to_FortranHost(sh_xy_DEV,             sh_xy_HOST->data);
    turbotmp::copy_array4_to_FortranHost(str_xy_DEV,            str_xy_HOST->data);
    turbotmp::copy_array4_to_FortranHost(bhstr_xy_DEV,          bhstr_xy_HOST->data);
    turbotmp::copy_array4_to_FortranHost(hq_DEV,                hq_HOST->data);
    if (MEKE_mom_src_allocated) turbotmp::copy_array4_to_FortranHost(MEKE_mom_src_DEV, MEKE_HOST->mom_src.data);
    if (MEKE_mom_src_bh_allocated)
        turbotmp::copy_array4_to_FortranHost(MEKE_mom_src_bh_DEV, MEKE_HOST->mom_src_bh.data);
    if (MEKE_GME_snk_allocated) turbotmp::copy_array4_to_FortranHost(MEKE_GME_snk_DEV, MEKE_HOST->GME_snk.data);
    if (skeb_diss_allocated) turbotmp::copy_array4_to_FortranHost(skeb_diss_DEV, STOCH_HOST->skeb_diss.data);

    /// Free everything we allocated (freeing a never-allocated container is a no-op)
    turbotmp::free_array4(mask2dBu_DEV);
    turbotmp::free_array4(IdxCu_DEV);
    turbotmp::free_array4(IdyCu_DEV);
    turbotmp::free_array4(IareaCu_DEV);
    turbotmp::free_array4(IdyCv_DEV);
    turbotmp::free_array4(IdxCv_DEV);
    turbotmp::free_array4(IareaCv_DEV);
    turbotmp::free_array4(IdxT_DEV);
    turbotmp::free_array4(IdyT_DEV);
    turbotmp::free_array4(IdyBu_DEV);
    turbotmp::free_array4(IdxBu_DEV);
    turbotmp::free_array4(IareaT_DEV);
    turbotmp::free_array4(dxCu_DEV);
    turbotmp::free_array4(dyCv_DEV);
    turbotmp::free_array4(CoriolisBu_DEV);
    turbotmp::free_array4(u_DEV);
    turbotmp::free_array4(v_DEV);
    turbotmp::free_array4(h_DEV);
    turbotmp::free_array4(uh_DEV);
    turbotmp::free_array4(vh_DEV);
    turbotmp::free_array4(diffu_DEV);
    turbotmp::free_array4(diffv_DEV);
    turbotmp::free_array4(str_xx_GME_DEV);
    turbotmp::free_array4(str_xy_GME_DEV);
    turbotmp::free_array4(FrictWork_DEV);
    turbotmp::free_array4(FrictWork_bh_DEV);
    turbotmp::free_array4(FrictWork_GME_DEV);
    turbotmp::free_array4(visc_limit_h_flag_DEV);
    turbotmp::free_array4(h_u_DEV);
    turbotmp::free_array4(h_v_DEV);
    turbotmp::free_array4(sh_xx_DEV);
    turbotmp::free_array4(str_xx_DEV);
    turbotmp::free_array4(bhstr_xx_DEV);
    turbotmp::free_array4(sh_xy_DEV);
    turbotmp::free_array4(str_xy_DEV);
    turbotmp::free_array4(bhstr_xy_DEV);
    turbotmp::free_array4(hq_DEV);
    turbotmp::free_array4(MEKE_mom_src_DEV);
    turbotmp::free_array4(MEKE_mom_src_bh_DEV);
    turbotmp::free_array4(MEKE_GME_snk_DEV);
    turbotmp::free_array4(skeb_diss_DEV);
    turbotmp::free_array4(reduction_xx_DEV);
    turbotmp::free_array4(reduction_xy_DEV);
    turbotmp::free_array4(dx2h_DEV);
    turbotmp::free_array4(dy2h_DEV);
    turbotmp::free_array4(dx2q_DEV);
    turbotmp::free_array4(dy2q_DEV);
}

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
/// @param bxHstr_HOST            The h points of this tile where the normal stresses are computed,
///                               [IscB:IecB+1, JscB:JecB+1, 1:kmax]
/// @param bxQstr_HOST            The q points of this tile where the shear stresses are computed,
///                               [isc-1:iec, jsc-1:jec, 1:kmax]
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
                                         const Box_C* bxHstr_HOST,
                                         const Box_C* bxQstr_HOST,
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
                                         RealArray_C* visc_bound_rem_HOST)
{
    (void)nkblock;

    /// Reconstruct the iteration Boxes (Fortran 1-based -> AMReX 0-based)
    Box bxH = turbotmp::to_box(bxH_HOST);
    Box bxQ = turbotmp::to_box(bxQ_HOST);
    Box bxKh = turbotmp::to_box(bxKh_HOST);
    Box bxHstr = turbotmp::to_box(bxHstr_HOST);
    Box bxQstr = turbotmp::to_box(bxQstr_HOST);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto IdyCu_DEV             = turbotmp::make_array4(&Gcore_HOST->IdyCu);
    auto IdxCv_DEV             = turbotmp::make_array4(&Gcore_HOST->IdxCv);
    auto IdyCv_DEV             = turbotmp::make_array4(&Gcore_HOST->IdyCv);
    auto IdxCu_DEV             = turbotmp::make_array4(&Gcore_HOST->IdxCu);
    auto mask2dBu_DEV          = turbotmp::make_array4(&Gcore_HOST->mask2dBu);
    auto mask2dCu_DEV          = turbotmp::make_array4(&Gcore_HOST->mask2dCu);
    auto mask2dCv_DEV          = turbotmp::make_array4(&Gcore_HOST->mask2dCv);
    auto u_DEV                 = turbotmp::make_array4(u_HOST);
    auto v_DEV                 = turbotmp::make_array4(v_HOST);
    auto h_DEV                 = turbotmp::make_array4(h_HOST);
    auto Ah_h_DEV              = turbotmp::make_array4(Ah_h_HOST);
    auto Kh_h_DEV              = turbotmp::make_array4(Kh_h_HOST);
    auto grid_Re_Ah_DEV        = turbotmp::make_array4(grid_Re_Ah_HOST);
    auto visc_limit_h_DEV      = turbotmp::make_array4(visc_limit_h_HOST);
    auto visc_limit_h_flag_DEV = turbotmp::make_array4(visc_limit_h_flag_HOST);
    auto visc_limit_h_frac_DEV = turbotmp::make_array4(visc_limit_h_frac_HOST);
    auto BS_coeff_h_DEV        = turbotmp::make_array4(BS_coeff_h_HOST);
    auto Del2u_DEV             = turbotmp::make_array4(Del2u_HOST);
    auto h_u_DEV               = turbotmp::make_array4(h_u_HOST);
    auto Del2v_DEV             = turbotmp::make_array4(Del2v_HOST);
    auto h_v_DEV               = turbotmp::make_array4(h_v_HOST);
    auto sh_xx_DEV             = turbotmp::make_array4(sh_xx_HOST);
    auto sh_xx_smooth_DEV      = turbotmp::make_array4(sh_xx_smooth_HOST);
    auto str_xx_DEV            = turbotmp::make_array4(str_xx_HOST);
    auto bhstr_xx_DEV          = turbotmp::make_array4(bhstr_xx_HOST);
    auto m_leithy_DEV          = turbotmp::make_array4(m_leithy_HOST);
    auto dDel2vdx_DEV          = turbotmp::make_array4(dDel2vdx_HOST);
    auto dDel2udy_DEV          = turbotmp::make_array4(dDel2udy_HOST);
    auto sh_xy_DEV             = turbotmp::make_array4(sh_xy_HOST);
    auto hq_DEV                = turbotmp::make_array4(hq_HOST);
    auto Ah_DEV                = turbotmp::make_array4(Ah_HOST);
    auto Kh_DEV                = turbotmp::make_array4(Kh_HOST);
    auto Shear_mag_DEV         = turbotmp::make_array4(Shear_mag_HOST);
    auto hrat_min_DEV          = turbotmp::make_array4(hrat_min_HOST);
    auto visc_bound_rem_DEV    = turbotmp::make_array4(visc_bound_rem_HOST);

    /// The MEKE, VarMix and metric arrays may be unallocated (data == nullptr); only allocate/copy them when
    /// allocated.
    const bool MEKE_Au_allocated = (MEKE_HOST->Au.data != nullptr);
    const bool MEKE_Ku_allocated = (MEKE_HOST->Ku.data != nullptr);
    const bool BS_struct_allocated = (VarMix_HOST->BS_struct.data != nullptr);
    const bool reduction_xx_allocated = (metrics_HOST->reduction_xx.data != nullptr);
    const bool Ah_Max_xx_allocated = (metrics_HOST->Ah_Max_xx.data != nullptr);
    const bool Ah_Max_xx_KS_allocated = (metrics_HOST->Ah_Max_xx_KS.data != nullptr);
    const bool grid_sp_h3_allocated = (metrics_HOST->grid_sp_h3.data != nullptr);
    const bool dx_dyT_allocated = (metrics_HOST->dx_dyT.data != nullptr);
    const bool dy_dxT_allocated = (metrics_HOST->dy_dxT.data != nullptr);
    const bool dx_dyBu_allocated = (metrics_HOST->dx_dyBu.data != nullptr);
    const bool dy_dxBu_allocated = (metrics_HOST->dy_dxBu.data != nullptr);
    const bool Re_Ah_const_xx_allocated = (metrics_HOST->Re_Ah_const_xx.data != nullptr);
    turbotmp::A4Box MEKE_Au_DEV{};
    turbotmp::A4Box MEKE_Ku_DEV{};
    turbotmp::A4Box BS_struct_DEV{};
    turbotmp::A4Box reduction_xx_DEV{};
    turbotmp::A4Box Ah_Max_xx_DEV{};
    turbotmp::A4Box Ah_Max_xx_KS_DEV{};
    turbotmp::A4Box grid_sp_h3_DEV{};
    turbotmp::A4Box dx_dyT_DEV{};
    turbotmp::A4Box dy_dxT_DEV{};
    turbotmp::A4Box dx_dyBu_DEV{};
    turbotmp::A4Box dy_dxBu_DEV{};
    turbotmp::A4Box Re_Ah_const_xx_DEV{};
    if (MEKE_Au_allocated) MEKE_Au_DEV = turbotmp::make_array4(&MEKE_HOST->Au);
    if (MEKE_Ku_allocated) MEKE_Ku_DEV = turbotmp::make_array4(&MEKE_HOST->Ku);
    if (BS_struct_allocated) BS_struct_DEV = turbotmp::make_array4(&VarMix_HOST->BS_struct);
    if (reduction_xx_allocated) reduction_xx_DEV = turbotmp::make_array4(&metrics_HOST->reduction_xx);
    if (Ah_Max_xx_allocated) Ah_Max_xx_DEV = turbotmp::make_array4(&metrics_HOST->Ah_Max_xx);
    if (Ah_Max_xx_KS_allocated) Ah_Max_xx_KS_DEV = turbotmp::make_array4(&metrics_HOST->Ah_Max_xx_KS);
    if (grid_sp_h3_allocated) grid_sp_h3_DEV = turbotmp::make_array4(&metrics_HOST->grid_sp_h3);
    if (dx_dyT_allocated) dx_dyT_DEV = turbotmp::make_array4(&metrics_HOST->dx_dyT);
    if (dy_dxT_allocated) dy_dxT_DEV = turbotmp::make_array4(&metrics_HOST->dy_dxT);
    if (dx_dyBu_allocated) dx_dyBu_DEV = turbotmp::make_array4(&metrics_HOST->dx_dyBu);
    if (dy_dxBu_allocated) dy_dxBu_DEV = turbotmp::make_array4(&metrics_HOST->dy_dxBu);
    if (Re_Ah_const_xx_allocated) Re_Ah_const_xx_DEV = turbotmp::make_array4(&metrics_HOST->Re_Ah_const_xx);

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdyCu.data,       IdyCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdxCv.data,       IdxCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdyCv.data,       IdyCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdxCu.data,       IdxCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dBu.data,    mask2dBu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dCu.data,    mask2dCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dCv.data,    mask2dCv_DEV);
    turbotmp::copy_FortranHost_to_array4(u_HOST->data,                 u_DEV);
    turbotmp::copy_FortranHost_to_array4(v_HOST->data,                 v_DEV);
    turbotmp::copy_FortranHost_to_array4(h_HOST->data,                 h_DEV);
    turbotmp::copy_FortranHost_to_array4(Ah_h_HOST->data,              Ah_h_DEV);
    turbotmp::copy_FortranHost_to_array4(Kh_h_HOST->data,              Kh_h_DEV);
    turbotmp::copy_FortranHost_to_array4(grid_Re_Ah_HOST->data,        grid_Re_Ah_DEV);
    turbotmp::copy_FortranHost_to_array4(visc_limit_h_HOST->data,      visc_limit_h_DEV);
    turbotmp::copy_FortranHost_to_array4(visc_limit_h_flag_HOST->data, visc_limit_h_flag_DEV);
    turbotmp::copy_FortranHost_to_array4(visc_limit_h_frac_HOST->data, visc_limit_h_frac_DEV);
    turbotmp::copy_FortranHost_to_array4(BS_coeff_h_HOST->data,        BS_coeff_h_DEV);
    turbotmp::copy_FortranHost_to_array4(Del2u_HOST->data,             Del2u_DEV);
    turbotmp::copy_FortranHost_to_array4(h_u_HOST->data,               h_u_DEV);
    turbotmp::copy_FortranHost_to_array4(Del2v_HOST->data,             Del2v_DEV);
    turbotmp::copy_FortranHost_to_array4(h_v_HOST->data,               h_v_DEV);
    turbotmp::copy_FortranHost_to_array4(sh_xx_HOST->data,             sh_xx_DEV);
    turbotmp::copy_FortranHost_to_array4(sh_xx_smooth_HOST->data,      sh_xx_smooth_DEV);
    turbotmp::copy_FortranHost_to_array4(str_xx_HOST->data,            str_xx_DEV);
    turbotmp::copy_FortranHost_to_array4(bhstr_xx_HOST->data,          bhstr_xx_DEV);
    turbotmp::copy_FortranHost_to_array4(m_leithy_HOST->data,          m_leithy_DEV);
    turbotmp::copy_FortranHost_to_array4(dDel2vdx_HOST->data,          dDel2vdx_DEV);
    turbotmp::copy_FortranHost_to_array4(dDel2udy_HOST->data,          dDel2udy_DEV);
    turbotmp::copy_FortranHost_to_array4(sh_xy_HOST->data,             sh_xy_DEV);
    turbotmp::copy_FortranHost_to_array4(hq_HOST->data,                hq_DEV);
    turbotmp::copy_FortranHost_to_array4(Ah_HOST->data,                Ah_DEV);
    turbotmp::copy_FortranHost_to_array4(Kh_HOST->data,                Kh_DEV);
    turbotmp::copy_FortranHost_to_array4(Shear_mag_HOST->data,         Shear_mag_DEV);
    turbotmp::copy_FortranHost_to_array4(hrat_min_HOST->data,          hrat_min_DEV);
    turbotmp::copy_FortranHost_to_array4(visc_bound_rem_HOST->data,    visc_bound_rem_DEV);
    if (MEKE_Au_allocated) turbotmp::copy_FortranHost_to_array4(MEKE_HOST->Au.data, MEKE_Au_DEV);
    if (MEKE_Ku_allocated) turbotmp::copy_FortranHost_to_array4(MEKE_HOST->Ku.data, MEKE_Ku_DEV);
    if (BS_struct_allocated) turbotmp::copy_FortranHost_to_array4(VarMix_HOST->BS_struct.data, BS_struct_DEV);
    if (reduction_xx_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->reduction_xx.data, reduction_xx_DEV);
    if (Ah_Max_xx_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Ah_Max_xx.data, Ah_Max_xx_DEV);
    if (Ah_Max_xx_KS_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Ah_Max_xx_KS.data, Ah_Max_xx_KS_DEV);
    if (grid_sp_h3_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->grid_sp_h3.data, grid_sp_h3_DEV);
    if (dx_dyT_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->dx_dyT.data, dx_dyT_DEV);
    if (dy_dxT_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->dy_dxT.data, dy_dxT_DEV);
    if (dx_dyBu_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->dx_dyBu.data, dx_dyBu_DEV);
    if (dy_dxBu_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->dy_dxBu.data, dy_dxBu_DEV);
    if (Re_Ah_const_xx_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->Re_Ah_const_xx.data, Re_Ah_const_xx_DEV);

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::hor_visc_h_coeffs_b(bxH, bxQ, bxKh, bxHstr, bxQstr, *CS_HOST,
                             MEKE_Au_DEV.arr, MEKE_Ku_DEV.arr, BS_struct_DEV.arr,
                             IdyCu_DEV.arr, IdxCv_DEV.arr, IdyCv_DEV.arr, IdxCu_DEV.arr,
                             mask2dBu_DEV.arr, mask2dCu_DEV.arr, mask2dCv_DEV.arr,
                             reduction_xx_DEV.arr, Ah_Max_xx_DEV.arr, Ah_Max_xx_KS_DEV.arr,
                             grid_sp_h3_DEV.arr, dx_dyT_DEV.arr, dy_dxT_DEV.arr,
                             dx_dyBu_DEV.arr, dy_dxBu_DEV.arr, Re_Ah_const_xx_DEV.arr,
                             OBC,
                             kstart, h_neglect, h_neglect3,
                             apply_OBC, use_kh_struct, use_MEKE_Au,
                             u_DEV.arr, v_DEV.arr, h_DEV.arr,
                             Ah_h_DEV.arr, Kh_h_DEV.arr, grid_Re_Ah_DEV.arr, visc_limit_h_DEV.arr,
                             visc_limit_h_flag_DEV.arr, visc_limit_h_frac_DEV.arr, BS_coeff_h_DEV.arr, Del2u_DEV.arr,
                             h_u_DEV.arr, Del2v_DEV.arr, h_v_DEV.arr, sh_xx_DEV.arr,
                             sh_xx_smooth_DEV.arr, str_xx_DEV.arr, bhstr_xx_DEV.arr, m_leithy_DEV.arr,
                             dDel2vdx_DEV.arr, dDel2udy_DEV.arr, sh_xy_DEV.arr, hq_DEV.arr,
                             Ah_DEV.arr, Kh_DEV.arr, Shear_mag_DEV.arr, hrat_min_DEV.arr,
                             visc_bound_rem_DEV.arr);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(Ah_h_DEV,              Ah_h_HOST->data);
    turbotmp::copy_array4_to_FortranHost(Kh_h_DEV,              Kh_h_HOST->data);
    turbotmp::copy_array4_to_FortranHost(grid_Re_Ah_DEV,        grid_Re_Ah_HOST->data);
    turbotmp::copy_array4_to_FortranHost(visc_limit_h_DEV,      visc_limit_h_HOST->data);
    turbotmp::copy_array4_to_FortranHost(visc_limit_h_flag_DEV, visc_limit_h_flag_HOST->data);
    turbotmp::copy_array4_to_FortranHost(visc_limit_h_frac_DEV, visc_limit_h_frac_HOST->data);
    turbotmp::copy_array4_to_FortranHost(BS_coeff_h_DEV,        BS_coeff_h_HOST->data);
    turbotmp::copy_array4_to_FortranHost(Del2u_DEV,             Del2u_HOST->data);
    turbotmp::copy_array4_to_FortranHost(h_u_DEV,               h_u_HOST->data);
    turbotmp::copy_array4_to_FortranHost(Del2v_DEV,             Del2v_HOST->data);
    turbotmp::copy_array4_to_FortranHost(h_v_DEV,               h_v_HOST->data);
    turbotmp::copy_array4_to_FortranHost(sh_xx_DEV,             sh_xx_HOST->data);
    turbotmp::copy_array4_to_FortranHost(sh_xx_smooth_DEV,      sh_xx_smooth_HOST->data);
    turbotmp::copy_array4_to_FortranHost(str_xx_DEV,            str_xx_HOST->data);
    turbotmp::copy_array4_to_FortranHost(bhstr_xx_DEV,          bhstr_xx_HOST->data);
    turbotmp::copy_array4_to_FortranHost(m_leithy_DEV,          m_leithy_HOST->data);
    turbotmp::copy_array4_to_FortranHost(dDel2vdx_DEV,          dDel2vdx_HOST->data);
    turbotmp::copy_array4_to_FortranHost(dDel2udy_DEV,          dDel2udy_HOST->data);
    turbotmp::copy_array4_to_FortranHost(sh_xy_DEV,             sh_xy_HOST->data);
    turbotmp::copy_array4_to_FortranHost(hq_DEV,                hq_HOST->data);
    turbotmp::copy_array4_to_FortranHost(Ah_DEV,                Ah_HOST->data);
    turbotmp::copy_array4_to_FortranHost(Kh_DEV,                Kh_HOST->data);
    turbotmp::copy_array4_to_FortranHost(Shear_mag_DEV,         Shear_mag_HOST->data);
    turbotmp::copy_array4_to_FortranHost(hrat_min_DEV,          hrat_min_HOST->data);
    turbotmp::copy_array4_to_FortranHost(visc_bound_rem_DEV,    visc_bound_rem_HOST->data);

    /// Free everything we allocated (freeing a never-allocated container is a no-op)
    turbotmp::free_array4(IdyCu_DEV);
    turbotmp::free_array4(IdxCv_DEV);
    turbotmp::free_array4(IdyCv_DEV);
    turbotmp::free_array4(IdxCu_DEV);
    turbotmp::free_array4(mask2dBu_DEV);
    turbotmp::free_array4(mask2dCu_DEV);
    turbotmp::free_array4(mask2dCv_DEV);
    turbotmp::free_array4(u_DEV);
    turbotmp::free_array4(v_DEV);
    turbotmp::free_array4(h_DEV);
    turbotmp::free_array4(Ah_h_DEV);
    turbotmp::free_array4(Kh_h_DEV);
    turbotmp::free_array4(grid_Re_Ah_DEV);
    turbotmp::free_array4(visc_limit_h_DEV);
    turbotmp::free_array4(visc_limit_h_flag_DEV);
    turbotmp::free_array4(visc_limit_h_frac_DEV);
    turbotmp::free_array4(BS_coeff_h_DEV);
    turbotmp::free_array4(Del2u_DEV);
    turbotmp::free_array4(h_u_DEV);
    turbotmp::free_array4(Del2v_DEV);
    turbotmp::free_array4(h_v_DEV);
    turbotmp::free_array4(sh_xx_DEV);
    turbotmp::free_array4(sh_xx_smooth_DEV);
    turbotmp::free_array4(str_xx_DEV);
    turbotmp::free_array4(bhstr_xx_DEV);
    turbotmp::free_array4(m_leithy_DEV);
    turbotmp::free_array4(dDel2vdx_DEV);
    turbotmp::free_array4(dDel2udy_DEV);
    turbotmp::free_array4(sh_xy_DEV);
    turbotmp::free_array4(hq_DEV);
    turbotmp::free_array4(Ah_DEV);
    turbotmp::free_array4(Kh_DEV);
    turbotmp::free_array4(Shear_mag_DEV);
    turbotmp::free_array4(hrat_min_DEV);
    turbotmp::free_array4(visc_bound_rem_DEV);
    turbotmp::free_array4(MEKE_Au_DEV);
    turbotmp::free_array4(MEKE_Ku_DEV);
    turbotmp::free_array4(BS_struct_DEV);
    turbotmp::free_array4(reduction_xx_DEV);
    turbotmp::free_array4(Ah_Max_xx_DEV);
    turbotmp::free_array4(Ah_Max_xx_KS_DEV);
    turbotmp::free_array4(grid_sp_h3_DEV);
    turbotmp::free_array4(dx_dyT_DEV);
    turbotmp::free_array4(dy_dxT_DEV);
    turbotmp::free_array4(dx_dyBu_DEV);
    turbotmp::free_array4(dy_dxBu_DEV);
    turbotmp::free_array4(Re_Ah_const_xx_DEV);
}

/// @brief Bridge for the MOM6 hor_visc_q_coeffs subroutine.
///
/// @param MEKE_HOST               MEKE fields related to Mesoscale Eddy Kinetic Energy; each field may be unallocated
///                                (data == nullptr)
/// @param VarMix_HOST             Variable mixing control structure; each field may be unallocated (data == nullptr)
/// @param Gcore_HOST              Persistent copies of the core ocean grid fields
/// @param CS_HOST                 Options used by the horizontal viscosity kernels
/// @param metrics_HOST            Persistent copies of the horizontal viscosity metric and coefficient arrays; each
///                                field may be unallocated (data == nullptr)
/// @param bxQstr_HOST             The q points of this tile where the shear stresses are computed,
///                                [isc-1:iec, jsc-1:jec, 1:kmax]
/// @param bxHgrad_HOST            The h points of this tile where the velocity gradients are computed,
///                                [IscB-1:IecB+2, JscB-1:JecB+2, 1:kmax]
/// @param bxQshear_HOST           The q points of this tile where the shearing strain is computed,
///                                [isc-2:iec+1, jsc-2:jec+1, 1:kmax]
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
                                       const Box_C* bxQstr_HOST,
                                       const Box_C* bxHgrad_HOST,
                                       const Box_C* bxQshear_HOST,
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
                                       RealArray_C* visc_bound_rem_HOST)
{
    (void)nkblock;
    (void)str_xx_HOST;

    /// Reconstruct the iteration Boxes (Fortran 1-based -> AMReX 0-based)
    Box bxQstr = turbotmp::to_box(bxQstr_HOST);
    Box bxHgrad = turbotmp::to_box(bxHgrad_HOST);
    Box bxQshear = turbotmp::to_box(bxQshear_HOST);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto mask2dBu_DEV           = turbotmp::make_array4(&Gcore_HOST->mask2dBu);
    auto u_DEV                  = turbotmp::make_array4(u_HOST);
    auto v_DEV                  = turbotmp::make_array4(v_HOST);
    auto sh_xx_bt_DEV           = turbotmp::make_array4(sh_xx_bt_HOST);
    auto sh_xy_bt_DEV           = turbotmp::make_array4(sh_xy_bt_HOST);
    auto GME_effic_h_DEV        = turbotmp::make_array4(GME_effic_h_HOST);
    auto GME_effic_q_DEV        = turbotmp::make_array4(GME_effic_q_HOST);
    auto KH_u_GME_DEV           = turbotmp::make_array4(KH_u_GME_HOST);
    auto KH_v_GME_DEV           = turbotmp::make_array4(KH_v_GME_HOST);
    auto GME_coeff_h_DEV        = turbotmp::make_array4(GME_coeff_h_HOST);
    auto GME_coeff_q_DEV        = turbotmp::make_array4(GME_coeff_q_HOST);
    auto str_xx_GME_DEV         = turbotmp::make_array4(str_xx_GME_HOST);
    auto str_xy_GME_DEV         = turbotmp::make_array4(str_xy_GME_HOST);
    auto Ah_h_DEV               = turbotmp::make_array4(Ah_h_HOST);
    auto Ah_q_DEV               = turbotmp::make_array4(Ah_q_HOST);
    auto Kh_h_DEV               = turbotmp::make_array4(Kh_h_HOST);
    auto Kh_q_DEV               = turbotmp::make_array4(Kh_q_HOST);
    auto vort_xy_q_DEV          = turbotmp::make_array4(vort_xy_q_HOST);
    auto sh_xy_q_DEV            = turbotmp::make_array4(sh_xy_q_HOST);
    auto visc_limit_q_DEV       = turbotmp::make_array4(visc_limit_q_HOST);
    auto visc_limit_q_flag_DEV  = turbotmp::make_array4(visc_limit_q_flag_HOST);
    auto visc_limit_q_frac_DEV  = turbotmp::make_array4(visc_limit_q_frac_HOST);
    auto BS_coeff_q_DEV         = turbotmp::make_array4(BS_coeff_q_HOST);
    auto sh_xx_DEV              = turbotmp::make_array4(sh_xx_HOST);
    auto dDel2vdx_DEV           = turbotmp::make_array4(dDel2vdx_HOST);
    auto dDel2udy_DEV           = turbotmp::make_array4(dDel2udy_HOST);
    auto sh_xy_DEV              = turbotmp::make_array4(sh_xy_HOST);
    auto sh_xy_smooth_DEV       = turbotmp::make_array4(sh_xy_smooth_HOST);
    auto str_xy_DEV             = turbotmp::make_array4(str_xy_HOST);
    auto bhstr_xy_DEV           = turbotmp::make_array4(bhstr_xy_HOST);
    auto vort_xy_DEV            = turbotmp::make_array4(vort_xy_HOST);
    auto grad_vort_mag_q_DEV    = turbotmp::make_array4(grad_vort_mag_q_HOST);
    auto grad_vort_mag_q_2d_DEV = turbotmp::make_array4(grad_vort_mag_q_2d_HOST);
    auto Del2vort_q_DEV         = turbotmp::make_array4(Del2vort_q_HOST);
    auto grad_div_mag_q_DEV     = turbotmp::make_array4(grad_div_mag_q_HOST);
    auto hq_DEV                 = turbotmp::make_array4(hq_HOST);
    auto Ah_DEV                 = turbotmp::make_array4(Ah_HOST);
    auto Kh_DEV                 = turbotmp::make_array4(Kh_HOST);
    auto Shear_mag_DEV          = turbotmp::make_array4(Shear_mag_HOST);
    auto vert_vort_mag_DEV      = turbotmp::make_array4(vert_vort_mag_HOST);
    auto hrat_min_DEV           = turbotmp::make_array4(hrat_min_HOST);
    auto visc_bound_rem_DEV     = turbotmp::make_array4(visc_bound_rem_HOST);

    /// The MEKE, VarMix and metric arrays may be unallocated (data == nullptr); only allocate/copy them when
    /// allocated.
    const bool MEKE_Ku_allocated = (MEKE_HOST->Ku.data != nullptr);
    const bool MEKE_Au_allocated = (MEKE_HOST->Au.data != nullptr);
    const bool Res_fn_q_allocated = (VarMix_HOST->Res_fn_q.data != nullptr);
    const bool BS_struct_allocated = (VarMix_HOST->BS_struct.data != nullptr);
    const bool Kh_bg_xy_allocated = (metrics_HOST->Kh_bg_xy.data != nullptr);
    const bool Ah_bg_xy_allocated = (metrics_HOST->Ah_bg_xy.data != nullptr);
    const bool reduction_xy_allocated = (metrics_HOST->reduction_xy.data != nullptr);
    const bool Kh_Max_xy_allocated = (metrics_HOST->Kh_Max_xy.data != nullptr);
    const bool Ah_Max_xy_allocated = (metrics_HOST->Ah_Max_xy.data != nullptr);
    const bool Ah_Max_xy_KS_allocated = (metrics_HOST->Ah_Max_xy_KS.data != nullptr);
    const bool n1n2_q_allocated = (metrics_HOST->n1n2_q.data != nullptr);
    const bool n1n1_m_n2n2_q_allocated = (metrics_HOST->n1n1_m_n2n2_q.data != nullptr);
    const bool Laplac2_const_xy_allocated = (metrics_HOST->Laplac2_const_xy.data != nullptr);
    const bool Biharm6_const_xy_allocated = (metrics_HOST->Biharm6_const_xy.data != nullptr);
    const bool Laplac3_const_xy_allocated = (metrics_HOST->Laplac3_const_xy.data != nullptr);
    const bool Biharm_const_xy_allocated = (metrics_HOST->Biharm_const_xy.data != nullptr);
    const bool Biharm_const2_xy_allocated = (metrics_HOST->Biharm_const2_xy.data != nullptr);
    const bool Re_Ah_const_xy_allocated = (metrics_HOST->Re_Ah_const_xy.data != nullptr);
    turbotmp::A4Box MEKE_Ku_DEV{};
    turbotmp::A4Box MEKE_Au_DEV{};
    turbotmp::A4Box Res_fn_q_DEV{};
    turbotmp::A4Box BS_struct_DEV{};
    turbotmp::A4Box Kh_bg_xy_DEV{};
    turbotmp::A4Box Ah_bg_xy_DEV{};
    turbotmp::A4Box reduction_xy_DEV{};
    turbotmp::A4Box Kh_Max_xy_DEV{};
    turbotmp::A4Box Ah_Max_xy_DEV{};
    turbotmp::A4Box Ah_Max_xy_KS_DEV{};
    turbotmp::A4Box n1n2_q_DEV{};
    turbotmp::A4Box n1n1_m_n2n2_q_DEV{};
    turbotmp::A4Box Laplac2_const_xy_DEV{};
    turbotmp::A4Box Biharm6_const_xy_DEV{};
    turbotmp::A4Box Laplac3_const_xy_DEV{};
    turbotmp::A4Box Biharm_const_xy_DEV{};
    turbotmp::A4Box Biharm_const2_xy_DEV{};
    turbotmp::A4Box Re_Ah_const_xy_DEV{};
    if (MEKE_Ku_allocated) MEKE_Ku_DEV = turbotmp::make_array4(&MEKE_HOST->Ku);
    if (MEKE_Au_allocated) MEKE_Au_DEV = turbotmp::make_array4(&MEKE_HOST->Au);
    if (Res_fn_q_allocated) Res_fn_q_DEV = turbotmp::make_array4(&VarMix_HOST->Res_fn_q);
    if (BS_struct_allocated) BS_struct_DEV = turbotmp::make_array4(&VarMix_HOST->BS_struct);
    if (Kh_bg_xy_allocated) Kh_bg_xy_DEV = turbotmp::make_array4(&metrics_HOST->Kh_bg_xy);
    if (Ah_bg_xy_allocated) Ah_bg_xy_DEV = turbotmp::make_array4(&metrics_HOST->Ah_bg_xy);
    if (reduction_xy_allocated) reduction_xy_DEV = turbotmp::make_array4(&metrics_HOST->reduction_xy);
    if (Kh_Max_xy_allocated) Kh_Max_xy_DEV = turbotmp::make_array4(&metrics_HOST->Kh_Max_xy);
    if (Ah_Max_xy_allocated) Ah_Max_xy_DEV = turbotmp::make_array4(&metrics_HOST->Ah_Max_xy);
    if (Ah_Max_xy_KS_allocated) Ah_Max_xy_KS_DEV = turbotmp::make_array4(&metrics_HOST->Ah_Max_xy_KS);
    if (n1n2_q_allocated) n1n2_q_DEV = turbotmp::make_array4(&metrics_HOST->n1n2_q);
    if (n1n1_m_n2n2_q_allocated) n1n1_m_n2n2_q_DEV = turbotmp::make_array4(&metrics_HOST->n1n1_m_n2n2_q);
    if (Laplac2_const_xy_allocated) Laplac2_const_xy_DEV = turbotmp::make_array4(&metrics_HOST->Laplac2_const_xy);
    if (Biharm6_const_xy_allocated) Biharm6_const_xy_DEV = turbotmp::make_array4(&metrics_HOST->Biharm6_const_xy);
    if (Laplac3_const_xy_allocated) Laplac3_const_xy_DEV = turbotmp::make_array4(&metrics_HOST->Laplac3_const_xy);
    if (Biharm_const_xy_allocated) Biharm_const_xy_DEV = turbotmp::make_array4(&metrics_HOST->Biharm_const_xy);
    if (Biharm_const2_xy_allocated) Biharm_const2_xy_DEV = turbotmp::make_array4(&metrics_HOST->Biharm_const2_xy);
    if (Re_Ah_const_xy_allocated) Re_Ah_const_xy_DEV = turbotmp::make_array4(&metrics_HOST->Re_Ah_const_xy);

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dBu.data,     mask2dBu_DEV);
    turbotmp::copy_FortranHost_to_array4(u_HOST->data,                  u_DEV);
    turbotmp::copy_FortranHost_to_array4(v_HOST->data,                  v_DEV);
    turbotmp::copy_FortranHost_to_array4(sh_xx_bt_HOST->data,           sh_xx_bt_DEV);
    turbotmp::copy_FortranHost_to_array4(sh_xy_bt_HOST->data,           sh_xy_bt_DEV);
    turbotmp::copy_FortranHost_to_array4(GME_effic_h_HOST->data,        GME_effic_h_DEV);
    turbotmp::copy_FortranHost_to_array4(GME_effic_q_HOST->data,        GME_effic_q_DEV);
    turbotmp::copy_FortranHost_to_array4(KH_u_GME_HOST->data,           KH_u_GME_DEV);
    turbotmp::copy_FortranHost_to_array4(KH_v_GME_HOST->data,           KH_v_GME_DEV);
    turbotmp::copy_FortranHost_to_array4(GME_coeff_h_HOST->data,        GME_coeff_h_DEV);
    turbotmp::copy_FortranHost_to_array4(GME_coeff_q_HOST->data,        GME_coeff_q_DEV);
    turbotmp::copy_FortranHost_to_array4(str_xx_GME_HOST->data,         str_xx_GME_DEV);
    turbotmp::copy_FortranHost_to_array4(str_xy_GME_HOST->data,         str_xy_GME_DEV);
    turbotmp::copy_FortranHost_to_array4(Ah_h_HOST->data,               Ah_h_DEV);
    turbotmp::copy_FortranHost_to_array4(Ah_q_HOST->data,               Ah_q_DEV);
    turbotmp::copy_FortranHost_to_array4(Kh_h_HOST->data,               Kh_h_DEV);
    turbotmp::copy_FortranHost_to_array4(Kh_q_HOST->data,               Kh_q_DEV);
    turbotmp::copy_FortranHost_to_array4(vort_xy_q_HOST->data,          vort_xy_q_DEV);
    turbotmp::copy_FortranHost_to_array4(sh_xy_q_HOST->data,            sh_xy_q_DEV);
    turbotmp::copy_FortranHost_to_array4(visc_limit_q_HOST->data,       visc_limit_q_DEV);
    turbotmp::copy_FortranHost_to_array4(visc_limit_q_flag_HOST->data,  visc_limit_q_flag_DEV);
    turbotmp::copy_FortranHost_to_array4(visc_limit_q_frac_HOST->data,  visc_limit_q_frac_DEV);
    turbotmp::copy_FortranHost_to_array4(BS_coeff_q_HOST->data,         BS_coeff_q_DEV);
    turbotmp::copy_FortranHost_to_array4(sh_xx_HOST->data,              sh_xx_DEV);
    turbotmp::copy_FortranHost_to_array4(dDel2vdx_HOST->data,           dDel2vdx_DEV);
    turbotmp::copy_FortranHost_to_array4(dDel2udy_HOST->data,           dDel2udy_DEV);
    turbotmp::copy_FortranHost_to_array4(sh_xy_HOST->data,              sh_xy_DEV);
    turbotmp::copy_FortranHost_to_array4(sh_xy_smooth_HOST->data,       sh_xy_smooth_DEV);
    turbotmp::copy_FortranHost_to_array4(str_xy_HOST->data,             str_xy_DEV);
    turbotmp::copy_FortranHost_to_array4(bhstr_xy_HOST->data,           bhstr_xy_DEV);
    turbotmp::copy_FortranHost_to_array4(vort_xy_HOST->data,            vort_xy_DEV);
    turbotmp::copy_FortranHost_to_array4(grad_vort_mag_q_HOST->data,    grad_vort_mag_q_DEV);
    turbotmp::copy_FortranHost_to_array4(grad_vort_mag_q_2d_HOST->data, grad_vort_mag_q_2d_DEV);
    turbotmp::copy_FortranHost_to_array4(Del2vort_q_HOST->data,         Del2vort_q_DEV);
    turbotmp::copy_FortranHost_to_array4(grad_div_mag_q_HOST->data,     grad_div_mag_q_DEV);
    turbotmp::copy_FortranHost_to_array4(hq_HOST->data,                 hq_DEV);
    turbotmp::copy_FortranHost_to_array4(Ah_HOST->data,                 Ah_DEV);
    turbotmp::copy_FortranHost_to_array4(Kh_HOST->data,                 Kh_DEV);
    turbotmp::copy_FortranHost_to_array4(Shear_mag_HOST->data,          Shear_mag_DEV);
    turbotmp::copy_FortranHost_to_array4(vert_vort_mag_HOST->data,      vert_vort_mag_DEV);
    turbotmp::copy_FortranHost_to_array4(hrat_min_HOST->data,           hrat_min_DEV);
    turbotmp::copy_FortranHost_to_array4(visc_bound_rem_HOST->data,     visc_bound_rem_DEV);
    if (MEKE_Ku_allocated) turbotmp::copy_FortranHost_to_array4(MEKE_HOST->Ku.data, MEKE_Ku_DEV);
    if (MEKE_Au_allocated) turbotmp::copy_FortranHost_to_array4(MEKE_HOST->Au.data, MEKE_Au_DEV);
    if (Res_fn_q_allocated) turbotmp::copy_FortranHost_to_array4(VarMix_HOST->Res_fn_q.data, Res_fn_q_DEV);
    if (BS_struct_allocated) turbotmp::copy_FortranHost_to_array4(VarMix_HOST->BS_struct.data, BS_struct_DEV);
    if (Kh_bg_xy_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Kh_bg_xy.data, Kh_bg_xy_DEV);
    if (Ah_bg_xy_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Ah_bg_xy.data, Ah_bg_xy_DEV);
    if (reduction_xy_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->reduction_xy.data, reduction_xy_DEV);
    if (Kh_Max_xy_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Kh_Max_xy.data, Kh_Max_xy_DEV);
    if (Ah_Max_xy_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Ah_Max_xy.data, Ah_Max_xy_DEV);
    if (Ah_Max_xy_KS_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Ah_Max_xy_KS.data, Ah_Max_xy_KS_DEV);
    if (n1n2_q_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->n1n2_q.data, n1n2_q_DEV);
    if (n1n1_m_n2n2_q_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->n1n1_m_n2n2_q.data, n1n1_m_n2n2_q_DEV);
    if (Laplac2_const_xy_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->Laplac2_const_xy.data, Laplac2_const_xy_DEV);
    if (Biharm6_const_xy_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->Biharm6_const_xy.data, Biharm6_const_xy_DEV);
    if (Laplac3_const_xy_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->Laplac3_const_xy.data, Laplac3_const_xy_DEV);
    if (Biharm_const_xy_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->Biharm_const_xy.data, Biharm_const_xy_DEV);
    if (Biharm_const2_xy_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->Biharm_const2_xy.data, Biharm_const2_xy_DEV);
    if (Re_Ah_const_xy_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->Re_Ah_const_xy.data, Re_Ah_const_xy_DEV);

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::hor_visc_q_coeffs(bxQstr, bxHgrad, bxQshear, *CS_HOST, MEKE_Ku_DEV.arr, MEKE_Au_DEV.arr, Res_fn_q_DEV.arr,
                           BS_struct_DEV.arr,
                           mask2dBu_DEV.arr, Kh_bg_xy_DEV.arr, Ah_bg_xy_DEV.arr, reduction_xy_DEV.arr,
                           Kh_Max_xy_DEV.arr, Ah_Max_xy_DEV.arr, Ah_Max_xy_KS_DEV.arr, n1n2_q_DEV.arr,
                           n1n1_m_n2n2_q_DEV.arr, Laplac2_const_xy_DEV.arr, Biharm6_const_xy_DEV.arr,
                           Laplac3_const_xy_DEV.arr, Biharm_const_xy_DEV.arr, Biharm_const2_xy_DEV.arr,
                           Re_Ah_const_xy_DEV.arr, kstart, inv_PI3, inv_PI6, find_FrictWork, use_kh_struct,
                           use_MEKE_Ku, use_MEKE_Au, rescale_Kh, u_DEV.arr, v_DEV.arr, sh_xx_bt_DEV.arr,
                           sh_xy_bt_DEV.arr, GME_effic_h_DEV.arr, GME_effic_q_DEV.arr, KH_u_GME_DEV.arr,
                           KH_v_GME_DEV.arr, GME_coeff_h_DEV.arr, GME_coeff_q_DEV.arr, str_xx_GME_DEV.arr,
                           str_xy_GME_DEV.arr, Ah_h_DEV.arr, Ah_q_DEV.arr, Kh_h_DEV.arr, Kh_q_DEV.arr,
                           vort_xy_q_DEV.arr, sh_xy_q_DEV.arr, visc_limit_q_DEV.arr, visc_limit_q_flag_DEV.arr,
                           visc_limit_q_frac_DEV.arr, BS_coeff_q_DEV.arr, sh_xx_DEV.arr, dDel2vdx_DEV.arr,
                           dDel2udy_DEV.arr, sh_xy_DEV.arr, sh_xy_smooth_DEV.arr, str_xy_DEV.arr, bhstr_xy_DEV.arr,
                           vort_xy_DEV.arr, grad_vort_mag_q_DEV.arr, grad_vort_mag_q_2d_DEV.arr, Del2vort_q_DEV.arr,
                           grad_div_mag_q_DEV.arr, hq_DEV.arr, Ah_DEV.arr, Kh_DEV.arr, Shear_mag_DEV.arr,
                           vert_vort_mag_DEV.arr, hrat_min_DEV.arr, visc_bound_rem_DEV.arr);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(GME_coeff_h_DEV,        GME_coeff_h_HOST->data);
    turbotmp::copy_array4_to_FortranHost(GME_coeff_q_DEV,        GME_coeff_q_HOST->data);
    turbotmp::copy_array4_to_FortranHost(str_xx_GME_DEV,         str_xx_GME_HOST->data);
    turbotmp::copy_array4_to_FortranHost(str_xy_GME_DEV,         str_xy_GME_HOST->data);
    turbotmp::copy_array4_to_FortranHost(Ah_h_DEV,               Ah_h_HOST->data);
    turbotmp::copy_array4_to_FortranHost(Ah_q_DEV,               Ah_q_HOST->data);
    turbotmp::copy_array4_to_FortranHost(Kh_h_DEV,               Kh_h_HOST->data);
    turbotmp::copy_array4_to_FortranHost(Kh_q_DEV,               Kh_q_HOST->data);
    turbotmp::copy_array4_to_FortranHost(vort_xy_q_DEV,          vort_xy_q_HOST->data);
    turbotmp::copy_array4_to_FortranHost(sh_xy_q_DEV,            sh_xy_q_HOST->data);
    turbotmp::copy_array4_to_FortranHost(visc_limit_q_DEV,       visc_limit_q_HOST->data);
    turbotmp::copy_array4_to_FortranHost(visc_limit_q_flag_DEV,  visc_limit_q_flag_HOST->data);
    turbotmp::copy_array4_to_FortranHost(visc_limit_q_frac_DEV,  visc_limit_q_frac_HOST->data);
    turbotmp::copy_array4_to_FortranHost(BS_coeff_q_DEV,         BS_coeff_q_HOST->data);
    turbotmp::copy_array4_to_FortranHost(sh_xx_DEV,              sh_xx_HOST->data);
    turbotmp::copy_array4_to_FortranHost(dDel2vdx_DEV,           dDel2vdx_HOST->data);
    turbotmp::copy_array4_to_FortranHost(dDel2udy_DEV,           dDel2udy_HOST->data);
    turbotmp::copy_array4_to_FortranHost(sh_xy_DEV,              sh_xy_HOST->data);
    turbotmp::copy_array4_to_FortranHost(sh_xy_smooth_DEV,       sh_xy_smooth_HOST->data);
    turbotmp::copy_array4_to_FortranHost(str_xy_DEV,             str_xy_HOST->data);
    turbotmp::copy_array4_to_FortranHost(bhstr_xy_DEV,           bhstr_xy_HOST->data);
    turbotmp::copy_array4_to_FortranHost(vort_xy_DEV,            vort_xy_HOST->data);
    turbotmp::copy_array4_to_FortranHost(grad_vort_mag_q_DEV,    grad_vort_mag_q_HOST->data);
    turbotmp::copy_array4_to_FortranHost(grad_vort_mag_q_2d_DEV, grad_vort_mag_q_2d_HOST->data);
    turbotmp::copy_array4_to_FortranHost(Del2vort_q_DEV,         Del2vort_q_HOST->data);
    turbotmp::copy_array4_to_FortranHost(grad_div_mag_q_DEV,     grad_div_mag_q_HOST->data);
    turbotmp::copy_array4_to_FortranHost(hq_DEV,                 hq_HOST->data);
    turbotmp::copy_array4_to_FortranHost(Ah_DEV,                 Ah_HOST->data);
    turbotmp::copy_array4_to_FortranHost(Kh_DEV,                 Kh_HOST->data);
    turbotmp::copy_array4_to_FortranHost(Shear_mag_DEV,          Shear_mag_HOST->data);
    turbotmp::copy_array4_to_FortranHost(vert_vort_mag_DEV,      vert_vort_mag_HOST->data);
    turbotmp::copy_array4_to_FortranHost(hrat_min_DEV,           hrat_min_HOST->data);
    turbotmp::copy_array4_to_FortranHost(visc_bound_rem_DEV,     visc_bound_rem_HOST->data);

    /// Free everything we allocated (freeing a never-allocated container is a no-op)
    turbotmp::free_array4(mask2dBu_DEV);
    turbotmp::free_array4(u_DEV);
    turbotmp::free_array4(v_DEV);
    turbotmp::free_array4(sh_xx_bt_DEV);
    turbotmp::free_array4(sh_xy_bt_DEV);
    turbotmp::free_array4(GME_effic_h_DEV);
    turbotmp::free_array4(GME_effic_q_DEV);
    turbotmp::free_array4(KH_u_GME_DEV);
    turbotmp::free_array4(KH_v_GME_DEV);
    turbotmp::free_array4(GME_coeff_h_DEV);
    turbotmp::free_array4(GME_coeff_q_DEV);
    turbotmp::free_array4(str_xx_GME_DEV);
    turbotmp::free_array4(str_xy_GME_DEV);
    turbotmp::free_array4(Ah_h_DEV);
    turbotmp::free_array4(Ah_q_DEV);
    turbotmp::free_array4(Kh_h_DEV);
    turbotmp::free_array4(Kh_q_DEV);
    turbotmp::free_array4(vort_xy_q_DEV);
    turbotmp::free_array4(sh_xy_q_DEV);
    turbotmp::free_array4(visc_limit_q_DEV);
    turbotmp::free_array4(visc_limit_q_flag_DEV);
    turbotmp::free_array4(visc_limit_q_frac_DEV);
    turbotmp::free_array4(BS_coeff_q_DEV);
    turbotmp::free_array4(sh_xx_DEV);
    turbotmp::free_array4(dDel2vdx_DEV);
    turbotmp::free_array4(dDel2udy_DEV);
    turbotmp::free_array4(sh_xy_DEV);
    turbotmp::free_array4(sh_xy_smooth_DEV);
    turbotmp::free_array4(str_xy_DEV);
    turbotmp::free_array4(bhstr_xy_DEV);
    turbotmp::free_array4(vort_xy_DEV);
    turbotmp::free_array4(grad_vort_mag_q_DEV);
    turbotmp::free_array4(grad_vort_mag_q_2d_DEV);
    turbotmp::free_array4(Del2vort_q_DEV);
    turbotmp::free_array4(grad_div_mag_q_DEV);
    turbotmp::free_array4(hq_DEV);
    turbotmp::free_array4(Ah_DEV);
    turbotmp::free_array4(Kh_DEV);
    turbotmp::free_array4(Shear_mag_DEV);
    turbotmp::free_array4(vert_vort_mag_DEV);
    turbotmp::free_array4(hrat_min_DEV);
    turbotmp::free_array4(visc_bound_rem_DEV);
    turbotmp::free_array4(MEKE_Ku_DEV);
    turbotmp::free_array4(MEKE_Au_DEV);
    turbotmp::free_array4(Res_fn_q_DEV);
    turbotmp::free_array4(BS_struct_DEV);
    turbotmp::free_array4(Kh_bg_xy_DEV);
    turbotmp::free_array4(Ah_bg_xy_DEV);
    turbotmp::free_array4(reduction_xy_DEV);
    turbotmp::free_array4(Kh_Max_xy_DEV);
    turbotmp::free_array4(Ah_Max_xy_DEV);
    turbotmp::free_array4(Ah_Max_xy_KS_DEV);
    turbotmp::free_array4(n1n2_q_DEV);
    turbotmp::free_array4(n1n1_m_n2n2_q_DEV);
    turbotmp::free_array4(Laplac2_const_xy_DEV);
    turbotmp::free_array4(Biharm6_const_xy_DEV);
    turbotmp::free_array4(Laplac3_const_xy_DEV);
    turbotmp::free_array4(Biharm_const_xy_DEV);
    turbotmp::free_array4(Biharm_const2_xy_DEV);
    turbotmp::free_array4(Re_Ah_const_xy_DEV);
}

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
                                             void* host)
{
    (void)host;
    (void)dz_HOST;
    (void)slope_x_HOST;
    (void)slope_y_HOST;

    /// Reconstruct the iteration Boxes (Fortran 1-based -> AMReX 0-based)
    Box bxH0 = turbotmp::to_box(bxH0_HOST);
    Box bxQ0 = turbotmp::to_box(bxQ0_HOST);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto mask2dT_DEV           = turbotmp::make_array4(&Gcore_HOST->mask2dT);
    auto IareaT_DEV            = turbotmp::make_array4(&Gcore_HOST->IareaT);
    auto IdxT_DEV              = turbotmp::make_array4(&Gcore_HOST->IdxT);
    auto IdyT_DEV              = turbotmp::make_array4(&Gcore_HOST->IdyT);
    auto mask2dCu_DEV          = turbotmp::make_array4(&Gcore_HOST->mask2dCu);
    auto IdxCu_DEV             = turbotmp::make_array4(&Gcore_HOST->IdxCu);
    auto dxCu_DEV              = turbotmp::make_array4(&Gcore_HOST->dxCu);
    auto mask2dCv_DEV          = turbotmp::make_array4(&Gcore_HOST->mask2dCv);
    auto IdyCv_DEV             = turbotmp::make_array4(&Gcore_HOST->IdyCv);
    auto dyCv_DEV              = turbotmp::make_array4(&Gcore_HOST->dyCv);
    auto CoriolisBu_DEV        = turbotmp::make_array4(&Gcore_HOST->CoriolisBu);
    auto IdyCu_DEV             = turbotmp::make_array4(&Gcore_HOST->IdyCu);
    auto IdxCv_DEV             = turbotmp::make_array4(&Gcore_HOST->IdxCv);
    auto mask2dBu_DEV          = turbotmp::make_array4(&Gcore_HOST->mask2dBu);
    auto dxBu_DEV              = turbotmp::make_array4(&Gcore_HOST->dxBu);
    auto dyBu_DEV              = turbotmp::make_array4(&Gcore_HOST->dyBu);
    auto IareaBu_DEV           = turbotmp::make_array4(&Gcore_HOST->IareaBu);
    auto df_dx_DEV             = turbotmp::make_array4(&Gcore_HOST->df_dx);
    auto df_dy_DEV             = turbotmp::make_array4(&Gcore_HOST->df_dy);
    auto IareaCu_DEV           = turbotmp::make_array4(&Gcore_HOST->IareaCu);
    auto IareaCv_DEV           = turbotmp::make_array4(&Gcore_HOST->IareaCv);
    auto IdxBu_DEV             = turbotmp::make_array4(&Gcore_HOST->IdxBu);
    auto IdyBu_DEV             = turbotmp::make_array4(&Gcore_HOST->IdyBu);
    auto u_DEV                 = turbotmp::make_array4(u_HOST);
    auto v_DEV                 = turbotmp::make_array4(v_HOST);
    auto h_DEV                 = turbotmp::make_array4(h_HOST);
    auto uh_DEV                = turbotmp::make_array4(uh_HOST);
    auto vh_DEV                = turbotmp::make_array4(vh_HOST);
    auto diffu_DEV             = turbotmp::make_array4(diffu_HOST);
    auto diffv_DEV             = turbotmp::make_array4(diffv_HOST);
    auto u_smooth_DEV          = turbotmp::make_array4(u_smooth_HOST);
    auto v_smooth_DEV          = turbotmp::make_array4(v_smooth_HOST);
    auto sh_xx_bt_DEV          = turbotmp::make_array4(sh_xx_bt_HOST);
    auto sh_xy_bt_DEV          = turbotmp::make_array4(sh_xy_bt_HOST);
    auto GME_effic_h_DEV       = turbotmp::make_array4(GME_effic_h_HOST);
    auto GME_effic_q_DEV       = turbotmp::make_array4(GME_effic_q_HOST);
    auto KH_u_GME_DEV          = turbotmp::make_array4(KH_u_GME_HOST);
    auto KH_v_GME_DEV          = turbotmp::make_array4(KH_v_GME_HOST);
    auto GME_coeff_h_DEV       = turbotmp::make_array4(GME_coeff_h_HOST);
    auto GME_coeff_q_DEV       = turbotmp::make_array4(GME_coeff_q_HOST);
    auto str_xx_GME_DEV        = turbotmp::make_array4(str_xx_GME_HOST);
    auto str_xy_GME_DEV        = turbotmp::make_array4(str_xy_GME_HOST);
    auto NoSt_DEV              = turbotmp::make_array4(NoSt_HOST);
    auto ShSt_DEV              = turbotmp::make_array4(ShSt_HOST);
    auto FrictWork_DEV         = turbotmp::make_array4(FrictWork_HOST);
    auto FrictWork_bh_DEV      = turbotmp::make_array4(FrictWork_bh_HOST);
    auto FrictWork_GME_DEV     = turbotmp::make_array4(FrictWork_GME_HOST);
    auto Ah_h_DEV              = turbotmp::make_array4(Ah_h_HOST);
    auto Ah_q_DEV              = turbotmp::make_array4(Ah_q_HOST);
    auto Kh_h_DEV              = turbotmp::make_array4(Kh_h_HOST);
    auto Kh_q_DEV              = turbotmp::make_array4(Kh_q_HOST);
    auto grid_Re_Kh_DEV        = turbotmp::make_array4(grid_Re_Kh_HOST);
    auto grid_Re_Ah_DEV        = turbotmp::make_array4(grid_Re_Ah_HOST);
    auto div_xx_h_DEV          = turbotmp::make_array4(div_xx_h_HOST);
    auto sh_xx_h_DEV           = turbotmp::make_array4(sh_xx_h_HOST);
    auto vort_xy_q_DEV         = turbotmp::make_array4(vort_xy_q_HOST);
    auto sh_xy_q_DEV           = turbotmp::make_array4(sh_xy_q_HOST);
    auto visc_limit_h_DEV      = turbotmp::make_array4(visc_limit_h_HOST);
    auto visc_limit_q_DEV      = turbotmp::make_array4(visc_limit_q_HOST);
    auto visc_limit_h_flag_DEV = turbotmp::make_array4(visc_limit_h_flag_HOST);
    auto visc_limit_q_flag_DEV = turbotmp::make_array4(visc_limit_q_flag_HOST);
    auto visc_limit_h_frac_DEV = turbotmp::make_array4(visc_limit_h_frac_HOST);
    auto visc_limit_q_frac_DEV = turbotmp::make_array4(visc_limit_q_frac_HOST);
    auto BS_coeff_h_DEV        = turbotmp::make_array4(BS_coeff_h_HOST);
    auto BS_coeff_q_DEV        = turbotmp::make_array4(BS_coeff_q_HOST);

    /// The continuity thicknesses and the MEKE, VarMix, STOCH and metric arrays may be unallocated
    /// (data == nullptr); only allocate/copy them when allocated.
    const bool hu_cont_allocated = (hu_cont_HOST->data != nullptr);
    const bool hv_cont_allocated = (hv_cont_HOST->data != nullptr);
    const bool MEKE_Ku_allocated = (MEKE_HOST->Ku.data != nullptr);
    const bool MEKE_Au_allocated = (MEKE_HOST->Au.data != nullptr);
    const bool MEKE_mom_src_allocated = (MEKE_HOST->mom_src.data != nullptr);
    const bool MEKE_mom_src_bh_allocated = (MEKE_HOST->mom_src_bh.data != nullptr);
    const bool MEKE_GME_snk_allocated = (MEKE_HOST->GME_snk.data != nullptr);
    const bool Res_fn_h_allocated = (VarMix_HOST->Res_fn_h.data != nullptr);
    const bool Res_fn_q_allocated = (VarMix_HOST->Res_fn_q.data != nullptr);
    const bool BS_struct_allocated = (VarMix_HOST->BS_struct.data != nullptr);
    const bool skeb_diss_allocated = (STOCH_HOST->skeb_diss.data != nullptr);
    const bool Kh_bg_xx_allocated = (metrics_HOST->Kh_bg_xx.data != nullptr);
    const bool Ah_bg_xx_allocated = (metrics_HOST->Ah_bg_xx.data != nullptr);
    const bool reduction_xx_allocated = (metrics_HOST->reduction_xx.data != nullptr);
    const bool Kh_Max_xx_allocated = (metrics_HOST->Kh_Max_xx.data != nullptr);
    const bool Ah_Max_xx_allocated = (metrics_HOST->Ah_Max_xx.data != nullptr);
    const bool Ah_Max_xx_KS_allocated = (metrics_HOST->Ah_Max_xx_KS.data != nullptr);
    const bool n1n2_h_allocated = (metrics_HOST->n1n2_h.data != nullptr);
    const bool n1n1_m_n2n2_h_allocated = (metrics_HOST->n1n1_m_n2n2_h.data != nullptr);
    const bool grid_sp_h2_allocated = (metrics_HOST->grid_sp_h2.data != nullptr);
    const bool grid_sp_h3_allocated = (metrics_HOST->grid_sp_h3.data != nullptr);
    const bool Kh_bg_xy_allocated = (metrics_HOST->Kh_bg_xy.data != nullptr);
    const bool Ah_bg_xy_allocated = (metrics_HOST->Ah_bg_xy.data != nullptr);
    const bool reduction_xy_allocated = (metrics_HOST->reduction_xy.data != nullptr);
    const bool Kh_Max_xy_allocated = (metrics_HOST->Kh_Max_xy.data != nullptr);
    const bool Ah_Max_xy_allocated = (metrics_HOST->Ah_Max_xy.data != nullptr);
    const bool Ah_Max_xy_KS_allocated = (metrics_HOST->Ah_Max_xy_KS.data != nullptr);
    const bool n1n2_q_allocated = (metrics_HOST->n1n2_q.data != nullptr);
    const bool n1n1_m_n2n2_q_allocated = (metrics_HOST->n1n1_m_n2n2_q.data != nullptr);
    const bool dx2h_allocated = (metrics_HOST->dx2h.data != nullptr);
    const bool dy2h_allocated = (metrics_HOST->dy2h.data != nullptr);
    const bool dx_dyT_allocated = (metrics_HOST->dx_dyT.data != nullptr);
    const bool dy_dxT_allocated = (metrics_HOST->dy_dxT.data != nullptr);
    const bool m_const_leithy_allocated = (metrics_HOST->m_const_leithy.data != nullptr);
    const bool m_leithy_max_allocated = (metrics_HOST->m_leithy_max.data != nullptr);
    const bool dx2q_allocated = (metrics_HOST->dx2q.data != nullptr);
    const bool dy2q_allocated = (metrics_HOST->dy2q.data != nullptr);
    const bool dx_dyBu_allocated = (metrics_HOST->dx_dyBu.data != nullptr);
    const bool dy_dxBu_allocated = (metrics_HOST->dy_dxBu.data != nullptr);
    const bool Idx2dyCu_allocated = (metrics_HOST->Idx2dyCu.data != nullptr);
    const bool Idxdy2u_allocated = (metrics_HOST->Idxdy2u.data != nullptr);
    const bool Idx2dyCv_allocated = (metrics_HOST->Idx2dyCv.data != nullptr);
    const bool Idxdy2v_allocated = (metrics_HOST->Idxdy2v.data != nullptr);
    const bool Laplac2_const_xx_allocated = (metrics_HOST->Laplac2_const_xx.data != nullptr);
    const bool Biharm6_const_xx_allocated = (metrics_HOST->Biharm6_const_xx.data != nullptr);
    const bool Laplac3_const_xx_allocated = (metrics_HOST->Laplac3_const_xx.data != nullptr);
    const bool Biharm_const_xx_allocated = (metrics_HOST->Biharm_const_xx.data != nullptr);
    const bool Biharm_const2_xx_allocated = (metrics_HOST->Biharm_const2_xx.data != nullptr);
    const bool Re_Ah_const_xx_allocated = (metrics_HOST->Re_Ah_const_xx.data != nullptr);
    const bool Laplac2_const_xy_allocated = (metrics_HOST->Laplac2_const_xy.data != nullptr);
    const bool Biharm6_const_xy_allocated = (metrics_HOST->Biharm6_const_xy.data != nullptr);
    const bool Laplac3_const_xy_allocated = (metrics_HOST->Laplac3_const_xy.data != nullptr);
    const bool Biharm_const_xy_allocated = (metrics_HOST->Biharm_const_xy.data != nullptr);
    const bool Biharm_const2_xy_allocated = (metrics_HOST->Biharm_const2_xy.data != nullptr);
    const bool Re_Ah_const_xy_allocated = (metrics_HOST->Re_Ah_const_xy.data != nullptr);
    turbotmp::A4Box hu_cont_DEV{};
    turbotmp::A4Box hv_cont_DEV{};
    turbotmp::A4Box MEKE_Ku_DEV{};
    turbotmp::A4Box MEKE_Au_DEV{};
    turbotmp::A4Box MEKE_mom_src_DEV{};
    turbotmp::A4Box MEKE_mom_src_bh_DEV{};
    turbotmp::A4Box MEKE_GME_snk_DEV{};
    turbotmp::A4Box Res_fn_h_DEV{};
    turbotmp::A4Box Res_fn_q_DEV{};
    turbotmp::A4Box BS_struct_DEV{};
    turbotmp::A4Box skeb_diss_DEV{};
    turbotmp::A4Box Kh_bg_xx_DEV{};
    turbotmp::A4Box Ah_bg_xx_DEV{};
    turbotmp::A4Box reduction_xx_DEV{};
    turbotmp::A4Box Kh_Max_xx_DEV{};
    turbotmp::A4Box Ah_Max_xx_DEV{};
    turbotmp::A4Box Ah_Max_xx_KS_DEV{};
    turbotmp::A4Box n1n2_h_DEV{};
    turbotmp::A4Box n1n1_m_n2n2_h_DEV{};
    turbotmp::A4Box grid_sp_h2_DEV{};
    turbotmp::A4Box grid_sp_h3_DEV{};
    turbotmp::A4Box Kh_bg_xy_DEV{};
    turbotmp::A4Box Ah_bg_xy_DEV{};
    turbotmp::A4Box reduction_xy_DEV{};
    turbotmp::A4Box Kh_Max_xy_DEV{};
    turbotmp::A4Box Ah_Max_xy_DEV{};
    turbotmp::A4Box Ah_Max_xy_KS_DEV{};
    turbotmp::A4Box n1n2_q_DEV{};
    turbotmp::A4Box n1n1_m_n2n2_q_DEV{};
    turbotmp::A4Box dx2h_DEV{};
    turbotmp::A4Box dy2h_DEV{};
    turbotmp::A4Box dx_dyT_DEV{};
    turbotmp::A4Box dy_dxT_DEV{};
    turbotmp::A4Box m_const_leithy_DEV{};
    turbotmp::A4Box m_leithy_max_DEV{};
    turbotmp::A4Box dx2q_DEV{};
    turbotmp::A4Box dy2q_DEV{};
    turbotmp::A4Box dx_dyBu_DEV{};
    turbotmp::A4Box dy_dxBu_DEV{};
    turbotmp::A4Box Idx2dyCu_DEV{};
    turbotmp::A4Box Idxdy2u_DEV{};
    turbotmp::A4Box Idx2dyCv_DEV{};
    turbotmp::A4Box Idxdy2v_DEV{};
    turbotmp::A4Box Laplac2_const_xx_DEV{};
    turbotmp::A4Box Biharm6_const_xx_DEV{};
    turbotmp::A4Box Laplac3_const_xx_DEV{};
    turbotmp::A4Box Biharm_const_xx_DEV{};
    turbotmp::A4Box Biharm_const2_xx_DEV{};
    turbotmp::A4Box Re_Ah_const_xx_DEV{};
    turbotmp::A4Box Laplac2_const_xy_DEV{};
    turbotmp::A4Box Biharm6_const_xy_DEV{};
    turbotmp::A4Box Laplac3_const_xy_DEV{};
    turbotmp::A4Box Biharm_const_xy_DEV{};
    turbotmp::A4Box Biharm_const2_xy_DEV{};
    turbotmp::A4Box Re_Ah_const_xy_DEV{};
    if (hu_cont_allocated) hu_cont_DEV = turbotmp::make_array4(hu_cont_HOST);
    if (hv_cont_allocated) hv_cont_DEV = turbotmp::make_array4(hv_cont_HOST);
    if (MEKE_Ku_allocated) MEKE_Ku_DEV = turbotmp::make_array4(&MEKE_HOST->Ku);
    if (MEKE_Au_allocated) MEKE_Au_DEV = turbotmp::make_array4(&MEKE_HOST->Au);
    if (MEKE_mom_src_allocated) MEKE_mom_src_DEV = turbotmp::make_array4(&MEKE_HOST->mom_src);
    if (MEKE_mom_src_bh_allocated) MEKE_mom_src_bh_DEV = turbotmp::make_array4(&MEKE_HOST->mom_src_bh);
    if (MEKE_GME_snk_allocated) MEKE_GME_snk_DEV = turbotmp::make_array4(&MEKE_HOST->GME_snk);
    if (Res_fn_h_allocated) Res_fn_h_DEV = turbotmp::make_array4(&VarMix_HOST->Res_fn_h);
    if (Res_fn_q_allocated) Res_fn_q_DEV = turbotmp::make_array4(&VarMix_HOST->Res_fn_q);
    if (BS_struct_allocated) BS_struct_DEV = turbotmp::make_array4(&VarMix_HOST->BS_struct);
    if (skeb_diss_allocated) skeb_diss_DEV = turbotmp::make_array4(&STOCH_HOST->skeb_diss);
    if (Kh_bg_xx_allocated) Kh_bg_xx_DEV = turbotmp::make_array4(&metrics_HOST->Kh_bg_xx);
    if (Ah_bg_xx_allocated) Ah_bg_xx_DEV = turbotmp::make_array4(&metrics_HOST->Ah_bg_xx);
    if (reduction_xx_allocated) reduction_xx_DEV = turbotmp::make_array4(&metrics_HOST->reduction_xx);
    if (Kh_Max_xx_allocated) Kh_Max_xx_DEV = turbotmp::make_array4(&metrics_HOST->Kh_Max_xx);
    if (Ah_Max_xx_allocated) Ah_Max_xx_DEV = turbotmp::make_array4(&metrics_HOST->Ah_Max_xx);
    if (Ah_Max_xx_KS_allocated) Ah_Max_xx_KS_DEV = turbotmp::make_array4(&metrics_HOST->Ah_Max_xx_KS);
    if (n1n2_h_allocated) n1n2_h_DEV = turbotmp::make_array4(&metrics_HOST->n1n2_h);
    if (n1n1_m_n2n2_h_allocated) n1n1_m_n2n2_h_DEV = turbotmp::make_array4(&metrics_HOST->n1n1_m_n2n2_h);
    if (grid_sp_h2_allocated) grid_sp_h2_DEV = turbotmp::make_array4(&metrics_HOST->grid_sp_h2);
    if (grid_sp_h3_allocated) grid_sp_h3_DEV = turbotmp::make_array4(&metrics_HOST->grid_sp_h3);
    if (Kh_bg_xy_allocated) Kh_bg_xy_DEV = turbotmp::make_array4(&metrics_HOST->Kh_bg_xy);
    if (Ah_bg_xy_allocated) Ah_bg_xy_DEV = turbotmp::make_array4(&metrics_HOST->Ah_bg_xy);
    if (reduction_xy_allocated) reduction_xy_DEV = turbotmp::make_array4(&metrics_HOST->reduction_xy);
    if (Kh_Max_xy_allocated) Kh_Max_xy_DEV = turbotmp::make_array4(&metrics_HOST->Kh_Max_xy);
    if (Ah_Max_xy_allocated) Ah_Max_xy_DEV = turbotmp::make_array4(&metrics_HOST->Ah_Max_xy);
    if (Ah_Max_xy_KS_allocated) Ah_Max_xy_KS_DEV = turbotmp::make_array4(&metrics_HOST->Ah_Max_xy_KS);
    if (n1n2_q_allocated) n1n2_q_DEV = turbotmp::make_array4(&metrics_HOST->n1n2_q);
    if (n1n1_m_n2n2_q_allocated) n1n1_m_n2n2_q_DEV = turbotmp::make_array4(&metrics_HOST->n1n1_m_n2n2_q);
    if (dx2h_allocated) dx2h_DEV = turbotmp::make_array4(&metrics_HOST->dx2h);
    if (dy2h_allocated) dy2h_DEV = turbotmp::make_array4(&metrics_HOST->dy2h);
    if (dx_dyT_allocated) dx_dyT_DEV = turbotmp::make_array4(&metrics_HOST->dx_dyT);
    if (dy_dxT_allocated) dy_dxT_DEV = turbotmp::make_array4(&metrics_HOST->dy_dxT);
    if (m_const_leithy_allocated) m_const_leithy_DEV = turbotmp::make_array4(&metrics_HOST->m_const_leithy);
    if (m_leithy_max_allocated) m_leithy_max_DEV = turbotmp::make_array4(&metrics_HOST->m_leithy_max);
    if (dx2q_allocated) dx2q_DEV = turbotmp::make_array4(&metrics_HOST->dx2q);
    if (dy2q_allocated) dy2q_DEV = turbotmp::make_array4(&metrics_HOST->dy2q);
    if (dx_dyBu_allocated) dx_dyBu_DEV = turbotmp::make_array4(&metrics_HOST->dx_dyBu);
    if (dy_dxBu_allocated) dy_dxBu_DEV = turbotmp::make_array4(&metrics_HOST->dy_dxBu);
    if (Idx2dyCu_allocated) Idx2dyCu_DEV = turbotmp::make_array4(&metrics_HOST->Idx2dyCu);
    if (Idxdy2u_allocated) Idxdy2u_DEV = turbotmp::make_array4(&metrics_HOST->Idxdy2u);
    if (Idx2dyCv_allocated) Idx2dyCv_DEV = turbotmp::make_array4(&metrics_HOST->Idx2dyCv);
    if (Idxdy2v_allocated) Idxdy2v_DEV = turbotmp::make_array4(&metrics_HOST->Idxdy2v);
    if (Laplac2_const_xx_allocated) Laplac2_const_xx_DEV = turbotmp::make_array4(&metrics_HOST->Laplac2_const_xx);
    if (Biharm6_const_xx_allocated) Biharm6_const_xx_DEV = turbotmp::make_array4(&metrics_HOST->Biharm6_const_xx);
    if (Laplac3_const_xx_allocated) Laplac3_const_xx_DEV = turbotmp::make_array4(&metrics_HOST->Laplac3_const_xx);
    if (Biharm_const_xx_allocated) Biharm_const_xx_DEV = turbotmp::make_array4(&metrics_HOST->Biharm_const_xx);
    if (Biharm_const2_xx_allocated) Biharm_const2_xx_DEV = turbotmp::make_array4(&metrics_HOST->Biharm_const2_xx);
    if (Re_Ah_const_xx_allocated) Re_Ah_const_xx_DEV = turbotmp::make_array4(&metrics_HOST->Re_Ah_const_xx);
    if (Laplac2_const_xy_allocated) Laplac2_const_xy_DEV = turbotmp::make_array4(&metrics_HOST->Laplac2_const_xy);
    if (Biharm6_const_xy_allocated) Biharm6_const_xy_DEV = turbotmp::make_array4(&metrics_HOST->Biharm6_const_xy);
    if (Laplac3_const_xy_allocated) Laplac3_const_xy_DEV = turbotmp::make_array4(&metrics_HOST->Laplac3_const_xy);
    if (Biharm_const_xy_allocated) Biharm_const_xy_DEV = turbotmp::make_array4(&metrics_HOST->Biharm_const_xy);
    if (Biharm_const2_xy_allocated) Biharm_const2_xy_DEV = turbotmp::make_array4(&metrics_HOST->Biharm_const2_xy);
    if (Re_Ah_const_xy_allocated) Re_Ah_const_xy_DEV = turbotmp::make_array4(&metrics_HOST->Re_Ah_const_xy);

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dT.data,     mask2dT_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IareaT.data,      IareaT_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdxT.data,        IdxT_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdyT.data,        IdyT_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dCu.data,    mask2dCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdxCu.data,       IdxCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->dxCu.data,        dxCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dCv.data,    mask2dCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdyCv.data,       IdyCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->dyCv.data,        dyCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->CoriolisBu.data,  CoriolisBu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdyCu.data,       IdyCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdxCv.data,       IdxCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dBu.data,    mask2dBu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->dxBu.data,        dxBu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->dyBu.data,        dyBu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IareaBu.data,     IareaBu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->df_dx.data,       df_dx_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->df_dy.data,       df_dy_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IareaCu.data,     IareaCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IareaCv.data,     IareaCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdxBu.data,       IdxBu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IdyBu.data,       IdyBu_DEV);
    turbotmp::copy_FortranHost_to_array4(u_HOST->data,                 u_DEV);
    turbotmp::copy_FortranHost_to_array4(v_HOST->data,                 v_DEV);
    turbotmp::copy_FortranHost_to_array4(h_HOST->data,                 h_DEV);
    turbotmp::copy_FortranHost_to_array4(uh_HOST->data,                uh_DEV);
    turbotmp::copy_FortranHost_to_array4(vh_HOST->data,                vh_DEV);
    turbotmp::copy_FortranHost_to_array4(diffu_HOST->data,             diffu_DEV);
    turbotmp::copy_FortranHost_to_array4(diffv_HOST->data,             diffv_DEV);
    turbotmp::copy_FortranHost_to_array4(u_smooth_HOST->data,          u_smooth_DEV);
    turbotmp::copy_FortranHost_to_array4(v_smooth_HOST->data,          v_smooth_DEV);
    turbotmp::copy_FortranHost_to_array4(sh_xx_bt_HOST->data,          sh_xx_bt_DEV);
    turbotmp::copy_FortranHost_to_array4(sh_xy_bt_HOST->data,          sh_xy_bt_DEV);
    turbotmp::copy_FortranHost_to_array4(GME_effic_h_HOST->data,       GME_effic_h_DEV);
    turbotmp::copy_FortranHost_to_array4(GME_effic_q_HOST->data,       GME_effic_q_DEV);
    turbotmp::copy_FortranHost_to_array4(KH_u_GME_HOST->data,          KH_u_GME_DEV);
    turbotmp::copy_FortranHost_to_array4(KH_v_GME_HOST->data,          KH_v_GME_DEV);
    turbotmp::copy_FortranHost_to_array4(GME_coeff_h_HOST->data,       GME_coeff_h_DEV);
    turbotmp::copy_FortranHost_to_array4(GME_coeff_q_HOST->data,       GME_coeff_q_DEV);
    turbotmp::copy_FortranHost_to_array4(str_xx_GME_HOST->data,        str_xx_GME_DEV);
    turbotmp::copy_FortranHost_to_array4(str_xy_GME_HOST->data,        str_xy_GME_DEV);
    turbotmp::copy_FortranHost_to_array4(NoSt_HOST->data,              NoSt_DEV);
    turbotmp::copy_FortranHost_to_array4(ShSt_HOST->data,              ShSt_DEV);
    turbotmp::copy_FortranHost_to_array4(FrictWork_HOST->data,         FrictWork_DEV);
    turbotmp::copy_FortranHost_to_array4(FrictWork_bh_HOST->data,      FrictWork_bh_DEV);
    turbotmp::copy_FortranHost_to_array4(FrictWork_GME_HOST->data,     FrictWork_GME_DEV);
    turbotmp::copy_FortranHost_to_array4(Ah_h_HOST->data,              Ah_h_DEV);
    turbotmp::copy_FortranHost_to_array4(Ah_q_HOST->data,              Ah_q_DEV);
    turbotmp::copy_FortranHost_to_array4(Kh_h_HOST->data,              Kh_h_DEV);
    turbotmp::copy_FortranHost_to_array4(Kh_q_HOST->data,              Kh_q_DEV);
    turbotmp::copy_FortranHost_to_array4(grid_Re_Kh_HOST->data,        grid_Re_Kh_DEV);
    turbotmp::copy_FortranHost_to_array4(grid_Re_Ah_HOST->data,        grid_Re_Ah_DEV);
    turbotmp::copy_FortranHost_to_array4(div_xx_h_HOST->data,          div_xx_h_DEV);
    turbotmp::copy_FortranHost_to_array4(sh_xx_h_HOST->data,           sh_xx_h_DEV);
    turbotmp::copy_FortranHost_to_array4(vort_xy_q_HOST->data,         vort_xy_q_DEV);
    turbotmp::copy_FortranHost_to_array4(sh_xy_q_HOST->data,           sh_xy_q_DEV);
    turbotmp::copy_FortranHost_to_array4(visc_limit_h_HOST->data,      visc_limit_h_DEV);
    turbotmp::copy_FortranHost_to_array4(visc_limit_q_HOST->data,      visc_limit_q_DEV);
    turbotmp::copy_FortranHost_to_array4(visc_limit_h_flag_HOST->data, visc_limit_h_flag_DEV);
    turbotmp::copy_FortranHost_to_array4(visc_limit_q_flag_HOST->data, visc_limit_q_flag_DEV);
    turbotmp::copy_FortranHost_to_array4(visc_limit_h_frac_HOST->data, visc_limit_h_frac_DEV);
    turbotmp::copy_FortranHost_to_array4(visc_limit_q_frac_HOST->data, visc_limit_q_frac_DEV);
    turbotmp::copy_FortranHost_to_array4(BS_coeff_h_HOST->data,        BS_coeff_h_DEV);
    turbotmp::copy_FortranHost_to_array4(BS_coeff_q_HOST->data,        BS_coeff_q_DEV);
    if (hu_cont_allocated) turbotmp::copy_FortranHost_to_array4(hu_cont_HOST->data, hu_cont_DEV);
    if (hv_cont_allocated) turbotmp::copy_FortranHost_to_array4(hv_cont_HOST->data, hv_cont_DEV);
    if (MEKE_Ku_allocated) turbotmp::copy_FortranHost_to_array4(MEKE_HOST->Ku.data, MEKE_Ku_DEV);
    if (MEKE_Au_allocated) turbotmp::copy_FortranHost_to_array4(MEKE_HOST->Au.data, MEKE_Au_DEV);
    if (MEKE_mom_src_allocated) turbotmp::copy_FortranHost_to_array4(MEKE_HOST->mom_src.data, MEKE_mom_src_DEV);
    if (MEKE_mom_src_bh_allocated)
        turbotmp::copy_FortranHost_to_array4(MEKE_HOST->mom_src_bh.data, MEKE_mom_src_bh_DEV);
    if (MEKE_GME_snk_allocated) turbotmp::copy_FortranHost_to_array4(MEKE_HOST->GME_snk.data, MEKE_GME_snk_DEV);
    if (Res_fn_h_allocated) turbotmp::copy_FortranHost_to_array4(VarMix_HOST->Res_fn_h.data, Res_fn_h_DEV);
    if (Res_fn_q_allocated) turbotmp::copy_FortranHost_to_array4(VarMix_HOST->Res_fn_q.data, Res_fn_q_DEV);
    if (BS_struct_allocated) turbotmp::copy_FortranHost_to_array4(VarMix_HOST->BS_struct.data, BS_struct_DEV);
    if (skeb_diss_allocated) turbotmp::copy_FortranHost_to_array4(STOCH_HOST->skeb_diss.data, skeb_diss_DEV);
    if (Kh_bg_xx_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Kh_bg_xx.data, Kh_bg_xx_DEV);
    if (Ah_bg_xx_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Ah_bg_xx.data, Ah_bg_xx_DEV);
    if (reduction_xx_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->reduction_xx.data, reduction_xx_DEV);
    if (Kh_Max_xx_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Kh_Max_xx.data, Kh_Max_xx_DEV);
    if (Ah_Max_xx_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Ah_Max_xx.data, Ah_Max_xx_DEV);
    if (Ah_Max_xx_KS_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Ah_Max_xx_KS.data, Ah_Max_xx_KS_DEV);
    if (n1n2_h_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->n1n2_h.data, n1n2_h_DEV);
    if (n1n1_m_n2n2_h_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->n1n1_m_n2n2_h.data, n1n1_m_n2n2_h_DEV);
    if (grid_sp_h2_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->grid_sp_h2.data, grid_sp_h2_DEV);
    if (grid_sp_h3_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->grid_sp_h3.data, grid_sp_h3_DEV);
    if (Kh_bg_xy_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Kh_bg_xy.data, Kh_bg_xy_DEV);
    if (Ah_bg_xy_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Ah_bg_xy.data, Ah_bg_xy_DEV);
    if (reduction_xy_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->reduction_xy.data, reduction_xy_DEV);
    if (Kh_Max_xy_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Kh_Max_xy.data, Kh_Max_xy_DEV);
    if (Ah_Max_xy_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Ah_Max_xy.data, Ah_Max_xy_DEV);
    if (Ah_Max_xy_KS_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Ah_Max_xy_KS.data, Ah_Max_xy_KS_DEV);
    if (n1n2_q_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->n1n2_q.data, n1n2_q_DEV);
    if (n1n1_m_n2n2_q_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->n1n1_m_n2n2_q.data, n1n1_m_n2n2_q_DEV);
    if (dx2h_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->dx2h.data, dx2h_DEV);
    if (dy2h_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->dy2h.data, dy2h_DEV);
    if (dx_dyT_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->dx_dyT.data, dx_dyT_DEV);
    if (dy_dxT_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->dy_dxT.data, dy_dxT_DEV);
    if (m_const_leithy_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->m_const_leithy.data, m_const_leithy_DEV);
    if (m_leithy_max_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->m_leithy_max.data, m_leithy_max_DEV);
    if (dx2q_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->dx2q.data, dx2q_DEV);
    if (dy2q_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->dy2q.data, dy2q_DEV);
    if (dx_dyBu_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->dx_dyBu.data, dx_dyBu_DEV);
    if (dy_dxBu_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->dy_dxBu.data, dy_dxBu_DEV);
    if (Idx2dyCu_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Idx2dyCu.data, Idx2dyCu_DEV);
    if (Idxdy2u_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Idxdy2u.data, Idxdy2u_DEV);
    if (Idx2dyCv_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Idx2dyCv.data, Idx2dyCv_DEV);
    if (Idxdy2v_allocated) turbotmp::copy_FortranHost_to_array4(metrics_HOST->Idxdy2v.data, Idxdy2v_DEV);
    if (Laplac2_const_xx_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->Laplac2_const_xx.data, Laplac2_const_xx_DEV);
    if (Biharm6_const_xx_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->Biharm6_const_xx.data, Biharm6_const_xx_DEV);
    if (Laplac3_const_xx_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->Laplac3_const_xx.data, Laplac3_const_xx_DEV);
    if (Biharm_const_xx_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->Biharm_const_xx.data, Biharm_const_xx_DEV);
    if (Biharm_const2_xx_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->Biharm_const2_xx.data, Biharm_const2_xx_DEV);
    if (Re_Ah_const_xx_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->Re_Ah_const_xx.data, Re_Ah_const_xx_DEV);
    if (Laplac2_const_xy_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->Laplac2_const_xy.data, Laplac2_const_xy_DEV);
    if (Biharm6_const_xy_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->Biharm6_const_xy.data, Biharm6_const_xy_DEV);
    if (Laplac3_const_xy_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->Laplac3_const_xy.data, Laplac3_const_xy_DEV);
    if (Biharm_const_xy_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->Biharm_const_xy.data, Biharm_const_xy_DEV);
    if (Biharm_const2_xy_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->Biharm_const2_xy.data, Biharm_const2_xy_DEV);
    if (Re_Ah_const_xy_allocated)
        turbotmp::copy_FortranHost_to_array4(metrics_HOST->Re_Ah_const_xy.data, Re_Ah_const_xy_DEV);

    /// Device views of the grid, metric and mixing fields used by the kernel
    MOM::hv_grid_dev G;
    G.mask2dT = mask2dT_DEV.arr;
    G.IareaT = IareaT_DEV.arr;
    G.IdxT = IdxT_DEV.arr;
    G.IdyT = IdyT_DEV.arr;
    G.mask2dCu = mask2dCu_DEV.arr;
    G.IdxCu = IdxCu_DEV.arr;
    G.dxCu = dxCu_DEV.arr;
    G.mask2dCv = mask2dCv_DEV.arr;
    G.IdyCv = IdyCv_DEV.arr;
    G.dyCv = dyCv_DEV.arr;
    G.CoriolisBu = CoriolisBu_DEV.arr;
    G.IdyCu = IdyCu_DEV.arr;
    G.IdxCv = IdxCv_DEV.arr;
    G.mask2dBu = mask2dBu_DEV.arr;
    G.dxBu = dxBu_DEV.arr;
    G.dyBu = dyBu_DEV.arr;
    G.IareaBu = IareaBu_DEV.arr;
    G.df_dx = df_dx_DEV.arr;
    G.df_dy = df_dy_DEV.arr;
    G.IareaCu = IareaCu_DEV.arr;
    G.IareaCv = IareaCv_DEV.arr;
    G.IdxBu = IdxBu_DEV.arr;
    G.IdyBu = IdyBu_DEV.arr;
    MOM::hv_metrics_dev M;
    M.Kh_bg_xx = Kh_bg_xx_DEV.arr;
    M.Ah_bg_xx = Ah_bg_xx_DEV.arr;
    M.reduction_xx = reduction_xx_DEV.arr;
    M.Kh_Max_xx = Kh_Max_xx_DEV.arr;
    M.Ah_Max_xx = Ah_Max_xx_DEV.arr;
    M.Ah_Max_xx_KS = Ah_Max_xx_KS_DEV.arr;
    M.n1n2_h = n1n2_h_DEV.arr;
    M.n1n1_m_n2n2_h = n1n1_m_n2n2_h_DEV.arr;
    M.grid_sp_h2 = grid_sp_h2_DEV.arr;
    M.grid_sp_h3 = grid_sp_h3_DEV.arr;
    M.Kh_bg_xy = Kh_bg_xy_DEV.arr;
    M.Ah_bg_xy = Ah_bg_xy_DEV.arr;
    M.reduction_xy = reduction_xy_DEV.arr;
    M.Kh_Max_xy = Kh_Max_xy_DEV.arr;
    M.Ah_Max_xy = Ah_Max_xy_DEV.arr;
    M.Ah_Max_xy_KS = Ah_Max_xy_KS_DEV.arr;
    M.n1n2_q = n1n2_q_DEV.arr;
    M.n1n1_m_n2n2_q = n1n1_m_n2n2_q_DEV.arr;
    M.dx2h = dx2h_DEV.arr;
    M.dy2h = dy2h_DEV.arr;
    M.dx_dyT = dx_dyT_DEV.arr;
    M.dy_dxT = dy_dxT_DEV.arr;
    M.m_const_leithy = m_const_leithy_DEV.arr;
    M.m_leithy_max = m_leithy_max_DEV.arr;
    M.dx2q = dx2q_DEV.arr;
    M.dy2q = dy2q_DEV.arr;
    M.dx_dyBu = dx_dyBu_DEV.arr;
    M.dy_dxBu = dy_dxBu_DEV.arr;
    M.Idx2dyCu = Idx2dyCu_DEV.arr;
    M.Idxdy2u = Idxdy2u_DEV.arr;
    M.Idx2dyCv = Idx2dyCv_DEV.arr;
    M.Idxdy2v = Idxdy2v_DEV.arr;
    M.Laplac2_const_xx = Laplac2_const_xx_DEV.arr;
    M.Biharm6_const_xx = Biharm6_const_xx_DEV.arr;
    M.Laplac3_const_xx = Laplac3_const_xx_DEV.arr;
    M.Biharm_const_xx = Biharm_const_xx_DEV.arr;
    M.Biharm_const2_xx = Biharm_const2_xx_DEV.arr;
    M.Re_Ah_const_xx = Re_Ah_const_xx_DEV.arr;
    M.Laplac2_const_xy = Laplac2_const_xy_DEV.arr;
    M.Biharm6_const_xy = Biharm6_const_xy_DEV.arr;
    M.Laplac3_const_xy = Laplac3_const_xy_DEV.arr;
    M.Biharm_const_xy = Biharm_const_xy_DEV.arr;
    M.Biharm_const2_xy = Biharm_const2_xy_DEV.arr;
    M.Re_Ah_const_xy = Re_Ah_const_xy_DEV.arr;
    MOM::hv_mixing_dev X;
    X.MEKE_Ku = MEKE_Ku_DEV.arr;
    X.MEKE_Au = MEKE_Au_DEV.arr;
    X.MEKE_mom_src = MEKE_mom_src_DEV.arr;
    X.MEKE_mom_src_bh = MEKE_mom_src_bh_DEV.arr;
    X.MEKE_GME_snk = MEKE_GME_snk_DEV.arr;
    X.Res_fn_h = Res_fn_h_DEV.arr;
    X.Res_fn_q = Res_fn_q_DEV.arr;
    X.BS_struct = BS_struct_DEV.arr;
    X.skeb_diss = skeb_diss_DEV.arr;
    X.backscatter_Ro_pow  = MEKE_HOST->backscatter_Ro_pow;
    X.backscatter_Ro_c    = MEKE_HOST->backscatter_Ro_c;
    X.use_variable_mixing = VarMix_HOST->use_variable_mixing;
    X.Resoln_scaled_Kh    = VarMix_HOST->Resoln_scaled_Kh;
    X.skeb_use_frict      = STOCH_HOST->skeb_use_frict;
    X.skeb_frict_coef     = STOCH_HOST->skeb_frict_coef;

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::horizontal_viscosity_tr(bxH0, bxQ0, *CS_HOST, G, M, X, OBC, nkblock, H_subroundoff, H_to_RZ, s_to_T,
                                 u_DEV.arr, v_DEV.arr, h_DEV.arr, uh_DEV.arr, vh_DEV.arr, diffu_DEV.arr,
                                 diffv_DEV.arr, hu_cont_DEV.arr, hv_cont_DEV.arr, u_smooth_DEV.arr, v_smooth_DEV.arr,
                                 sh_xx_bt_DEV.arr, sh_xy_bt_DEV.arr, GME_effic_h_DEV.arr, GME_effic_q_DEV.arr,
                                 KH_u_GME_DEV.arr, KH_v_GME_DEV.arr, GME_coeff_h_DEV.arr, GME_coeff_q_DEV.arr,
                                 str_xx_GME_DEV.arr, str_xy_GME_DEV.arr, NoSt_DEV.arr, ShSt_DEV.arr,
                                 FrictWork_DEV.arr, FrictWork_bh_DEV.arr, FrictWork_GME_DEV.arr, Ah_h_DEV.arr,
                                 Ah_q_DEV.arr, Kh_h_DEV.arr, Kh_q_DEV.arr, grid_Re_Kh_DEV.arr, grid_Re_Ah_DEV.arr,
                                 div_xx_h_DEV.arr, sh_xx_h_DEV.arr, vort_xy_q_DEV.arr, sh_xy_q_DEV.arr,
                                 visc_limit_h_DEV.arr, visc_limit_q_DEV.arr, visc_limit_h_flag_DEV.arr,
                                 visc_limit_q_flag_DEV.arr, visc_limit_h_frac_DEV.arr, visc_limit_q_frac_DEV.arr,
                                 BS_coeff_h_DEV.arr, BS_coeff_q_DEV.arr);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(diffu_DEV,             diffu_HOST->data);
    turbotmp::copy_array4_to_FortranHost(diffv_DEV,             diffv_HOST->data);
    turbotmp::copy_array4_to_FortranHost(GME_coeff_h_DEV,       GME_coeff_h_HOST->data);
    turbotmp::copy_array4_to_FortranHost(GME_coeff_q_DEV,       GME_coeff_q_HOST->data);
    turbotmp::copy_array4_to_FortranHost(str_xx_GME_DEV,        str_xx_GME_HOST->data);
    turbotmp::copy_array4_to_FortranHost(str_xy_GME_DEV,        str_xy_GME_HOST->data);
    turbotmp::copy_array4_to_FortranHost(NoSt_DEV,              NoSt_HOST->data);
    turbotmp::copy_array4_to_FortranHost(ShSt_DEV,              ShSt_HOST->data);
    turbotmp::copy_array4_to_FortranHost(FrictWork_DEV,         FrictWork_HOST->data);
    turbotmp::copy_array4_to_FortranHost(FrictWork_bh_DEV,      FrictWork_bh_HOST->data);
    turbotmp::copy_array4_to_FortranHost(FrictWork_GME_DEV,     FrictWork_GME_HOST->data);
    turbotmp::copy_array4_to_FortranHost(Ah_h_DEV,              Ah_h_HOST->data);
    turbotmp::copy_array4_to_FortranHost(Ah_q_DEV,              Ah_q_HOST->data);
    turbotmp::copy_array4_to_FortranHost(Kh_h_DEV,              Kh_h_HOST->data);
    turbotmp::copy_array4_to_FortranHost(Kh_q_DEV,              Kh_q_HOST->data);
    turbotmp::copy_array4_to_FortranHost(grid_Re_Kh_DEV,        grid_Re_Kh_HOST->data);
    turbotmp::copy_array4_to_FortranHost(grid_Re_Ah_DEV,        grid_Re_Ah_HOST->data);
    turbotmp::copy_array4_to_FortranHost(div_xx_h_DEV,          div_xx_h_HOST->data);
    turbotmp::copy_array4_to_FortranHost(sh_xx_h_DEV,           sh_xx_h_HOST->data);
    turbotmp::copy_array4_to_FortranHost(vort_xy_q_DEV,         vort_xy_q_HOST->data);
    turbotmp::copy_array4_to_FortranHost(sh_xy_q_DEV,           sh_xy_q_HOST->data);
    turbotmp::copy_array4_to_FortranHost(visc_limit_h_DEV,      visc_limit_h_HOST->data);
    turbotmp::copy_array4_to_FortranHost(visc_limit_q_DEV,      visc_limit_q_HOST->data);
    turbotmp::copy_array4_to_FortranHost(visc_limit_h_flag_DEV, visc_limit_h_flag_HOST->data);
    turbotmp::copy_array4_to_FortranHost(visc_limit_q_flag_DEV, visc_limit_q_flag_HOST->data);
    turbotmp::copy_array4_to_FortranHost(visc_limit_h_frac_DEV, visc_limit_h_frac_HOST->data);
    turbotmp::copy_array4_to_FortranHost(visc_limit_q_frac_DEV, visc_limit_q_frac_HOST->data);
    turbotmp::copy_array4_to_FortranHost(BS_coeff_h_DEV,        BS_coeff_h_HOST->data);
    turbotmp::copy_array4_to_FortranHost(BS_coeff_q_DEV,        BS_coeff_q_HOST->data);
    if (hu_cont_allocated) turbotmp::copy_array4_to_FortranHost(hu_cont_DEV, hu_cont_HOST->data);
    if (hv_cont_allocated) turbotmp::copy_array4_to_FortranHost(hv_cont_DEV, hv_cont_HOST->data);
    if (MEKE_mom_src_allocated) turbotmp::copy_array4_to_FortranHost(MEKE_mom_src_DEV, MEKE_HOST->mom_src.data);
    if (MEKE_mom_src_bh_allocated)
        turbotmp::copy_array4_to_FortranHost(MEKE_mom_src_bh_DEV, MEKE_HOST->mom_src_bh.data);
    if (MEKE_GME_snk_allocated) turbotmp::copy_array4_to_FortranHost(MEKE_GME_snk_DEV, MEKE_HOST->GME_snk.data);
    if (skeb_diss_allocated) turbotmp::copy_array4_to_FortranHost(skeb_diss_DEV, STOCH_HOST->skeb_diss.data);

    /// Free everything we allocated (freeing a never-allocated container is a no-op)
    turbotmp::free_array4(mask2dT_DEV);
    turbotmp::free_array4(IareaT_DEV);
    turbotmp::free_array4(IdxT_DEV);
    turbotmp::free_array4(IdyT_DEV);
    turbotmp::free_array4(mask2dCu_DEV);
    turbotmp::free_array4(IdxCu_DEV);
    turbotmp::free_array4(dxCu_DEV);
    turbotmp::free_array4(mask2dCv_DEV);
    turbotmp::free_array4(IdyCv_DEV);
    turbotmp::free_array4(dyCv_DEV);
    turbotmp::free_array4(CoriolisBu_DEV);
    turbotmp::free_array4(IdyCu_DEV);
    turbotmp::free_array4(IdxCv_DEV);
    turbotmp::free_array4(mask2dBu_DEV);
    turbotmp::free_array4(dxBu_DEV);
    turbotmp::free_array4(dyBu_DEV);
    turbotmp::free_array4(IareaBu_DEV);
    turbotmp::free_array4(df_dx_DEV);
    turbotmp::free_array4(df_dy_DEV);
    turbotmp::free_array4(IareaCu_DEV);
    turbotmp::free_array4(IareaCv_DEV);
    turbotmp::free_array4(IdxBu_DEV);
    turbotmp::free_array4(IdyBu_DEV);
    turbotmp::free_array4(u_DEV);
    turbotmp::free_array4(v_DEV);
    turbotmp::free_array4(h_DEV);
    turbotmp::free_array4(uh_DEV);
    turbotmp::free_array4(vh_DEV);
    turbotmp::free_array4(diffu_DEV);
    turbotmp::free_array4(diffv_DEV);
    turbotmp::free_array4(u_smooth_DEV);
    turbotmp::free_array4(v_smooth_DEV);
    turbotmp::free_array4(sh_xx_bt_DEV);
    turbotmp::free_array4(sh_xy_bt_DEV);
    turbotmp::free_array4(GME_effic_h_DEV);
    turbotmp::free_array4(GME_effic_q_DEV);
    turbotmp::free_array4(KH_u_GME_DEV);
    turbotmp::free_array4(KH_v_GME_DEV);
    turbotmp::free_array4(GME_coeff_h_DEV);
    turbotmp::free_array4(GME_coeff_q_DEV);
    turbotmp::free_array4(str_xx_GME_DEV);
    turbotmp::free_array4(str_xy_GME_DEV);
    turbotmp::free_array4(NoSt_DEV);
    turbotmp::free_array4(ShSt_DEV);
    turbotmp::free_array4(FrictWork_DEV);
    turbotmp::free_array4(FrictWork_bh_DEV);
    turbotmp::free_array4(FrictWork_GME_DEV);
    turbotmp::free_array4(Ah_h_DEV);
    turbotmp::free_array4(Ah_q_DEV);
    turbotmp::free_array4(Kh_h_DEV);
    turbotmp::free_array4(Kh_q_DEV);
    turbotmp::free_array4(grid_Re_Kh_DEV);
    turbotmp::free_array4(grid_Re_Ah_DEV);
    turbotmp::free_array4(div_xx_h_DEV);
    turbotmp::free_array4(sh_xx_h_DEV);
    turbotmp::free_array4(vort_xy_q_DEV);
    turbotmp::free_array4(sh_xy_q_DEV);
    turbotmp::free_array4(visc_limit_h_DEV);
    turbotmp::free_array4(visc_limit_q_DEV);
    turbotmp::free_array4(visc_limit_h_flag_DEV);
    turbotmp::free_array4(visc_limit_q_flag_DEV);
    turbotmp::free_array4(visc_limit_h_frac_DEV);
    turbotmp::free_array4(visc_limit_q_frac_DEV);
    turbotmp::free_array4(BS_coeff_h_DEV);
    turbotmp::free_array4(BS_coeff_q_DEV);
    turbotmp::free_array4(hu_cont_DEV);
    turbotmp::free_array4(hv_cont_DEV);
    turbotmp::free_array4(MEKE_Ku_DEV);
    turbotmp::free_array4(MEKE_Au_DEV);
    turbotmp::free_array4(MEKE_mom_src_DEV);
    turbotmp::free_array4(MEKE_mom_src_bh_DEV);
    turbotmp::free_array4(MEKE_GME_snk_DEV);
    turbotmp::free_array4(Res_fn_h_DEV);
    turbotmp::free_array4(Res_fn_q_DEV);
    turbotmp::free_array4(BS_struct_DEV);
    turbotmp::free_array4(skeb_diss_DEV);
    turbotmp::free_array4(Kh_bg_xx_DEV);
    turbotmp::free_array4(Ah_bg_xx_DEV);
    turbotmp::free_array4(reduction_xx_DEV);
    turbotmp::free_array4(Kh_Max_xx_DEV);
    turbotmp::free_array4(Ah_Max_xx_DEV);
    turbotmp::free_array4(Ah_Max_xx_KS_DEV);
    turbotmp::free_array4(n1n2_h_DEV);
    turbotmp::free_array4(n1n1_m_n2n2_h_DEV);
    turbotmp::free_array4(grid_sp_h2_DEV);
    turbotmp::free_array4(grid_sp_h3_DEV);
    turbotmp::free_array4(Kh_bg_xy_DEV);
    turbotmp::free_array4(Ah_bg_xy_DEV);
    turbotmp::free_array4(reduction_xy_DEV);
    turbotmp::free_array4(Kh_Max_xy_DEV);
    turbotmp::free_array4(Ah_Max_xy_DEV);
    turbotmp::free_array4(Ah_Max_xy_KS_DEV);
    turbotmp::free_array4(n1n2_q_DEV);
    turbotmp::free_array4(n1n1_m_n2n2_q_DEV);
    turbotmp::free_array4(dx2h_DEV);
    turbotmp::free_array4(dy2h_DEV);
    turbotmp::free_array4(dx_dyT_DEV);
    turbotmp::free_array4(dy_dxT_DEV);
    turbotmp::free_array4(m_const_leithy_DEV);
    turbotmp::free_array4(m_leithy_max_DEV);
    turbotmp::free_array4(dx2q_DEV);
    turbotmp::free_array4(dy2q_DEV);
    turbotmp::free_array4(dx_dyBu_DEV);
    turbotmp::free_array4(dy_dxBu_DEV);
    turbotmp::free_array4(Idx2dyCu_DEV);
    turbotmp::free_array4(Idxdy2u_DEV);
    turbotmp::free_array4(Idx2dyCv_DEV);
    turbotmp::free_array4(Idxdy2v_DEV);
    turbotmp::free_array4(Laplac2_const_xx_DEV);
    turbotmp::free_array4(Biharm6_const_xx_DEV);
    turbotmp::free_array4(Laplac3_const_xx_DEV);
    turbotmp::free_array4(Biharm_const_xx_DEV);
    turbotmp::free_array4(Biharm_const2_xx_DEV);
    turbotmp::free_array4(Re_Ah_const_xx_DEV);
    turbotmp::free_array4(Laplac2_const_xy_DEV);
    turbotmp::free_array4(Biharm6_const_xy_DEV);
    turbotmp::free_array4(Laplac3_const_xy_DEV);
    turbotmp::free_array4(Biharm_const_xy_DEV);
    turbotmp::free_array4(Biharm_const2_xy_DEV);
    turbotmp::free_array4(Re_Ah_const_xy_DEV);
}
