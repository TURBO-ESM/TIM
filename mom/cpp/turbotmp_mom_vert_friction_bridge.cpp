// SKILLS: 0.3.3
/// @file turbotmp_mom_vert_friction_bridge.cpp
/// @brief Marshalling bridges between the MOM6 Fortran shims in MOM_vert_friction
///        (mode = AMREX) and the AMReX kernels in mom_vert_friction.cpp.

#include <AMReX_Arena.H>
#include <AMReX_Gpu.H>

#include "mom_vert_friction.hpp"
#include "turbotmp_helper.hpp"
#include "turbotmp_mom_vert_friction_bridge.h"

using namespace amrex;

/// @brief Bridge for the MOM6 vertvisc_limit_vel subroutine.
///
/// @param bxU0_HOST           The u-point box of the whole column, [IscB:IecB, jsc:jec, 1:nz]
/// @param bxV0_HOST           The v-point box of the whole column, [isc:iec, JscB:JecB, 1:nz]
/// @param u_HOST              Zonal velocity [L T-1 ~> m s-1]
/// @param v_HOST              Meridional velocity [L T-1 ~> m s-1]
/// @param dt                  Time increment [T ~> s]
/// @param Gcore_HOST          Persistent copies of the core ocean grid fields
/// @param vgrid_HOST          The vertical grid and unit-scaling constants
/// @param CS_HOST             Vertical viscosity control structure
/// @param coupling_state_HOST The viscous coupling coefficients and effective layer thicknesses
///                            at velocity points
/// @param report_HOST         The velocities and flags used to report the columns in which
///                            velocity truncations occur; each field may be unallocated
///                            (data == nullptr)
/// @param ntrunc_sum          The number of velocity truncations, added to the value on entry
///                            [nondim]
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
                                        int* ntrunc_sum)
{
    /// Reconstruct the iteration Boxes (Fortran 1-based -> AMReX 0-based)
    Box bxU0 = turbotmp::to_box(bxU0_HOST);
    Box bxV0 = turbotmp::to_box(bxV0_HOST);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto u_DEV      = turbotmp::make_array4(u_HOST);
    auto v_DEV      = turbotmp::make_array4(v_HOST);
    auto areaT_DEV  = turbotmp::make_array4(&Gcore_HOST->areaT);
    auto IareaT_DEV = turbotmp::make_array4(&Gcore_HOST->IareaT);
    auto dx_Cv_DEV  = turbotmp::make_array4(&Gcore_HOST->dx_Cv);
    auto dy_Cu_DEV  = turbotmp::make_array4(&Gcore_HOST->dy_Cu);
    auto h_u_DEV    = turbotmp::make_array4(&coupling_state_HOST->h_u);
    auto h_v_DEV    = turbotmp::make_array4(&coupling_state_HOST->h_v);

    /// The report fields may be unallocated (data == nullptr); only allocate/copy them when allocated.
    const bool u_old_allocated        = (report_HOST->u_old.data != nullptr);
    const bool v_old_allocated        = (report_HOST->v_old.data != nullptr);
    const bool vel_report_u_allocated = (report_HOST->vel_report_u.data != nullptr);
    const bool vel_report_v_allocated = (report_HOST->vel_report_v.data != nullptr);
    const bool dowrite_u_allocated    = (report_HOST->dowrite_u.data != nullptr);
    const bool dowrite_v_allocated    = (report_HOST->dowrite_v.data != nullptr);
    turbotmp::A4Box    u_old_DEV{};
    turbotmp::A4Box    v_old_DEV{};
    turbotmp::A4Box    vel_report_u_DEV{};
    turbotmp::A4Box    vel_report_v_DEV{};
    turbotmp::IntA4Box dowrite_u_DEV{};
    turbotmp::IntA4Box dowrite_v_DEV{};
    if (u_old_allocated)        u_old_DEV        = turbotmp::make_array4(&report_HOST->u_old);
    if (v_old_allocated)        v_old_DEV        = turbotmp::make_array4(&report_HOST->v_old);
    if (vel_report_u_allocated) vel_report_u_DEV = turbotmp::make_array4(&report_HOST->vel_report_u);
    if (vel_report_v_allocated) vel_report_v_DEV = turbotmp::make_array4(&report_HOST->vel_report_v);
    if (dowrite_u_allocated)    dowrite_u_DEV    = turbotmp::make_int_array4(&report_HOST->dowrite_u);
    if (dowrite_v_allocated)    dowrite_v_DEV    = turbotmp::make_int_array4(&report_HOST->dowrite_v);

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(u_HOST->data,                   u_DEV);
    turbotmp::copy_FortranHost_to_array4(v_HOST->data,                   v_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->areaT.data,         areaT_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IareaT.data,        IareaT_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->dx_Cv.data,         dx_Cv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->dy_Cu.data,         dy_Cu_DEV);
    turbotmp::copy_FortranHost_to_array4(coupling_state_HOST->h_u.data,  h_u_DEV);
    turbotmp::copy_FortranHost_to_array4(coupling_state_HOST->h_v.data,  h_v_DEV);
    if (u_old_allocated)        turbotmp::copy_FortranHost_to_array4(report_HOST->u_old.data,        u_old_DEV);
    if (v_old_allocated)        turbotmp::copy_FortranHost_to_array4(report_HOST->v_old.data,        v_old_DEV);
    if (vel_report_u_allocated) turbotmp::copy_FortranHost_to_array4(report_HOST->vel_report_u.data, vel_report_u_DEV);
    if (vel_report_v_allocated) turbotmp::copy_FortranHost_to_array4(report_HOST->vel_report_v.data, vel_report_v_DEV);
    if (dowrite_u_allocated)    turbotmp::copy_FortranHost_to_int_array4(report_HOST->dowrite_u.data, dowrite_u_DEV);
    if (dowrite_v_allocated)    turbotmp::copy_FortranHost_to_int_array4(report_HOST->dowrite_v.data, dowrite_v_DEV);

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::vertvisc_limit_vel(bxU0, bxV0,
                            u_DEV.arr, v_DEV.arr, dt,
                            areaT_DEV.arr, IareaT_DEV.arr, dx_Cv_DEV.arr, dy_Cu_DEV.arr,
                            vgrid_HOST->ke, vgrid_HOST->Angstrom_H, vgrid_HOST->m_s_to_L_T,
                            *CS_HOST,
                            h_u_DEV.arr, h_v_DEV.arr,
                            u_old_DEV.arr, v_old_DEV.arr,
                            vel_report_u_DEV.arr, vel_report_v_DEV.arr,
                            dowrite_u_DEV.arr, dowrite_v_DEV.arr,
                            *ntrunc_sum);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(u_DEV, u_HOST->data);
    turbotmp::copy_array4_to_FortranHost(v_DEV, v_HOST->data);
    if (u_old_allocated)        turbotmp::copy_array4_to_FortranHost(u_old_DEV,        report_HOST->u_old.data);
    if (v_old_allocated)        turbotmp::copy_array4_to_FortranHost(v_old_DEV,        report_HOST->v_old.data);
    if (vel_report_u_allocated) turbotmp::copy_array4_to_FortranHost(vel_report_u_DEV, report_HOST->vel_report_u.data);
    if (vel_report_v_allocated) turbotmp::copy_array4_to_FortranHost(vel_report_v_DEV, report_HOST->vel_report_v.data);
    if (dowrite_u_allocated)    turbotmp::copy_int_array4_to_FortranHost(dowrite_u_DEV, report_HOST->dowrite_u.data);
    if (dowrite_v_allocated)    turbotmp::copy_int_array4_to_FortranHost(dowrite_v_DEV, report_HOST->dowrite_v.data);

    /// Free everything we allocated (freeing a never-allocated report container is a no-op)
    turbotmp::free_array4(u_DEV);
    turbotmp::free_array4(v_DEV);
    turbotmp::free_array4(areaT_DEV);
    turbotmp::free_array4(IareaT_DEV);
    turbotmp::free_array4(dx_Cv_DEV);
    turbotmp::free_array4(dy_Cu_DEV);
    turbotmp::free_array4(h_u_DEV);
    turbotmp::free_array4(h_v_DEV);
    turbotmp::free_array4(u_old_DEV);
    turbotmp::free_array4(v_old_DEV);
    turbotmp::free_array4(vel_report_u_DEV);
    turbotmp::free_array4(vel_report_v_DEV);
    turbotmp::free_int_array4(dowrite_u_DEV);
    turbotmp::free_int_array4(dowrite_v_DEV);
}

/// @brief Bridge for the MOM6 vertvisc_remnant_TR subroutine.
///
/// @param bxU0_HOST           The u-point box of the whole column, [IscB:IecB, jsc:jec, 1:nz]
/// @param bxV0_HOST           The v-point box of the whole column, [isc:iec, JscB:JecB, 1:nz]
/// @param visc_HOST           Copies of the viscosities and bottom drag; Ray_u and Ray_v may be
///                            unallocated (data == nullptr)
/// @param visc_rem_u_HOST     Fraction of a time-step's worth of a barotropic acceleration that a
///                            layer experiences after viscosity is applied in the zonal direction
///                            [nondim]
/// @param visc_rem_v_HOST     Fraction of a time-step's worth of a barotropic acceleration that a
///                            layer experiences after viscosity is applied in the meridional
///                            direction [nondim]
/// @param dt                  Time increment [T ~> s]
/// @param Gcore_HOST          Persistent copies of the core ocean grid fields
/// @param ke                  The number of layers/levels in the vertical
/// @param coupling_state_HOST The viscous coupling coefficients and effective layer thicknesses
///                            at velocity points
void turbotmp_vertvisc_remnant_tr_bridge(const Box_C* bxU0_HOST,
                                         const Box_C* bxV0_HOST,
                                         const vertvisc_container_C* visc_HOST,
                                         RealArray_C* visc_rem_u_HOST,
                                         RealArray_C* visc_rem_v_HOST,
                                         const double dt,
                                         const grid_core_C* Gcore_HOST,
                                         const int ke,
                                         const coupling_state_C* coupling_state_HOST)
{
    /// Reconstruct the iteration Boxes (Fortran 1-based -> AMReX 0-based)
    Box bxU0 = turbotmp::to_box(bxU0_HOST);
    Box bxV0 = turbotmp::to_box(bxV0_HOST);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto visc_rem_u_DEV = turbotmp::make_array4(visc_rem_u_HOST);
    auto visc_rem_v_DEV = turbotmp::make_array4(visc_rem_v_HOST);
    auto mask2dCu_DEV   = turbotmp::make_array4(&Gcore_HOST->mask2dCu);
    auto mask2dCv_DEV   = turbotmp::make_array4(&Gcore_HOST->mask2dCv);
    auto a_u_DEV        = turbotmp::make_array4(&coupling_state_HOST->a_u);
    auto a_v_DEV        = turbotmp::make_array4(&coupling_state_HOST->a_v);
    auto h_u_DEV        = turbotmp::make_array4(&coupling_state_HOST->h_u);
    auto h_v_DEV        = turbotmp::make_array4(&coupling_state_HOST->h_v);

    /// Ray_u / Ray_v may be unallocated (data == nullptr); only allocate/copy them when allocated.
    const bool Ray_u_allocated = (visc_HOST->Ray_u.data != nullptr);
    const bool Ray_v_allocated = (visc_HOST->Ray_v.data != nullptr);
    turbotmp::A4Box Ray_u_DEV{};
    turbotmp::A4Box Ray_v_DEV{};
    if (Ray_u_allocated) Ray_u_DEV = turbotmp::make_array4(&visc_HOST->Ray_u);
    if (Ray_v_allocated) Ray_v_DEV = turbotmp::make_array4(&visc_HOST->Ray_v);

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(visc_rem_u_HOST->data,          visc_rem_u_DEV);
    turbotmp::copy_FortranHost_to_array4(visc_rem_v_HOST->data,          visc_rem_v_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dCu.data,      mask2dCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dCv.data,      mask2dCv_DEV);
    turbotmp::copy_FortranHost_to_array4(coupling_state_HOST->a_u.data,  a_u_DEV);
    turbotmp::copy_FortranHost_to_array4(coupling_state_HOST->a_v.data,  a_v_DEV);
    turbotmp::copy_FortranHost_to_array4(coupling_state_HOST->h_u.data,  h_u_DEV);
    turbotmp::copy_FortranHost_to_array4(coupling_state_HOST->h_v.data,  h_v_DEV);
    if (Ray_u_allocated) turbotmp::copy_FortranHost_to_array4(visc_HOST->Ray_u.data, Ray_u_DEV);
    if (Ray_v_allocated) turbotmp::copy_FortranHost_to_array4(visc_HOST->Ray_v.data, Ray_v_DEV);

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::vertvisc_remnant_tr(bxU0, bxV0,
                             Ray_u_DEV.arr, Ray_v_DEV.arr,
                             visc_rem_u_DEV.arr, visc_rem_v_DEV.arr,
                             dt,
                             mask2dCu_DEV.arr, mask2dCv_DEV.arr,
                             ke,
                             a_u_DEV.arr, a_v_DEV.arr,
                             h_u_DEV.arr, h_v_DEV.arr);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(visc_rem_u_DEV, visc_rem_u_HOST->data);
    turbotmp::copy_array4_to_FortranHost(visc_rem_v_DEV, visc_rem_v_HOST->data);

    /// Free everything we allocated (freeing a never-allocated Ray container is a no-op)
    turbotmp::free_array4(visc_rem_u_DEV);
    turbotmp::free_array4(visc_rem_v_DEV);
    turbotmp::free_array4(mask2dCu_DEV);
    turbotmp::free_array4(mask2dCv_DEV);
    turbotmp::free_array4(a_u_DEV);
    turbotmp::free_array4(a_v_DEV);
    turbotmp::free_array4(h_u_DEV);
    turbotmp::free_array4(h_v_DEV);
    turbotmp::free_array4(Ray_u_DEV);
    turbotmp::free_array4(Ray_v_DEV);
}

/// @brief Bridge for the MOM6 vertFPmix_TR subroutine.
///
/// @param bxU0_HOST           The u-point box of the whole column, [IscB:IecB, jsc:jec, 1:nz]
/// @param bxV0_HOST           The v-point box of the whole column, [isc:iec, JscB:JecB, 1:nz]
/// @param ui_HOST             Zonal velocity after vertvisc [L T-1 ~> m s-1]
/// @param vi_HOST             Meridional velocity after vertvisc [L T-1 ~> m s-1]
/// @param hbl_h_HOST          boundary layer depth [H ~> m]
/// @param forces_HOST         Copies of the driving mechanical forces
/// @param Us_x_HOST           3d zonal Stokes drift profile [L T-1 ~> m s-1]
/// @param Us_y_HOST           3d meridional Stokes drift profile [L T-1 ~> m s-1]
/// @param dt                  Time increment [T ~> s]
/// @param Cemp_NL             empirical coefficient of non-local momentum mixing [nondim]
/// @param Irho0               Inverse of the mean density rescaled to [Z L-1 R-1 ~> m3 kg-1]
/// @param Gcore_HOST          Persistent copies of the core ocean grid fields
/// @param coupling_state_HOST The viscous coupling coefficients and effective layer thicknesses
///                            at velocity points
/// @param FPdiag_HOST         The Eulerian velocities and the nonlocal stress increments for the
///                            diagnostics
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
                                  FPmix_diag_C* FPdiag_HOST)
{
    /// Reconstruct the iteration Boxes (Fortran 1-based -> AMReX 0-based)
    Box bxU0 = turbotmp::to_box(bxU0_HOST);
    Box bxV0 = turbotmp::to_box(bxV0_HOST);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto ui_DEV       = turbotmp::make_array4(ui_HOST);
    auto vi_DEV       = turbotmp::make_array4(vi_HOST);
    auto hbl_h_DEV    = turbotmp::make_array4(hbl_h_HOST);
    auto taux_DEV     = turbotmp::make_array4(&forces_HOST->taux);
    auto tauy_DEV     = turbotmp::make_array4(&forces_HOST->tauy);
    auto Us_x_DEV     = turbotmp::make_array4(Us_x_HOST);
    auto Us_y_DEV     = turbotmp::make_array4(Us_y_HOST);
    auto mask2dCu_DEV = turbotmp::make_array4(&Gcore_HOST->mask2dCu);
    auto mask2dCv_DEV = turbotmp::make_array4(&Gcore_HOST->mask2dCv);
    auto mask2dT_DEV  = turbotmp::make_array4(&Gcore_HOST->mask2dT);
    auto h_u_DEV      = turbotmp::make_array4(&coupling_state_HOST->h_u);
    auto h_v_DEV      = turbotmp::make_array4(&coupling_state_HOST->h_v);
    auto uE_u_DEV     = turbotmp::make_array4(&FPdiag_HOST->uE_u);
    auto vE_v_DEV     = turbotmp::make_array4(&FPdiag_HOST->vE_v);
    auto uInc_u_DEV   = turbotmp::make_array4(&FPdiag_HOST->uInc_u);
    auto vInc_v_DEV   = turbotmp::make_array4(&FPdiag_HOST->vInc_v);

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(ui_HOST->data,                 ui_DEV);
    turbotmp::copy_FortranHost_to_array4(vi_HOST->data,                 vi_DEV);
    turbotmp::copy_FortranHost_to_array4(hbl_h_HOST->data,              hbl_h_DEV);
    turbotmp::copy_FortranHost_to_array4(forces_HOST->taux.data,        taux_DEV);
    turbotmp::copy_FortranHost_to_array4(forces_HOST->tauy.data,        tauy_DEV);
    turbotmp::copy_FortranHost_to_array4(Us_x_HOST->data,               Us_x_DEV);
    turbotmp::copy_FortranHost_to_array4(Us_y_HOST->data,               Us_y_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dCu.data,     mask2dCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dCv.data,     mask2dCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dT.data,      mask2dT_DEV);
    turbotmp::copy_FortranHost_to_array4(coupling_state_HOST->h_u.data, h_u_DEV);
    turbotmp::copy_FortranHost_to_array4(coupling_state_HOST->h_v.data, h_v_DEV);
    turbotmp::copy_FortranHost_to_array4(FPdiag_HOST->uE_u.data,        uE_u_DEV);
    turbotmp::copy_FortranHost_to_array4(FPdiag_HOST->vE_v.data,        vE_v_DEV);
    turbotmp::copy_FortranHost_to_array4(FPdiag_HOST->uInc_u.data,      uInc_u_DEV);
    turbotmp::copy_FortranHost_to_array4(FPdiag_HOST->vInc_v.data,      vInc_v_DEV);

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::vertfpmix_tr(bxU0, bxV0,
                      ui_DEV.arr, vi_DEV.arr,
                      hbl_h_DEV.arr,
                      taux_DEV.arr, tauy_DEV.arr,
                      Us_x_DEV.arr, Us_y_DEV.arr,
                      dt, Cemp_NL, Irho0,
                      mask2dCu_DEV.arr, mask2dCv_DEV.arr, mask2dT_DEV.arr,
                      h_u_DEV.arr, h_v_DEV.arr,
                      uE_u_DEV.arr, vE_v_DEV.arr,
                      uInc_u_DEV.arr, vInc_v_DEV.arr);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(ui_DEV,     ui_HOST->data);
    turbotmp::copy_array4_to_FortranHost(vi_DEV,     vi_HOST->data);
    turbotmp::copy_array4_to_FortranHost(uE_u_DEV,   FPdiag_HOST->uE_u.data);
    turbotmp::copy_array4_to_FortranHost(vE_v_DEV,   FPdiag_HOST->vE_v.data);
    turbotmp::copy_array4_to_FortranHost(uInc_u_DEV, FPdiag_HOST->uInc_u.data);
    turbotmp::copy_array4_to_FortranHost(vInc_v_DEV, FPdiag_HOST->vInc_v.data);

    /// Free everything we allocated
    turbotmp::free_array4(ui_DEV);
    turbotmp::free_array4(vi_DEV);
    turbotmp::free_array4(hbl_h_DEV);
    turbotmp::free_array4(taux_DEV);
    turbotmp::free_array4(tauy_DEV);
    turbotmp::free_array4(Us_x_DEV);
    turbotmp::free_array4(Us_y_DEV);
    turbotmp::free_array4(mask2dCu_DEV);
    turbotmp::free_array4(mask2dCv_DEV);
    turbotmp::free_array4(mask2dT_DEV);
    turbotmp::free_array4(h_u_DEV);
    turbotmp::free_array4(h_v_DEV);
    turbotmp::free_array4(uE_u_DEV);
    turbotmp::free_array4(vE_v_DEV);
    turbotmp::free_array4(uInc_u_DEV);
    turbotmp::free_array4(vInc_v_DEV);
}

/// @brief Bridge for the MOM6 vertvisc_coef_TR subroutine.
///
/// @param bxU0_HOST           The u-point box of the whole column, [IscB:IecB, jsc:jec, 1:nz]
/// @param bxV0_HOST           The v-point box of the whole column, [isc:iec, JscB:JecB, 1:nz]
/// @param u_HOST              Zonal velocity [L T-1 ~> m s-1]
/// @param v_HOST              Meridional velocity [L T-1 ~> m s-1]
/// @param h_HOST              Layer thickness [H ~> m or kg m-2]
/// @param dz_HOST             Vertical distance across layers [Z ~> m]
/// @param forces_HOST         Copies of the driving mechanical forces; frac_shelf_u and
///                            frac_shelf_v may be unallocated (data == nullptr)
/// @param Ustar_2d_HOST       The wind friction velocity, calculated using the Boussinesq
///                            reference density or the time-evolving surface density in
///                            non-Boussinesq mode [Z T-1 ~> m s-1]; may be unallocated
/// @param visc_HOST           Copies of the viscosities and bottom drag; each field may be
///                            unallocated (data == nullptr)
/// @param SpV_avg_HOST        The layer averaged in situ specific volume [R-1 ~> m3 kg-1]; may be
///                            unallocated
/// @param dt                  Time increment [T ~> s]
/// @param Gcore_HOST          Persistent copies of the core ocean grid fields
/// @param vgrid_HOST          The vertical grid and unit-scaling constants; g_prime_a may be
///                            unallocated
/// @param CS_HOST             Vertical viscosity control structure; GL90_CS.kappa_gl90_2d may be
///                            unallocated
/// @param coupling_state_HOST The viscous coupling coefficients and effective layer thicknesses
///                            at velocity points; a_u_gl90, a_v_gl90, a1_shelf_u and a1_shelf_v
///                            may be unallocated
/// @param OBC                 Open boundary condition structure
/// @param kdgl90_struct_HOST  Vertical structure function used in GL90 diffusivity [nondim]; may
///                            be unallocated
/// @param coef_diag_HOST      The diagnostic arrays filled here, each allocated only if it is
///                            needed (data == nullptr otherwise)
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
                                      vertvisc_coef_diag_C* coef_diag_HOST)
{
    /// Reconstruct the iteration Boxes (Fortran 1-based -> AMReX 0-based)
    Box bxU0 = turbotmp::to_box(bxU0_HOST);
    Box bxV0 = turbotmp::to_box(bxV0_HOST);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto u_DEV                 = turbotmp::make_array4(u_HOST);
    auto v_DEV                 = turbotmp::make_array4(v_HOST);
    auto h_DEV                 = turbotmp::make_array4(h_HOST);
    auto dz_DEV                = turbotmp::make_array4(dz_HOST);
    auto mask2dCu_DEV          = turbotmp::make_array4(&Gcore_HOST->mask2dCu);
    auto mask2dCv_DEV          = turbotmp::make_array4(&Gcore_HOST->mask2dCv);
    auto bathyT_DEV            = turbotmp::make_array4(&Gcore_HOST->bathyT);
    auto CoriolisBu_DEV        = turbotmp::make_array4(&Gcore_HOST->CoriolisBu);
    auto a_u_DEV               = turbotmp::make_array4(&coupling_state_HOST->a_u);
    auto h_u_DEV               = turbotmp::make_array4(&coupling_state_HOST->h_u);
    auto a_v_DEV               = turbotmp::make_array4(&coupling_state_HOST->a_v);
    auto h_v_DEV               = turbotmp::make_array4(&coupling_state_HOST->h_v);

    /// Fields that may be unallocated (data == nullptr); only allocate/copy them when allocated.
    const bool frac_shelf_u_allocated      = (forces_HOST->frac_shelf_u.data != nullptr);
    const bool frac_shelf_v_allocated      = (forces_HOST->frac_shelf_v.data != nullptr);
    const bool Ustar_2d_allocated          = (Ustar_2d_HOST->data != nullptr);
    const bool Kv_bbl_u_allocated          = (visc_HOST->kv_bbl_u.data != nullptr);
    const bool Kv_bbl_v_allocated          = (visc_HOST->kv_bbl_v.data != nullptr);
    const bool bbl_thick_u_allocated       = (visc_HOST->bbl_thick_u.data != nullptr);
    const bool bbl_thick_v_allocated       = (visc_HOST->bbl_thick_v.data != nullptr);
    const bool tbl_thick_shelf_u_allocated = (visc_HOST->tbl_thick_shelf_u.data != nullptr);
    const bool tbl_thick_shelf_v_allocated = (visc_HOST->tbl_thick_shelf_v.data != nullptr);
    const bool Kv_tbl_shelf_u_allocated    = (visc_HOST->kv_tbl_shelf_u.data != nullptr);
    const bool Kv_tbl_shelf_v_allocated    = (visc_HOST->kv_tbl_shelf_v.data != nullptr);
    const bool nkml_visc_u_allocated       = (visc_HOST->nkml_visc_u.data != nullptr);
    const bool nkml_visc_v_allocated       = (visc_HOST->nkml_visc_v.data != nullptr);
    const bool Kv_shear_allocated          = (visc_HOST->Kv_shear.data != nullptr);
    const bool Kv_shear_Bu_allocated       = (visc_HOST->Kv_shear_Bu.data != nullptr);
    const bool SpV_avg_allocated           = (SpV_avg_HOST->data != nullptr);
    const bool kappa_gl90_2d_allocated     = (CS_HOST->GL90_CS.kappa_gl90_2d.data != nullptr);
    const bool kdgl90_struct_allocated     = (kdgl90_struct_HOST->data != nullptr);
    const bool a_u_gl90_allocated          = (coupling_state_HOST->a_u_gl90.data != nullptr);
    const bool a_v_gl90_allocated          = (coupling_state_HOST->a_v_gl90.data != nullptr);
    const bool a1_shelf_u_allocated        = (coupling_state_HOST->a1_shelf_u.data != nullptr);
    const bool a1_shelf_v_allocated        = (coupling_state_HOST->a1_shelf_v.data != nullptr);
    const bool hML_u_allocated             = (coef_diag_HOST->hML_u.data != nullptr);
    const bool hML_v_allocated             = (coef_diag_HOST->hML_v.data != nullptr);
    const bool Kv_u_allocated              = (coef_diag_HOST->Kv_u.data != nullptr);
    const bool Kv_v_allocated              = (coef_diag_HOST->Kv_v.data != nullptr);
    const bool Kv_gl90_u_allocated         = (coef_diag_HOST->Kv_gl90_u.data != nullptr);
    const bool Kv_gl90_v_allocated         = (coef_diag_HOST->Kv_gl90_v.data != nullptr);
    turbotmp::A4Box frac_shelf_u_DEV{};
    turbotmp::A4Box frac_shelf_v_DEV{};
    turbotmp::A4Box Ustar_2d_DEV{};
    turbotmp::A4Box Kv_bbl_u_DEV{};
    turbotmp::A4Box Kv_bbl_v_DEV{};
    turbotmp::A4Box bbl_thick_u_DEV{};
    turbotmp::A4Box bbl_thick_v_DEV{};
    turbotmp::A4Box tbl_thick_shelf_u_DEV{};
    turbotmp::A4Box tbl_thick_shelf_v_DEV{};
    turbotmp::A4Box Kv_tbl_shelf_u_DEV{};
    turbotmp::A4Box Kv_tbl_shelf_v_DEV{};
    turbotmp::A4Box nkml_visc_u_DEV{};
    turbotmp::A4Box nkml_visc_v_DEV{};
    turbotmp::A4Box Kv_shear_DEV{};
    turbotmp::A4Box Kv_shear_Bu_DEV{};
    turbotmp::A4Box SpV_avg_DEV{};
    turbotmp::A4Box kappa_gl90_2d_DEV{};
    turbotmp::A4Box kdgl90_struct_DEV{};
    turbotmp::A4Box a_u_gl90_DEV{};
    turbotmp::A4Box a_v_gl90_DEV{};
    turbotmp::A4Box a1_shelf_u_DEV{};
    turbotmp::A4Box a1_shelf_v_DEV{};
    turbotmp::A4Box hML_u_DEV{};
    turbotmp::A4Box hML_v_DEV{};
    turbotmp::A4Box Kv_u_DEV{};
    turbotmp::A4Box Kv_v_DEV{};
    turbotmp::A4Box Kv_gl90_u_DEV{};
    turbotmp::A4Box Kv_gl90_v_DEV{};
    if (frac_shelf_u_allocated) frac_shelf_u_DEV      = turbotmp::make_array4(&forces_HOST->frac_shelf_u);
    if (frac_shelf_v_allocated) frac_shelf_v_DEV      = turbotmp::make_array4(&forces_HOST->frac_shelf_v);
    if (Ustar_2d_allocated) Ustar_2d_DEV          = turbotmp::make_array4(Ustar_2d_HOST);
    if (Kv_bbl_u_allocated) Kv_bbl_u_DEV          = turbotmp::make_array4(&visc_HOST->kv_bbl_u);
    if (Kv_bbl_v_allocated) Kv_bbl_v_DEV          = turbotmp::make_array4(&visc_HOST->kv_bbl_v);
    if (bbl_thick_u_allocated) bbl_thick_u_DEV       = turbotmp::make_array4(&visc_HOST->bbl_thick_u);
    if (bbl_thick_v_allocated) bbl_thick_v_DEV       = turbotmp::make_array4(&visc_HOST->bbl_thick_v);
    if (tbl_thick_shelf_u_allocated) tbl_thick_shelf_u_DEV = turbotmp::make_array4(&visc_HOST->tbl_thick_shelf_u);
    if (tbl_thick_shelf_v_allocated) tbl_thick_shelf_v_DEV = turbotmp::make_array4(&visc_HOST->tbl_thick_shelf_v);
    if (Kv_tbl_shelf_u_allocated) Kv_tbl_shelf_u_DEV    = turbotmp::make_array4(&visc_HOST->kv_tbl_shelf_u);
    if (Kv_tbl_shelf_v_allocated) Kv_tbl_shelf_v_DEV    = turbotmp::make_array4(&visc_HOST->kv_tbl_shelf_v);
    if (nkml_visc_u_allocated) nkml_visc_u_DEV       = turbotmp::make_array4(&visc_HOST->nkml_visc_u);
    if (nkml_visc_v_allocated) nkml_visc_v_DEV       = turbotmp::make_array4(&visc_HOST->nkml_visc_v);
    if (Kv_shear_allocated) Kv_shear_DEV          = turbotmp::make_array4(&visc_HOST->Kv_shear);
    if (Kv_shear_Bu_allocated) Kv_shear_Bu_DEV       = turbotmp::make_array4(&visc_HOST->Kv_shear_Bu);
    if (SpV_avg_allocated) SpV_avg_DEV           = turbotmp::make_array4(SpV_avg_HOST);
    if (kappa_gl90_2d_allocated) kappa_gl90_2d_DEV     = turbotmp::make_array4(&CS_HOST->GL90_CS.kappa_gl90_2d);
    if (kdgl90_struct_allocated) kdgl90_struct_DEV     = turbotmp::make_array4(kdgl90_struct_HOST);
    if (a_u_gl90_allocated) a_u_gl90_DEV          = turbotmp::make_array4(&coupling_state_HOST->a_u_gl90);
    if (a_v_gl90_allocated) a_v_gl90_DEV          = turbotmp::make_array4(&coupling_state_HOST->a_v_gl90);
    if (a1_shelf_u_allocated) a1_shelf_u_DEV        = turbotmp::make_array4(&coupling_state_HOST->a1_shelf_u);
    if (a1_shelf_v_allocated) a1_shelf_v_DEV        = turbotmp::make_array4(&coupling_state_HOST->a1_shelf_v);
    if (hML_u_allocated) hML_u_DEV             = turbotmp::make_array4(&coef_diag_HOST->hML_u);
    if (hML_v_allocated) hML_v_DEV             = turbotmp::make_array4(&coef_diag_HOST->hML_v);
    if (Kv_u_allocated) Kv_u_DEV              = turbotmp::make_array4(&coef_diag_HOST->Kv_u);
    if (Kv_v_allocated) Kv_v_DEV              = turbotmp::make_array4(&coef_diag_HOST->Kv_v);
    if (Kv_gl90_u_allocated) Kv_gl90_u_DEV         = turbotmp::make_array4(&coef_diag_HOST->Kv_gl90_u);
    if (Kv_gl90_v_allocated) Kv_gl90_v_DEV         = turbotmp::make_array4(&coef_diag_HOST->Kv_gl90_v);

    const bool g_prime_allocated = (vgrid_HOST->g_prime_a.data != nullptr);
    Real* g_prime_data = nullptr;
    Array4<const Real> g_prime_arr{};
    if (g_prime_allocated) {
        const int n_g  = vgrid_HOST->g_prime_a.shape[0];
        const int lb_g = vgrid_HOST->g_prime_a.lb[0];
        g_prime_data = static_cast<Real*>(The_Arena()->alloc(n_g * sizeof(Real)));
        Gpu::htod_memcpy(g_prime_data, vgrid_HOST->g_prime_a.data, n_g * sizeof(Real));
        g_prime_arr = Array4<const Real>(g_prime_data, Dim3{lb_g-1, 0, 0}, Dim3{lb_g-1+n_g, 1, 1}, 1);
    }

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(u_HOST->data,                                       u_DEV);
    turbotmp::copy_FortranHost_to_array4(v_HOST->data,                                       v_DEV);
    turbotmp::copy_FortranHost_to_array4(h_HOST->data,                                       h_DEV);
    turbotmp::copy_FortranHost_to_array4(dz_HOST->data,                                      dz_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dCu.data,                          mask2dCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dCv.data,                          mask2dCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->bathyT.data,                            bathyT_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->CoriolisBu.data,                        CoriolisBu_DEV);
    turbotmp::copy_FortranHost_to_array4(coupling_state_HOST->a_u.data,                      a_u_DEV);
    turbotmp::copy_FortranHost_to_array4(coupling_state_HOST->h_u.data,                      h_u_DEV);
    turbotmp::copy_FortranHost_to_array4(coupling_state_HOST->a_v.data,                      a_v_DEV);
    turbotmp::copy_FortranHost_to_array4(coupling_state_HOST->h_v.data,                      h_v_DEV);
    if (frac_shelf_u_allocated) turbotmp::copy_FortranHost_to_array4(forces_HOST->frac_shelf_u.data, frac_shelf_u_DEV);
    if (frac_shelf_v_allocated) turbotmp::copy_FortranHost_to_array4(forces_HOST->frac_shelf_v.data, frac_shelf_v_DEV);
    if (Ustar_2d_allocated) turbotmp::copy_FortranHost_to_array4(Ustar_2d_HOST->data, Ustar_2d_DEV);
    if (Kv_bbl_u_allocated) turbotmp::copy_FortranHost_to_array4(visc_HOST->kv_bbl_u.data, Kv_bbl_u_DEV);
    if (Kv_bbl_v_allocated) turbotmp::copy_FortranHost_to_array4(visc_HOST->kv_bbl_v.data, Kv_bbl_v_DEV);
    if (bbl_thick_u_allocated) turbotmp::copy_FortranHost_to_array4(visc_HOST->bbl_thick_u.data, bbl_thick_u_DEV);
    if (bbl_thick_v_allocated) turbotmp::copy_FortranHost_to_array4(visc_HOST->bbl_thick_v.data, bbl_thick_v_DEV);
    if (tbl_thick_shelf_u_allocated) turbotmp::copy_FortranHost_to_array4(visc_HOST->tbl_thick_shelf_u.data, tbl_thick_shelf_u_DEV);
    if (tbl_thick_shelf_v_allocated) turbotmp::copy_FortranHost_to_array4(visc_HOST->tbl_thick_shelf_v.data, tbl_thick_shelf_v_DEV);
    if (Kv_tbl_shelf_u_allocated) turbotmp::copy_FortranHost_to_array4(visc_HOST->kv_tbl_shelf_u.data, Kv_tbl_shelf_u_DEV);
    if (Kv_tbl_shelf_v_allocated) turbotmp::copy_FortranHost_to_array4(visc_HOST->kv_tbl_shelf_v.data, Kv_tbl_shelf_v_DEV);
    if (nkml_visc_u_allocated) turbotmp::copy_FortranHost_to_array4(visc_HOST->nkml_visc_u.data, nkml_visc_u_DEV);
    if (nkml_visc_v_allocated) turbotmp::copy_FortranHost_to_array4(visc_HOST->nkml_visc_v.data, nkml_visc_v_DEV);
    if (Kv_shear_allocated) turbotmp::copy_FortranHost_to_array4(visc_HOST->Kv_shear.data, Kv_shear_DEV);
    if (Kv_shear_Bu_allocated) turbotmp::copy_FortranHost_to_array4(visc_HOST->Kv_shear_Bu.data, Kv_shear_Bu_DEV);
    if (SpV_avg_allocated) turbotmp::copy_FortranHost_to_array4(SpV_avg_HOST->data, SpV_avg_DEV);
    if (kappa_gl90_2d_allocated) turbotmp::copy_FortranHost_to_array4(CS_HOST->GL90_CS.kappa_gl90_2d.data, kappa_gl90_2d_DEV);
    if (kdgl90_struct_allocated) turbotmp::copy_FortranHost_to_array4(kdgl90_struct_HOST->data, kdgl90_struct_DEV);
    if (a_u_gl90_allocated) turbotmp::copy_FortranHost_to_array4(coupling_state_HOST->a_u_gl90.data, a_u_gl90_DEV);
    if (a_v_gl90_allocated) turbotmp::copy_FortranHost_to_array4(coupling_state_HOST->a_v_gl90.data, a_v_gl90_DEV);
    if (a1_shelf_u_allocated) turbotmp::copy_FortranHost_to_array4(coupling_state_HOST->a1_shelf_u.data, a1_shelf_u_DEV);
    if (a1_shelf_v_allocated) turbotmp::copy_FortranHost_to_array4(coupling_state_HOST->a1_shelf_v.data, a1_shelf_v_DEV);
    if (hML_u_allocated) turbotmp::copy_FortranHost_to_array4(coef_diag_HOST->hML_u.data, hML_u_DEV);
    if (hML_v_allocated) turbotmp::copy_FortranHost_to_array4(coef_diag_HOST->hML_v.data, hML_v_DEV);
    if (Kv_u_allocated) turbotmp::copy_FortranHost_to_array4(coef_diag_HOST->Kv_u.data, Kv_u_DEV);
    if (Kv_v_allocated) turbotmp::copy_FortranHost_to_array4(coef_diag_HOST->Kv_v.data, Kv_v_DEV);
    if (Kv_gl90_u_allocated) turbotmp::copy_FortranHost_to_array4(coef_diag_HOST->Kv_gl90_u.data, Kv_gl90_u_DEV);
    if (Kv_gl90_v_allocated) turbotmp::copy_FortranHost_to_array4(coef_diag_HOST->Kv_gl90_v.data, Kv_gl90_v_DEV);

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::vertvisc_coef_tr(bxU0, bxV0,
                          u_DEV.arr, v_DEV.arr, h_DEV.arr, dz_DEV.arr,
                          frac_shelf_u_DEV.arr, frac_shelf_v_DEV.arr,
                          Ustar_2d_DEV.arr,
                          Kv_bbl_u_DEV.arr, Kv_bbl_v_DEV.arr,
                          bbl_thick_u_DEV.arr, bbl_thick_v_DEV.arr,
                          tbl_thick_shelf_u_DEV.arr, tbl_thick_shelf_v_DEV.arr,
                          Kv_tbl_shelf_u_DEV.arr, Kv_tbl_shelf_v_DEV.arr,
                          nkml_visc_u_DEV.arr, nkml_visc_v_DEV.arr,
                          Kv_shear_DEV.arr, Kv_shear_Bu_DEV.arr,
                          SpV_avg_DEV.arr,
                          dt,
                          mask2dCu_DEV.arr, mask2dCv_DEV.arr, bathyT_DEV.arr, CoriolisBu_DEV.arr,
                          vgrid_HOST->ke, vgrid_HOST->Boussinesq, vgrid_HOST->nkml,
                          vgrid_HOST->dZ_subroundoff, vgrid_HOST->H_subroundoff,
                          vgrid_HOST->H_to_m, vgrid_HOST->m_to_H, vgrid_HOST->RZ_to_H,
                          vgrid_HOST->Z_to_H, vgrid_HOST->T_to_s,
                          g_prime_arr,
                          *CS_HOST,
                          kappa_gl90_2d_DEV.arr,
                          a_u_DEV.arr, a_u_gl90_DEV.arr, h_u_DEV.arr,
                          a_v_DEV.arr, a_v_gl90_DEV.arr, h_v_DEV.arr,
                          a1_shelf_u_DEV.arr, a1_shelf_v_DEV.arr,
                          OBC,
                          kdgl90_struct_DEV.arr,
                          hML_u_DEV.arr, hML_v_DEV.arr,
                          Kv_u_DEV.arr, Kv_v_DEV.arr,
                          Kv_gl90_u_DEV.arr, Kv_gl90_v_DEV.arr);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(a_u_DEV, coupling_state_HOST->a_u.data);
    turbotmp::copy_array4_to_FortranHost(h_u_DEV, coupling_state_HOST->h_u.data);
    turbotmp::copy_array4_to_FortranHost(a_v_DEV, coupling_state_HOST->a_v.data);
    turbotmp::copy_array4_to_FortranHost(h_v_DEV, coupling_state_HOST->h_v.data);
    if (a_u_gl90_allocated) turbotmp::copy_array4_to_FortranHost(a_u_gl90_DEV, coupling_state_HOST->a_u_gl90.data);
    if (a_v_gl90_allocated) turbotmp::copy_array4_to_FortranHost(a_v_gl90_DEV, coupling_state_HOST->a_v_gl90.data);
    if (a1_shelf_u_allocated) turbotmp::copy_array4_to_FortranHost(a1_shelf_u_DEV, coupling_state_HOST->a1_shelf_u.data);
    if (a1_shelf_v_allocated) turbotmp::copy_array4_to_FortranHost(a1_shelf_v_DEV, coupling_state_HOST->a1_shelf_v.data);
    if (hML_u_allocated) turbotmp::copy_array4_to_FortranHost(hML_u_DEV, coef_diag_HOST->hML_u.data);
    if (hML_v_allocated) turbotmp::copy_array4_to_FortranHost(hML_v_DEV, coef_diag_HOST->hML_v.data);
    if (Kv_u_allocated) turbotmp::copy_array4_to_FortranHost(Kv_u_DEV, coef_diag_HOST->Kv_u.data);
    if (Kv_v_allocated) turbotmp::copy_array4_to_FortranHost(Kv_v_DEV, coef_diag_HOST->Kv_v.data);
    if (Kv_gl90_u_allocated) turbotmp::copy_array4_to_FortranHost(Kv_gl90_u_DEV, coef_diag_HOST->Kv_gl90_u.data);
    if (Kv_gl90_v_allocated) turbotmp::copy_array4_to_FortranHost(Kv_gl90_v_DEV, coef_diag_HOST->Kv_gl90_v.data);

    /// Free everything we allocated (freeing a never-allocated container is a no-op)
    turbotmp::free_array4(u_DEV);
    turbotmp::free_array4(v_DEV);
    turbotmp::free_array4(h_DEV);
    turbotmp::free_array4(dz_DEV);
    turbotmp::free_array4(mask2dCu_DEV);
    turbotmp::free_array4(mask2dCv_DEV);
    turbotmp::free_array4(bathyT_DEV);
    turbotmp::free_array4(CoriolisBu_DEV);
    turbotmp::free_array4(a_u_DEV);
    turbotmp::free_array4(h_u_DEV);
    turbotmp::free_array4(a_v_DEV);
    turbotmp::free_array4(h_v_DEV);
    turbotmp::free_array4(frac_shelf_u_DEV);
    turbotmp::free_array4(frac_shelf_v_DEV);
    turbotmp::free_array4(Ustar_2d_DEV);
    turbotmp::free_array4(Kv_bbl_u_DEV);
    turbotmp::free_array4(Kv_bbl_v_DEV);
    turbotmp::free_array4(bbl_thick_u_DEV);
    turbotmp::free_array4(bbl_thick_v_DEV);
    turbotmp::free_array4(tbl_thick_shelf_u_DEV);
    turbotmp::free_array4(tbl_thick_shelf_v_DEV);
    turbotmp::free_array4(Kv_tbl_shelf_u_DEV);
    turbotmp::free_array4(Kv_tbl_shelf_v_DEV);
    turbotmp::free_array4(nkml_visc_u_DEV);
    turbotmp::free_array4(nkml_visc_v_DEV);
    turbotmp::free_array4(Kv_shear_DEV);
    turbotmp::free_array4(Kv_shear_Bu_DEV);
    turbotmp::free_array4(SpV_avg_DEV);
    turbotmp::free_array4(kappa_gl90_2d_DEV);
    turbotmp::free_array4(kdgl90_struct_DEV);
    turbotmp::free_array4(a_u_gl90_DEV);
    turbotmp::free_array4(a_v_gl90_DEV);
    turbotmp::free_array4(a1_shelf_u_DEV);
    turbotmp::free_array4(a1_shelf_v_DEV);
    turbotmp::free_array4(hML_u_DEV);
    turbotmp::free_array4(hML_v_DEV);
    turbotmp::free_array4(Kv_u_DEV);
    turbotmp::free_array4(Kv_v_DEV);
    turbotmp::free_array4(Kv_gl90_u_DEV);
    turbotmp::free_array4(Kv_gl90_v_DEV);
    if (g_prime_data != nullptr) The_Arena()->free(g_prime_data);
}

/// @brief Bridge for the MOM6 vertvisc_TR subroutine.
///
/// @param bxU0_HOST           The u-point box of the whole column, [IscB:IecB, jsc:jec, 1:nz]
/// @param bxV0_HOST           The v-point box of the whole column, [isc:iec, JscB:JecB, 1:nz]
/// @param u_HOST              Zonal velocity [L T-1 ~> m s-1]
/// @param v_HOST              Meridional velocity [L T-1 ~> m s-1]
/// @param h_HOST              Layer thickness [H ~> m or kg m-2]
/// @param forces_HOST         Copies of the driving mechanical forces; each field may be
///                            unallocated (data == nullptr)
/// @param visc_HOST           Copies of the viscosities and bottom drag; each field may be
///                            unallocated (data == nullptr)
/// @param dt                  Time increment [T ~> s]
/// @param OBC                 Open boundary condition structure
/// @param AD_HOST             Copies of the accelerations in the momentum equations for
///                            diagnostics; each field may be unallocated (data == nullptr)
/// @param Gcore_HOST          Persistent copies of the core ocean grid fields
/// @param vgrid_HOST          The vertical grid and unit-scaling constants
/// @param CS_HOST             Vertical viscosity control structure
/// @param coupling_state_HOST The viscous coupling coefficients and effective layer thicknesses
///                            at velocity points; a_u_gl90, a_v_gl90, a1_shelf_u and a1_shelf_v
///                            may be unallocated (data == nullptr)
/// @param taux_bot_HOST       Zonal bottom stress from ocean to rock [R L Z T-2 ~> Pa]; may be
///                            unallocated (data == nullptr)
/// @param tauy_bot_HOST       Meridional bottom stress from ocean to rock [R L Z T-2 ~> Pa]; may
///                            be unallocated (data == nullptr)
/// @param KE_u_HOST           The area integral of a KE term in a layer at u-points
///                            [H L4 T-3 ~> m5 s-3 or kg m2 s-3]; may be unallocated
/// @param KE_v_HOST           The area integral of a KE term in a layer at v-points
///                            [H L4 T-3 ~> m5 s-3 or kg m2 s-3]; may be unallocated
/// @param report_HOST         The velocities and flags used to report the columns in which
///                            velocity truncations occur; each field may be unallocated
///                            (data == nullptr)
/// @param ntrunc_sum          The number of velocity truncations, added to the value on entry
///                            [nondim]
/// @param fpmix               fpmix along Eulerian shear
/// @param Us_x_HOST           3d zonal Stokes drift profile [L T-1 ~> m s-1]; may be unallocated
/// @param Us_y_HOST           3d meridional Stokes drift profile [L T-1 ~> m s-1]; may be
///                            unallocated
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
                                 const RealArray_C* Us_y_HOST)
{
    /// Reconstruct the iteration Boxes (Fortran 1-based -> AMReX 0-based)
    Box bxU0 = turbotmp::to_box(bxU0_HOST);
    Box bxV0 = turbotmp::to_box(bxV0_HOST);

    /// Allocate one A4Box per array, sized by RealArray_C::shape[] and
    /// positioned by RealArray_C::lb[]
    auto u_DEV               = turbotmp::make_array4(u_HOST);
    auto v_DEV               = turbotmp::make_array4(v_HOST);
    auto h_DEV               = turbotmp::make_array4(h_HOST);
    auto mask2dCu_DEV        = turbotmp::make_array4(&Gcore_HOST->mask2dCu);
    auto mask2dCv_DEV        = turbotmp::make_array4(&Gcore_HOST->mask2dCv);
    auto areaCu_DEV          = turbotmp::make_array4(&Gcore_HOST->areaCu);
    auto areaCv_DEV          = turbotmp::make_array4(&Gcore_HOST->areaCv);
    auto areaT_DEV           = turbotmp::make_array4(&Gcore_HOST->areaT);
    auto IareaT_DEV          = turbotmp::make_array4(&Gcore_HOST->IareaT);
    auto dx_Cv_DEV           = turbotmp::make_array4(&Gcore_HOST->dx_Cv);
    auto dy_Cu_DEV           = turbotmp::make_array4(&Gcore_HOST->dy_Cu);
    auto a_u_DEV             = turbotmp::make_array4(&coupling_state_HOST->a_u);
    auto h_u_DEV             = turbotmp::make_array4(&coupling_state_HOST->h_u);
    auto a_v_DEV             = turbotmp::make_array4(&coupling_state_HOST->a_v);
    auto h_v_DEV             = turbotmp::make_array4(&coupling_state_HOST->h_v);

    /// Fields that may be unallocated (data == nullptr); only allocate/copy them when allocated.
    const bool taux_allocated            = (forces_HOST->taux.data != nullptr);
    const bool tauy_allocated            = (forces_HOST->tauy.data != nullptr);
    const bool Ray_u_allocated           = (visc_HOST->Ray_u.data != nullptr);
    const bool Ray_v_allocated           = (visc_HOST->Ray_v.data != nullptr);
    const bool a_u_gl90_allocated        = (coupling_state_HOST->a_u_gl90.data != nullptr);
    const bool a_v_gl90_allocated        = (coupling_state_HOST->a_v_gl90.data != nullptr);
    const bool a1_shelf_u_allocated      = (coupling_state_HOST->a1_shelf_u.data != nullptr);
    const bool a1_shelf_v_allocated      = (coupling_state_HOST->a1_shelf_v.data != nullptr);
    const bool Us_x_allocated            = (Us_x_HOST->data != nullptr);
    const bool Us_y_allocated            = (Us_y_HOST->data != nullptr);
    const bool taux_shelf_allocated      = (visc_HOST->taux_shelf.data != nullptr);
    const bool tauy_shelf_allocated      = (visc_HOST->tauy_shelf.data != nullptr);
    const bool du_dt_visc_allocated      = (AD_HOST->du_dt_visc.data != nullptr);
    const bool dv_dt_visc_allocated      = (AD_HOST->dv_dt_visc.data != nullptr);
    const bool du_dt_visc_gl90_allocated = (AD_HOST->du_dt_visc_gl90.data != nullptr);
    const bool dv_dt_visc_gl90_allocated = (AD_HOST->dv_dt_visc_gl90.data != nullptr);
    const bool du_dt_str_allocated       = (AD_HOST->du_dt_str.data != nullptr);
    const bool dv_dt_str_allocated       = (AD_HOST->dv_dt_str.data != nullptr);
    const bool taux_bot_allocated        = (taux_bot_HOST->data != nullptr);
    const bool tauy_bot_allocated        = (tauy_bot_HOST->data != nullptr);
    const bool KE_u_allocated            = (KE_u_HOST->data != nullptr);
    const bool KE_v_allocated            = (KE_v_HOST->data != nullptr);
    const bool u_old_allocated           = (report_HOST->u_old.data != nullptr);
    const bool v_old_allocated           = (report_HOST->v_old.data != nullptr);
    const bool vel_report_u_allocated    = (report_HOST->vel_report_u.data != nullptr);
    const bool vel_report_v_allocated    = (report_HOST->vel_report_v.data != nullptr);
    const bool dowrite_u_allocated       = (report_HOST->dowrite_u.data != nullptr);
    const bool dowrite_v_allocated       = (report_HOST->dowrite_v.data != nullptr);
    turbotmp::A4Box taux_DEV{};
    turbotmp::A4Box tauy_DEV{};
    turbotmp::A4Box Ray_u_DEV{};
    turbotmp::A4Box Ray_v_DEV{};
    turbotmp::A4Box a_u_gl90_DEV{};
    turbotmp::A4Box a_v_gl90_DEV{};
    turbotmp::A4Box a1_shelf_u_DEV{};
    turbotmp::A4Box a1_shelf_v_DEV{};
    turbotmp::A4Box Us_x_DEV{};
    turbotmp::A4Box Us_y_DEV{};
    turbotmp::A4Box taux_shelf_DEV{};
    turbotmp::A4Box tauy_shelf_DEV{};
    turbotmp::A4Box du_dt_visc_DEV{};
    turbotmp::A4Box dv_dt_visc_DEV{};
    turbotmp::A4Box du_dt_visc_gl90_DEV{};
    turbotmp::A4Box dv_dt_visc_gl90_DEV{};
    turbotmp::A4Box du_dt_str_DEV{};
    turbotmp::A4Box dv_dt_str_DEV{};
    turbotmp::A4Box taux_bot_DEV{};
    turbotmp::A4Box tauy_bot_DEV{};
    turbotmp::A4Box KE_u_DEV{};
    turbotmp::A4Box KE_v_DEV{};
    turbotmp::A4Box u_old_DEV{};
    turbotmp::A4Box v_old_DEV{};
    turbotmp::A4Box vel_report_u_DEV{};
    turbotmp::A4Box vel_report_v_DEV{};
    turbotmp::IntA4Box dowrite_u_DEV{};
    turbotmp::IntA4Box dowrite_v_DEV{};
    if (taux_allocated) taux_DEV = turbotmp::make_array4(&forces_HOST->taux);
    if (tauy_allocated) tauy_DEV = turbotmp::make_array4(&forces_HOST->tauy);
    if (Ray_u_allocated) Ray_u_DEV = turbotmp::make_array4(&visc_HOST->Ray_u);
    if (Ray_v_allocated) Ray_v_DEV = turbotmp::make_array4(&visc_HOST->Ray_v);
    if (a_u_gl90_allocated) a_u_gl90_DEV = turbotmp::make_array4(&coupling_state_HOST->a_u_gl90);
    if (a_v_gl90_allocated) a_v_gl90_DEV = turbotmp::make_array4(&coupling_state_HOST->a_v_gl90);
    if (a1_shelf_u_allocated) a1_shelf_u_DEV = turbotmp::make_array4(&coupling_state_HOST->a1_shelf_u);
    if (a1_shelf_v_allocated) a1_shelf_v_DEV = turbotmp::make_array4(&coupling_state_HOST->a1_shelf_v);
    if (Us_x_allocated) Us_x_DEV = turbotmp::make_array4(Us_x_HOST);
    if (Us_y_allocated) Us_y_DEV = turbotmp::make_array4(Us_y_HOST);
    if (taux_shelf_allocated) taux_shelf_DEV = turbotmp::make_array4(&visc_HOST->taux_shelf);
    if (tauy_shelf_allocated) tauy_shelf_DEV = turbotmp::make_array4(&visc_HOST->tauy_shelf);
    if (du_dt_visc_allocated) du_dt_visc_DEV = turbotmp::make_array4(&AD_HOST->du_dt_visc);
    if (dv_dt_visc_allocated) dv_dt_visc_DEV = turbotmp::make_array4(&AD_HOST->dv_dt_visc);
    if (du_dt_visc_gl90_allocated) du_dt_visc_gl90_DEV = turbotmp::make_array4(&AD_HOST->du_dt_visc_gl90);
    if (dv_dt_visc_gl90_allocated) dv_dt_visc_gl90_DEV = turbotmp::make_array4(&AD_HOST->dv_dt_visc_gl90);
    if (du_dt_str_allocated) du_dt_str_DEV = turbotmp::make_array4(&AD_HOST->du_dt_str);
    if (dv_dt_str_allocated) dv_dt_str_DEV = turbotmp::make_array4(&AD_HOST->dv_dt_str);
    if (taux_bot_allocated) taux_bot_DEV = turbotmp::make_array4(taux_bot_HOST);
    if (tauy_bot_allocated) tauy_bot_DEV = turbotmp::make_array4(tauy_bot_HOST);
    if (KE_u_allocated) KE_u_DEV = turbotmp::make_array4(KE_u_HOST);
    if (KE_v_allocated) KE_v_DEV = turbotmp::make_array4(KE_v_HOST);
    if (u_old_allocated) u_old_DEV = turbotmp::make_array4(&report_HOST->u_old);
    if (v_old_allocated) v_old_DEV = turbotmp::make_array4(&report_HOST->v_old);
    if (vel_report_u_allocated) vel_report_u_DEV = turbotmp::make_array4(&report_HOST->vel_report_u);
    if (vel_report_v_allocated) vel_report_v_DEV = turbotmp::make_array4(&report_HOST->vel_report_v);
    if (dowrite_u_allocated) dowrite_u_DEV = turbotmp::make_int_array4(&report_HOST->dowrite_u);
    if (dowrite_v_allocated) dowrite_v_DEV = turbotmp::make_int_array4(&report_HOST->dowrite_v);

    /// Copy host -> device for every array, including the inouts
    turbotmp::copy_FortranHost_to_array4(u_HOST->data, u_DEV);
    turbotmp::copy_FortranHost_to_array4(v_HOST->data, v_DEV);
    turbotmp::copy_FortranHost_to_array4(h_HOST->data, h_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dCu.data, mask2dCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->mask2dCv.data, mask2dCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->areaCu.data, areaCu_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->areaCv.data, areaCv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->areaT.data, areaT_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->IareaT.data, IareaT_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->dx_Cv.data, dx_Cv_DEV);
    turbotmp::copy_FortranHost_to_array4(Gcore_HOST->dy_Cu.data, dy_Cu_DEV);
    turbotmp::copy_FortranHost_to_array4(coupling_state_HOST->a_u.data, a_u_DEV);
    turbotmp::copy_FortranHost_to_array4(coupling_state_HOST->h_u.data, h_u_DEV);
    turbotmp::copy_FortranHost_to_array4(coupling_state_HOST->a_v.data, a_v_DEV);
    turbotmp::copy_FortranHost_to_array4(coupling_state_HOST->h_v.data, h_v_DEV);
    if (taux_allocated) turbotmp::copy_FortranHost_to_array4(forces_HOST->taux.data, taux_DEV);
    if (tauy_allocated) turbotmp::copy_FortranHost_to_array4(forces_HOST->tauy.data, tauy_DEV);
    if (Ray_u_allocated) turbotmp::copy_FortranHost_to_array4(visc_HOST->Ray_u.data, Ray_u_DEV);
    if (Ray_v_allocated) turbotmp::copy_FortranHost_to_array4(visc_HOST->Ray_v.data, Ray_v_DEV);
    if (a_u_gl90_allocated) turbotmp::copy_FortranHost_to_array4(coupling_state_HOST->a_u_gl90.data, a_u_gl90_DEV);
    if (a_v_gl90_allocated) turbotmp::copy_FortranHost_to_array4(coupling_state_HOST->a_v_gl90.data, a_v_gl90_DEV);
    if (a1_shelf_u_allocated) turbotmp::copy_FortranHost_to_array4(coupling_state_HOST->a1_shelf_u.data, a1_shelf_u_DEV);
    if (a1_shelf_v_allocated) turbotmp::copy_FortranHost_to_array4(coupling_state_HOST->a1_shelf_v.data, a1_shelf_v_DEV);
    if (Us_x_allocated) turbotmp::copy_FortranHost_to_array4(Us_x_HOST->data, Us_x_DEV);
    if (Us_y_allocated) turbotmp::copy_FortranHost_to_array4(Us_y_HOST->data, Us_y_DEV);
    if (taux_shelf_allocated) turbotmp::copy_FortranHost_to_array4(visc_HOST->taux_shelf.data, taux_shelf_DEV);
    if (tauy_shelf_allocated) turbotmp::copy_FortranHost_to_array4(visc_HOST->tauy_shelf.data, tauy_shelf_DEV);
    if (du_dt_visc_allocated) turbotmp::copy_FortranHost_to_array4(AD_HOST->du_dt_visc.data, du_dt_visc_DEV);
    if (dv_dt_visc_allocated) turbotmp::copy_FortranHost_to_array4(AD_HOST->dv_dt_visc.data, dv_dt_visc_DEV);
    if (du_dt_visc_gl90_allocated) turbotmp::copy_FortranHost_to_array4(AD_HOST->du_dt_visc_gl90.data, du_dt_visc_gl90_DEV);
    if (dv_dt_visc_gl90_allocated) turbotmp::copy_FortranHost_to_array4(AD_HOST->dv_dt_visc_gl90.data, dv_dt_visc_gl90_DEV);
    if (du_dt_str_allocated) turbotmp::copy_FortranHost_to_array4(AD_HOST->du_dt_str.data, du_dt_str_DEV);
    if (dv_dt_str_allocated) turbotmp::copy_FortranHost_to_array4(AD_HOST->dv_dt_str.data, dv_dt_str_DEV);
    if (taux_bot_allocated) turbotmp::copy_FortranHost_to_array4(taux_bot_HOST->data, taux_bot_DEV);
    if (tauy_bot_allocated) turbotmp::copy_FortranHost_to_array4(tauy_bot_HOST->data, tauy_bot_DEV);
    if (KE_u_allocated) turbotmp::copy_FortranHost_to_array4(KE_u_HOST->data, KE_u_DEV);
    if (KE_v_allocated) turbotmp::copy_FortranHost_to_array4(KE_v_HOST->data, KE_v_DEV);
    if (u_old_allocated) turbotmp::copy_FortranHost_to_array4(report_HOST->u_old.data, u_old_DEV);
    if (v_old_allocated) turbotmp::copy_FortranHost_to_array4(report_HOST->v_old.data, v_old_DEV);
    if (vel_report_u_allocated) turbotmp::copy_FortranHost_to_array4(report_HOST->vel_report_u.data, vel_report_u_DEV);
    if (vel_report_v_allocated) turbotmp::copy_FortranHost_to_array4(report_HOST->vel_report_v.data, vel_report_v_DEV);
    if (dowrite_u_allocated) turbotmp::copy_FortranHost_to_int_array4(report_HOST->dowrite_u.data, dowrite_u_DEV);
    if (dowrite_v_allocated) turbotmp::copy_FortranHost_to_int_array4(report_HOST->dowrite_v.data, dowrite_v_DEV);

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::vertvisc_tr(bxU0, bxV0,
                     u_DEV.arr, v_DEV.arr, h_DEV.arr,
                     taux_DEV.arr, tauy_DEV.arr,
                     Ray_u_DEV.arr, Ray_v_DEV.arr,
                     taux_shelf_DEV.arr, tauy_shelf_DEV.arr,
                     dt,
                     OBC,
                     du_dt_visc_DEV.arr, dv_dt_visc_DEV.arr,
                     du_dt_visc_gl90_DEV.arr, dv_dt_visc_gl90_DEV.arr,
                     du_dt_str_DEV.arr, dv_dt_str_DEV.arr,
                     mask2dCu_DEV.arr, mask2dCv_DEV.arr,
                     areaCu_DEV.arr, areaCv_DEV.arr,
                     areaT_DEV.arr, IareaT_DEV.arr,
                     dx_Cv_DEV.arr, dy_Cu_DEV.arr,
                     vgrid_HOST->ke, vgrid_HOST->H_subroundoff, vgrid_HOST->Angstrom_H,
                     vgrid_HOST->H_to_RZ, vgrid_HOST->m_s_to_L_T,
                     *CS_HOST,
                     a_u_DEV.arr, a_u_gl90_DEV.arr, h_u_DEV.arr,
                     a_v_DEV.arr, a_v_gl90_DEV.arr, h_v_DEV.arr,
                     a1_shelf_u_DEV.arr, a1_shelf_v_DEV.arr,
                     taux_bot_DEV.arr, tauy_bot_DEV.arr,
                     KE_u_DEV.arr, KE_v_DEV.arr,
                     u_old_DEV.arr, v_old_DEV.arr,
                     vel_report_u_DEV.arr, vel_report_v_DEV.arr,
                     dowrite_u_DEV.arr, dowrite_v_DEV.arr,
                     *ntrunc_sum,
                     fpmix,
                     Us_x_DEV.arr, Us_y_DEV.arr);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(u_DEV, u_HOST->data);
    turbotmp::copy_array4_to_FortranHost(v_DEV, v_HOST->data);
    if (taux_shelf_allocated) turbotmp::copy_array4_to_FortranHost(taux_shelf_DEV, visc_HOST->taux_shelf.data);
    if (tauy_shelf_allocated) turbotmp::copy_array4_to_FortranHost(tauy_shelf_DEV, visc_HOST->tauy_shelf.data);
    if (du_dt_visc_allocated) turbotmp::copy_array4_to_FortranHost(du_dt_visc_DEV, AD_HOST->du_dt_visc.data);
    if (dv_dt_visc_allocated) turbotmp::copy_array4_to_FortranHost(dv_dt_visc_DEV, AD_HOST->dv_dt_visc.data);
    if (du_dt_visc_gl90_allocated) turbotmp::copy_array4_to_FortranHost(du_dt_visc_gl90_DEV, AD_HOST->du_dt_visc_gl90.data);
    if (dv_dt_visc_gl90_allocated) turbotmp::copy_array4_to_FortranHost(dv_dt_visc_gl90_DEV, AD_HOST->dv_dt_visc_gl90.data);
    if (du_dt_str_allocated) turbotmp::copy_array4_to_FortranHost(du_dt_str_DEV, AD_HOST->du_dt_str.data);
    if (dv_dt_str_allocated) turbotmp::copy_array4_to_FortranHost(dv_dt_str_DEV, AD_HOST->dv_dt_str.data);
    if (taux_bot_allocated) turbotmp::copy_array4_to_FortranHost(taux_bot_DEV, taux_bot_HOST->data);
    if (tauy_bot_allocated) turbotmp::copy_array4_to_FortranHost(tauy_bot_DEV, tauy_bot_HOST->data);
    if (KE_u_allocated) turbotmp::copy_array4_to_FortranHost(KE_u_DEV, KE_u_HOST->data);
    if (KE_v_allocated) turbotmp::copy_array4_to_FortranHost(KE_v_DEV, KE_v_HOST->data);
    if (u_old_allocated) turbotmp::copy_array4_to_FortranHost(u_old_DEV, report_HOST->u_old.data);
    if (v_old_allocated) turbotmp::copy_array4_to_FortranHost(v_old_DEV, report_HOST->v_old.data);
    if (vel_report_u_allocated) turbotmp::copy_array4_to_FortranHost(vel_report_u_DEV, report_HOST->vel_report_u.data);
    if (vel_report_v_allocated) turbotmp::copy_array4_to_FortranHost(vel_report_v_DEV, report_HOST->vel_report_v.data);
    if (dowrite_u_allocated) turbotmp::copy_int_array4_to_FortranHost(dowrite_u_DEV, report_HOST->dowrite_u.data);
    if (dowrite_v_allocated) turbotmp::copy_int_array4_to_FortranHost(dowrite_v_DEV, report_HOST->dowrite_v.data);

    /// Free everything we allocated (freeing a never-allocated container is a no-op)
    turbotmp::free_array4(u_DEV);
    turbotmp::free_array4(v_DEV);
    turbotmp::free_array4(h_DEV);
    turbotmp::free_array4(mask2dCu_DEV);
    turbotmp::free_array4(mask2dCv_DEV);
    turbotmp::free_array4(areaCu_DEV);
    turbotmp::free_array4(areaCv_DEV);
    turbotmp::free_array4(areaT_DEV);
    turbotmp::free_array4(IareaT_DEV);
    turbotmp::free_array4(dx_Cv_DEV);
    turbotmp::free_array4(dy_Cu_DEV);
    turbotmp::free_array4(a_u_DEV);
    turbotmp::free_array4(h_u_DEV);
    turbotmp::free_array4(a_v_DEV);
    turbotmp::free_array4(h_v_DEV);
    turbotmp::free_array4(taux_DEV);
    turbotmp::free_array4(tauy_DEV);
    turbotmp::free_array4(Ray_u_DEV);
    turbotmp::free_array4(Ray_v_DEV);
    turbotmp::free_array4(a_u_gl90_DEV);
    turbotmp::free_array4(a_v_gl90_DEV);
    turbotmp::free_array4(a1_shelf_u_DEV);
    turbotmp::free_array4(a1_shelf_v_DEV);
    turbotmp::free_array4(Us_x_DEV);
    turbotmp::free_array4(Us_y_DEV);
    turbotmp::free_array4(taux_shelf_DEV);
    turbotmp::free_array4(tauy_shelf_DEV);
    turbotmp::free_array4(du_dt_visc_DEV);
    turbotmp::free_array4(dv_dt_visc_DEV);
    turbotmp::free_array4(du_dt_visc_gl90_DEV);
    turbotmp::free_array4(dv_dt_visc_gl90_DEV);
    turbotmp::free_array4(du_dt_str_DEV);
    turbotmp::free_array4(dv_dt_str_DEV);
    turbotmp::free_array4(taux_bot_DEV);
    turbotmp::free_array4(tauy_bot_DEV);
    turbotmp::free_array4(KE_u_DEV);
    turbotmp::free_array4(KE_v_DEV);
    turbotmp::free_array4(u_old_DEV);
    turbotmp::free_array4(v_old_DEV);
    turbotmp::free_array4(vel_report_u_DEV);
    turbotmp::free_array4(vel_report_v_DEV);
    turbotmp::free_int_array4(dowrite_u_DEV);
    turbotmp::free_int_array4(dowrite_v_DEV);
}
