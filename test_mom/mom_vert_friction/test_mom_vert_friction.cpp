// Unit tests for MOM::vertvisc_limit_vel / vertvisc_remnant_tr / vertfpmix_tr /
// vertvisc_coef_tr / vertvisc_tr.
//
// Each test loads a captured Fortran (input, expected-output) pair from
// <data-dir>/<name>.{bin,meta}, runs the C++ kernel over equivalent AMReX
// containers, and compares the result against the captured "after" arrays.

// SKILLS: 0.3.3

#include <gtest/gtest.h>

#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>

#include "amrex_assertions.hpp"
#include "captured_io.hpp"
#include "data_dir.hpp"
#include "mom_vert_friction.hpp"

using test_mom::expect_arrays_equal;
using test_mom::expect_int_arrays_equal;
using test_mom::to_host_fab;
using test_mom::to_host_int_fab;

namespace {

// A truncation-report array that is recorded only when it was allocated at
// capture time ("<prefix>_allocated"). When allocated, its pre- and post-call
// states are "<prefix>_before" / "<prefix>_after"; when unallocated, `arr`
// stays default-constructed (.p == nullptr).
struct AllocatableReport {
    amrex::FArrayBox before_fab;
    amrex::FArrayBox after_fab;
    amrex::Array4<amrex::Real> arr{};
    bool allocated = false;
};

AllocatableReport bind_allocatable_report(const test_mom::CapturedFile& captured,
                                    const std::string& prefix) {
    AllocatableReport d;
    d.allocated = captured.logical(prefix + "_allocated");
    if (d.allocated) {
        d.before_fab = captured.fab_device(prefix + "_before");
        d.arr = d.before_fab.array();
        d.after_fab = captured.fab_host(prefix + "_after");
    }
    return d;
}

// The integer-encoded (0/1) counterpart of AllocatableReport, for the column
// flags captured as LogicalArray_t.
struct AllocatableIntReport {
    amrex::IArrayBox before_fab;
    amrex::IArrayBox after_fab;
    amrex::Array4<int> arr{};
    bool allocated = false;
};

AllocatableIntReport bind_allocatable_int_report(const test_mom::CapturedFile& captured,
                                           const std::string& prefix) {
    AllocatableIntReport d;
    d.allocated = captured.logical(prefix + "_allocated");
    if (d.allocated) {
        d.before_fab = captured.int_fab_device(prefix + "_before");
        d.arr = d.before_fab.array();
        d.after_fab = captured.int_fab_host(prefix + "_after");
    }
    return d;
}

} // namespace

// -------------------------------------------------------------------------
// vertvisc_limit_vel
// -------------------------------------------------------------------------
// Field mapping: MOM6 src/parameterizations/vertical/MOM_vert_friction.F90:6422-6451 (vertvisc_limit_vel capture block).
TEST(VertviscLimitVel, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "vertvisc_limit_vel");

    const auto   bxU0       = captured.box("_bxU0");
    const auto   bxV0       = captured.box("_bxV0");
    auto         u          = captured.fab_device("_u_before");
    auto         v          = captured.fab_device("_v_before");
    const auto   u_after    = captured.fab_host("_u_after");
    const auto   v_after    = captured.fab_host("_v_after");
    const double dt         = captured.real64("_dt");
    const auto   areaT      = captured.fab_device("_areaT");
    const auto   IareaT     = captured.fab_device("_IareaT");
    const auto   dx_Cv      = captured.fab_device("_dx_Cv");
    const auto   dy_Cu      = captured.fab_device("_dy_Cu");
    const int    ke         = captured.integer("_vgrid_ke");
    const double Angstrom_H = captured.real64("_vgrid_Angstrom_H");
    const double m_s_to_L_T = captured.real64("_vgrid_m_s_to_L_T");

    vertvisc_opts_CS_C CS{};
    CS.Hmix_stress        = captured.real64("_CS_Hmix_stress");
    CS.direct_stress      = captured.logical("_CS_direct_stress");
    CS.StokesMixing       = captured.logical("_CS_StokesMixing");
    CS.vel_underflow      = captured.real64("_CS_vel_underflow");
    CS.CFL_trunc          = captured.real64("_CS_CFL_trunc");
    CS.CFL_report         = captured.real64("_CS_CFL_report");
    CS.do_report_u        = captured.logical("_CS_do_report_u");
    CS.do_report_v        = captured.logical("_CS_do_report_v");
    CS.do_GLwork          = captured.logical("_CS_do_GLwork");
    CS.do_du_dt_visc_gl90 = captured.logical("_CS_do_du_dt_visc_gl90");
    CS.do_dv_dt_visc_gl90 = captured.logical("_CS_do_dv_dt_visc_gl90");

    // The effective layer thicknesses at u- and v-points [H ~> m or kg m-2]
    ASSERT_TRUE(captured.logical("_coupling_state_h_u_allocated"));
    ASSERT_TRUE(captured.logical("_coupling_state_h_v_allocated"));
    const auto   h_u        = captured.fab_device("_coupling_state_h_u_before");
    const auto   h_v        = captured.fab_device("_coupling_state_h_v_before");

    auto u_old        = bind_allocatable_report(captured, "_report_u_old");
    auto v_old        = bind_allocatable_report(captured, "_report_v_old");
    auto vel_report_u = bind_allocatable_report(captured, "_report_vel_report_u");
    auto vel_report_v = bind_allocatable_report(captured, "_report_vel_report_v");
    auto dowrite_u    = bind_allocatable_int_report(captured, "_report_dowrite_u");
    auto dowrite_v    = bind_allocatable_int_report(captured, "_report_dowrite_v");

    int          ntrunc_sum       = captured.integer("_ntrunc_sum_before");
    const int    ntrunc_sum_after = captured.integer("_ntrunc_sum_after");

    MOM::vertvisc_limit_vel(bxU0, bxV0,
                            u.array(),
                            v.array(),
                            dt,
                            areaT.const_array(),
                            IareaT.const_array(),
                            dx_Cv.const_array(),
                            dy_Cu.const_array(),
                            ke,
                            Angstrom_H,
                            m_s_to_L_T,
                            CS,
                            h_u.const_array(),
                            h_v.const_array(),
                            u_old.arr,
                            v_old.arr,
                            vel_report_u.arr,
                            vel_report_v.arr,
                            dowrite_u.arr,
                            dowrite_v.arr,
                            ntrunc_sum);
    amrex::Gpu::synchronize();

    expect_arrays_equal(u_after, to_host_fab(u), "u");
    expect_arrays_equal(v_after, to_host_fab(v), "v");
    EXPECT_EQ(ntrunc_sum, ntrunc_sum_after) << "ntrunc_sum";
    if (u_old.allocated) {
        expect_arrays_equal(u_old.after_fab, to_host_fab(u_old.before_fab), "u_old");
    }
    if (v_old.allocated) {
        expect_arrays_equal(v_old.after_fab, to_host_fab(v_old.before_fab), "v_old");
    }
    if (vel_report_u.allocated) {
        expect_arrays_equal(vel_report_u.after_fab, to_host_fab(vel_report_u.before_fab),
                            "vel_report_u");
    }
    if (vel_report_v.allocated) {
        expect_arrays_equal(vel_report_v.after_fab, to_host_fab(vel_report_v.before_fab),
                            "vel_report_v");
    }
    if (dowrite_u.allocated) {
        expect_int_arrays_equal(dowrite_u.after_fab, to_host_int_fab(dowrite_u.before_fab),
                                "dowrite_u");
    }
    if (dowrite_v.allocated) {
        expect_int_arrays_equal(dowrite_v.after_fab, to_host_int_fab(dowrite_v.before_fab),
                                "dowrite_v");
    }
}

// -------------------------------------------------------------------------
// vertvisc_remnant_TR
// -------------------------------------------------------------------------
// Field mapping: MOM6 src/parameterizations/vertical/MOM_vert_friction.F90:3709-3739 (vertvisc_remnant_TR capture block).
TEST(VertviscRemnantTr, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "vertvisc_remnant_tr");

    const auto   bxU0             = captured.box("_bxU0");
    const auto   bxV0             = captured.box("_bxV0");

    // The Rayleigh drag velocities at u- and v-points [H T-1 ~> m s-1 or Pa s m-1],
    // recorded only when allocated at capture time
    const bool   Ray_u_allocated    = captured.logical("_visc_Ray_u_allocated");
    const bool   Ray_v_allocated    = captured.logical("_visc_Ray_v_allocated");
    amrex::FArrayBox Ray_u_fab;
    amrex::FArrayBox Ray_v_fab;
    amrex::Array4<const amrex::Real> Ray_u{};
    amrex::Array4<const amrex::Real> Ray_v{};
    if (Ray_u_allocated) {
        Ray_u_fab = captured.fab_device("_visc_Ray_u_before");
        Ray_u = Ray_u_fab.const_array();
    }
    if (Ray_v_allocated) {
        Ray_v_fab = captured.fab_device("_visc_Ray_v_before");
        Ray_v = Ray_v_fab.const_array();
    }

    auto         visc_rem_u       = captured.fab_device("_visc_rem_u_before");
    auto         visc_rem_v       = captured.fab_device("_visc_rem_v_before");
    const auto   visc_rem_u_after = captured.fab_host("_visc_rem_u_after");
    const auto   visc_rem_v_after = captured.fab_host("_visc_rem_v_after");
    const double dt               = captured.real64("_dt");
    const auto   mask2dCu         = captured.fab_device("_mask2dCu");
    const auto   mask2dCv         = captured.fab_device("_mask2dCv");
    const int    ke               = captured.integer("_ke");

    // The drag coefficients across interfaces [H T-1 ~> m s-1 or Pa s m-1] and the
    // effective layer thicknesses [H ~> m or kg m-2] at u- and v-points
    ASSERT_TRUE(captured.logical("_coupling_state_a_u_allocated"));
    ASSERT_TRUE(captured.logical("_coupling_state_a_v_allocated"));
    ASSERT_TRUE(captured.logical("_coupling_state_h_u_allocated"));
    ASSERT_TRUE(captured.logical("_coupling_state_h_v_allocated"));
    const auto   a_u              = captured.fab_device("_coupling_state_a_u_before");
    const auto   a_v              = captured.fab_device("_coupling_state_a_v_before");
    const auto   h_u              = captured.fab_device("_coupling_state_h_u_before");
    const auto   h_v              = captured.fab_device("_coupling_state_h_v_before");

    MOM::vertvisc_remnant_tr(bxU0, bxV0,
                             Ray_u,
                             Ray_v,
                             visc_rem_u.array(),
                             visc_rem_v.array(),
                             dt,
                             mask2dCu.const_array(),
                             mask2dCv.const_array(),
                             ke,
                             a_u.const_array(),
                             a_v.const_array(),
                             h_u.const_array(),
                             h_v.const_array());
    amrex::Gpu::synchronize();

    expect_arrays_equal(visc_rem_u_after, to_host_fab(visc_rem_u), "visc_rem_u");
    expect_arrays_equal(visc_rem_v_after, to_host_fab(visc_rem_v), "visc_rem_v");
}

// -------------------------------------------------------------------------
// vertFPmix_TR
// -------------------------------------------------------------------------
// Field mapping: MOM6 src/parameterizations/vertical/MOM_vert_friction.F90:1171-1209 (vertFPmix_TR capture block).
TEST(VertfpmixTr, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "vertfpmix_tr");

    const auto   bxU0         = captured.box("_bxU0");
    const auto   bxV0         = captured.box("_bxV0");
    auto         ui           = captured.fab_device("_ui_before");
    auto         vi           = captured.fab_device("_vi_before");
    const auto   ui_after     = captured.fab_host("_ui_after");
    const auto   vi_after     = captured.fab_host("_vi_after");
    const auto   hbl_h        = captured.fab_device("_hbl_h");

    // The zonal and meridional wind stresses [R L Z T-2 ~> Pa]
    ASSERT_TRUE(captured.logical("_forces_taux_allocated"));
    ASSERT_TRUE(captured.logical("_forces_tauy_allocated"));
    const auto   taux         = captured.fab_device("_forces_taux_before");
    const auto   tauy         = captured.fab_device("_forces_tauy_before");

    const auto   Us_x         = captured.fab_device("_Us_x");
    const auto   Us_y         = captured.fab_device("_Us_y");
    const double dt           = captured.real64("_dt");
    const double Cemp_NL      = captured.real64("_Cemp_NL");
    const double Irho0        = captured.real64("_Irho0");
    const auto   mask2dCu     = captured.fab_device("_mask2dCu");
    const auto   mask2dCv     = captured.fab_device("_mask2dCv");
    const auto   mask2dT      = captured.fab_device("_mask2dT");

    // The effective layer thicknesses at u- and v-points [H ~> m or kg m-2]
    ASSERT_TRUE(captured.logical("_coupling_state_h_u_allocated"));
    ASSERT_TRUE(captured.logical("_coupling_state_h_v_allocated"));
    const auto   h_u          = captured.fab_device("_coupling_state_h_u_before");
    const auto   h_v          = captured.fab_device("_coupling_state_h_v_before");

    // The Eulerian velocities and the nonlocal increments at u- and v-points [L T-1 ~> m s-1]
    ASSERT_TRUE(captured.logical("_FPdiag_uE_u_allocated"));
    ASSERT_TRUE(captured.logical("_FPdiag_vE_v_allocated"));
    ASSERT_TRUE(captured.logical("_FPdiag_uInc_u_allocated"));
    ASSERT_TRUE(captured.logical("_FPdiag_vInc_v_allocated"));
    auto         uE_u         = captured.fab_device("_FPdiag_uE_u_before");
    auto         vE_v         = captured.fab_device("_FPdiag_vE_v_before");
    auto         uInc_u       = captured.fab_device("_FPdiag_uInc_u_before");
    auto         vInc_v       = captured.fab_device("_FPdiag_vInc_v_before");
    const auto   uE_u_after   = captured.fab_host("_FPdiag_uE_u_after");
    const auto   vE_v_after   = captured.fab_host("_FPdiag_vE_v_after");
    const auto   uInc_u_after = captured.fab_host("_FPdiag_uInc_u_after");
    const auto   vInc_v_after = captured.fab_host("_FPdiag_vInc_v_after");

    MOM::vertfpmix_tr(bxU0, bxV0,
                      ui.array(),
                      vi.array(),
                      hbl_h.const_array(),
                      taux.const_array(),
                      tauy.const_array(),
                      Us_x.const_array(),
                      Us_y.const_array(),
                      dt,
                      Cemp_NL,
                      Irho0,
                      mask2dCu.const_array(),
                      mask2dCv.const_array(),
                      mask2dT.const_array(),
                      h_u.const_array(),
                      h_v.const_array(),
                      uE_u.array(),
                      vE_v.array(),
                      uInc_u.array(),
                      vInc_v.array());
    amrex::Gpu::synchronize();

    expect_arrays_equal(ui_after,     to_host_fab(ui),     "ui");
    expect_arrays_equal(vi_after,     to_host_fab(vi),     "vi");
    expect_arrays_equal(uE_u_after,   to_host_fab(uE_u),   "uE_u");
    expect_arrays_equal(vE_v_after,   to_host_fab(vE_v),   "vE_v");
    expect_arrays_equal(uInc_u_after, to_host_fab(uInc_u), "uInc_u");
    expect_arrays_equal(vInc_v_after, to_host_fab(vInc_v), "vInc_v");
}

namespace {

// An input array that is recorded only when it was allocated at capture time
// (flag entry "<allocated_flag>"). When allocated, its value is the entry
// "<value_name>"; when unallocated, `arr` stays default-constructed
// (.p == nullptr).
struct AllocatableInput {
    amrex::FArrayBox fab;
    amrex::Array4<const amrex::Real> arr{};
    bool allocated = false;
};

AllocatableInput bind_allocatable_input(const test_mom::CapturedFile& captured,
                                        const std::string& allocated_flag,
                                        const std::string& value_name) {
    AllocatableInput d;
    d.allocated = captured.logical(allocated_flag);
    if (d.allocated) {
        d.fab = captured.fab_device(value_name);
        d.arr = d.fab.const_array();
    }
    return d;
}

} // namespace

// -------------------------------------------------------------------------
// vertvisc_coef_TR
// -------------------------------------------------------------------------
// Field mapping: MOM6 src/parameterizations/vertical/MOM_vert_friction.F90:4855-4901 (vertvisc_coef_TR capture block).
TEST(VertviscCoefTr, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "vertvisc_coef_tr");

    // Open boundary condition structure
    if (captured.logical("_OBC_associated")) {
        GTEST_SKIP() << "vertvisc_coef_tr capture has an associated OBC; OBC-inactive only";
    }

    const auto   bxU0       = captured.box("_bxU0");
    const auto   bxV0       = captured.box("_bxV0");
    const auto   u          = captured.fab_device("_u");
    const auto   v          = captured.fab_device("_v");
    const auto   h          = captured.fab_device("_h");
    const auto   dz         = captured.fab_device("_dz");

    auto frac_shelf_u      = bind_allocatable_input(captured, "_forces_frac_shelf_u_allocated", "_forces_frac_shelf_u_before");
    auto frac_shelf_v      = bind_allocatable_input(captured, "_forces_frac_shelf_v_allocated", "_forces_frac_shelf_v_before");
    auto Kv_bbl_u          = bind_allocatable_input(captured, "_visc_kv_bbl_u_allocated", "_visc_kv_bbl_u_before");
    auto Kv_bbl_v          = bind_allocatable_input(captured, "_visc_kv_bbl_v_allocated", "_visc_kv_bbl_v_before");
    auto bbl_thick_u       = bind_allocatable_input(captured, "_visc_bbl_thick_u_allocated", "_visc_bbl_thick_u_before");
    auto bbl_thick_v       = bind_allocatable_input(captured, "_visc_bbl_thick_v_allocated", "_visc_bbl_thick_v_before");
    auto tbl_thick_shelf_u = bind_allocatable_input(captured, "_visc_tbl_thick_shelf_u_allocated", "_visc_tbl_thick_shelf_u_before");
    auto tbl_thick_shelf_v = bind_allocatable_input(captured, "_visc_tbl_thick_shelf_v_allocated", "_visc_tbl_thick_shelf_v_before");
    auto Kv_tbl_shelf_u    = bind_allocatable_input(captured, "_visc_kv_tbl_shelf_u_allocated", "_visc_kv_tbl_shelf_u_before");
    auto Kv_tbl_shelf_v    = bind_allocatable_input(captured, "_visc_kv_tbl_shelf_v_allocated", "_visc_kv_tbl_shelf_v_before");
    auto nkml_visc_u       = bind_allocatable_input(captured, "_visc_nkml_visc_u_allocated", "_visc_nkml_visc_u_before");
    auto nkml_visc_v       = bind_allocatable_input(captured, "_visc_nkml_visc_v_allocated", "_visc_nkml_visc_v_before");
    auto Kv_shear          = bind_allocatable_input(captured, "_visc_Kv_shear_allocated", "_visc_Kv_shear_before");
    auto Kv_shear_Bu       = bind_allocatable_input(captured, "_visc_Kv_shear_Bu_allocated", "_visc_Kv_shear_Bu_before");

    // The wind friction velocity [Z T-1 ~> m s-1]
    amrex::FArrayBox Ustar_2d_fab;
    amrex::Array4<const amrex::Real> Ustar_2d{};
    if (captured.is_associated("_Ustar_2d")) {
        Ustar_2d_fab = captured.fab_device("_Ustar_2d");
        Ustar_2d = Ustar_2d_fab.const_array();
    }

    auto SpV_avg       = bind_allocatable_input(captured, "_SpV_avg_allocated", "_SpV_avg");
    auto kdgl90_struct = bind_allocatable_input(captured, "_kdgl90_struct_allocated", "_kdgl90_struct");
    auto kappa_gl90_2d = bind_allocatable_input(captured, "_CS_GL90_CS_kappa_gl90_2d_allocated",
                                                "_CS_GL90_CS_kappa_gl90_2d");

    const double dt         = captured.real64("_dt");
    const auto   mask2dCu   = captured.fab_device("_mask2dCu");
    const auto   mask2dCv   = captured.fab_device("_mask2dCv");
    const auto   bathyT     = captured.fab_device("_bathyT");
    const auto   CoriolisBu = captured.fab_device("_CoriolisBu");

    // The vertical grid and unit-scaling constants
    const int    ke             = captured.integer("_vgrid_ke");
    const bool   Boussinesq     = captured.logical("_vgrid_Boussinesq");
    const int    nkml           = captured.integer("_vgrid_nkml");
    const double dZ_subroundoff = captured.real64("_vgrid_dZ_subroundoff");
    const double H_subroundoff  = captured.real64("_vgrid_H_subroundoff");
    const double H_to_m         = captured.real64("_vgrid_H_to_m");
    const double m_to_H         = captured.real64("_vgrid_m_to_H");
    const double RZ_to_H        = captured.real64("_vgrid_RZ_to_H");
    const double Z_to_H         = captured.real64("_vgrid_Z_to_H");
    const double T_to_s         = captured.real64("_vgrid_T_to_s");

    // The reduced gravity at each interface [L2 Z-1 T-2 ~> m s-2]
    amrex::FArrayBox g_prime_fab;
    amrex::Array4<const amrex::Real> g_prime{};
    if (captured.is_associated("_vgrid_g_prime")) {
        g_prime_fab = captured.fab_device_1d("_vgrid_g_prime");
        g_prime = g_prime_fab.const_array();
    }

    // Vertical viscosity control structure
    coef_opts_CS_C CS{};
    CS.coupling_CS.Hmix               = captured.real64("_CS_coupling_CS_Hmix");
    CS.coupling_CS.Kvml_invZ2         = captured.real64("_CS_coupling_CS_Kvml_invZ2");
    CS.coupling_CS.Kv                 = captured.real64("_CS_coupling_CS_Kv");
    CS.coupling_CS.Kv_extra_bbl       = captured.real64("_CS_coupling_CS_Kv_extra_bbl");
    CS.coupling_CS.vonKar             = captured.real64("_CS_coupling_CS_vonKar");
    CS.coupling_CS.bottomdraglaw      = captured.logical("_CS_coupling_CS_bottomdraglaw");
    CS.coupling_CS.dynamic_viscous_ML = captured.logical("_CS_coupling_CS_dynamic_viscous_ML");
    CS.coupling_CS.fixed_LOTW_ML      = captured.logical("_CS_coupling_CS_fixed_LOTW_ML");
    CS.coupling_CS.apply_LOTW_floor   = captured.logical("_CS_coupling_CS_apply_LOTW_floor");
    CS.coupling_CS.answer_date        = captured.integer("_CS_coupling_CS_answer_date");
    CS.Hbbl                           = captured.real64("_CS_Hbbl");
    CS.Hbbl_gl90                      = captured.real64("_CS_Hbbl_gl90");
    CS.use_GL90_in_SSW                = captured.logical("_CS_use_GL90_in_SSW");
    CS.GL90_CS.use_GL90_N2            = captured.logical("_CS_GL90_CS_use_GL90_N2");
    CS.GL90_CS.kappa_gl90             = captured.real64("_CS_GL90_CS_kappa_gl90");
    CS.GL90_CS.read_kappa_gl90        = captured.logical("_CS_GL90_CS_read_kappa_gl90");
    CS.GL90_CS.alpha_gl90             = captured.real64("_CS_GL90_CS_alpha_gl90");
    CS.harmonic_visc                  = captured.logical("_CS_harmonic_visc");
    CS.harm_BL_val                    = captured.real64("_CS_harm_BL_val");

    // The viscous coupling coefficients and effective layer thicknesses at velocity points
    ASSERT_TRUE(captured.logical("_coupling_state_a_u_allocated"));
    ASSERT_TRUE(captured.logical("_coupling_state_h_u_allocated"));
    ASSERT_TRUE(captured.logical("_coupling_state_a_v_allocated"));
    ASSERT_TRUE(captured.logical("_coupling_state_h_v_allocated"));

    auto a_u               = bind_allocatable_report(captured, "_coupling_state_a_u");
    auto a_u_gl90          = bind_allocatable_report(captured, "_coupling_state_a_u_gl90");
    auto h_u               = bind_allocatable_report(captured, "_coupling_state_h_u");
    auto a_v               = bind_allocatable_report(captured, "_coupling_state_a_v");
    auto a_v_gl90          = bind_allocatable_report(captured, "_coupling_state_a_v_gl90");
    auto h_v               = bind_allocatable_report(captured, "_coupling_state_h_v");
    auto a1_shelf_u        = bind_allocatable_report(captured, "_coupling_state_a1_shelf_u");
    auto a1_shelf_v        = bind_allocatable_report(captured, "_coupling_state_a1_shelf_v");
    auto hML_u             = bind_allocatable_report(captured, "_coef_diag_hML_u");
    auto hML_v             = bind_allocatable_report(captured, "_coef_diag_hML_v");
    auto Kv_u              = bind_allocatable_report(captured, "_coef_diag_Kv_u");
    auto Kv_v              = bind_allocatable_report(captured, "_coef_diag_Kv_v");
    auto Kv_gl90_u         = bind_allocatable_report(captured, "_coef_diag_Kv_gl90_u");
    auto Kv_gl90_v         = bind_allocatable_report(captured, "_coef_diag_Kv_gl90_v");

    MOM::vertvisc_coef_tr(bxU0, bxV0,
                          u.const_array(), v.const_array(), h.const_array(), dz.const_array(),
                          frac_shelf_u.arr, frac_shelf_v.arr,
                          Ustar_2d,
                          Kv_bbl_u.arr, Kv_bbl_v.arr,
                          bbl_thick_u.arr, bbl_thick_v.arr,
                          tbl_thick_shelf_u.arr, tbl_thick_shelf_v.arr,
                          Kv_tbl_shelf_u.arr, Kv_tbl_shelf_v.arr,
                          nkml_visc_u.arr, nkml_visc_v.arr,
                          Kv_shear.arr, Kv_shear_Bu.arr,
                          SpV_avg.arr,
                          dt,
                          mask2dCu.const_array(), mask2dCv.const_array(),
                          bathyT.const_array(), CoriolisBu.const_array(),
                          ke, Boussinesq, nkml,
                          dZ_subroundoff, H_subroundoff,
                          H_to_m, m_to_H, RZ_to_H, Z_to_H, T_to_s,
                          g_prime,
                          CS,
                          kappa_gl90_2d.arr,
                          a_u.arr, a_u_gl90.arr, h_u.arr,
                          a_v.arr, a_v_gl90.arr, h_v.arr,
                          a1_shelf_u.arr, a1_shelf_v.arr,
                          nullptr,
                          kdgl90_struct.arr,
                          hML_u.arr, hML_v.arr,
                          Kv_u.arr, Kv_v.arr,
                          Kv_gl90_u.arr, Kv_gl90_v.arr);
    amrex::Gpu::synchronize();

    if (a_u.allocated) {
        expect_arrays_equal(a_u.after_fab, to_host_fab(a_u.before_fab), "a_u");
    }
    if (a_u_gl90.allocated) {
        expect_arrays_equal(a_u_gl90.after_fab, to_host_fab(a_u_gl90.before_fab), "a_u_gl90");
    }
    if (h_u.allocated) {
        expect_arrays_equal(h_u.after_fab, to_host_fab(h_u.before_fab), "h_u");
    }
    if (a_v.allocated) {
        expect_arrays_equal(a_v.after_fab, to_host_fab(a_v.before_fab), "a_v");
    }
    if (a_v_gl90.allocated) {
        expect_arrays_equal(a_v_gl90.after_fab, to_host_fab(a_v_gl90.before_fab), "a_v_gl90");
    }
    if (h_v.allocated) {
        expect_arrays_equal(h_v.after_fab, to_host_fab(h_v.before_fab), "h_v");
    }
    if (a1_shelf_u.allocated) {
        expect_arrays_equal(a1_shelf_u.after_fab, to_host_fab(a1_shelf_u.before_fab), "a1_shelf_u");
    }
    if (a1_shelf_v.allocated) {
        expect_arrays_equal(a1_shelf_v.after_fab, to_host_fab(a1_shelf_v.before_fab), "a1_shelf_v");
    }
    if (hML_u.allocated) {
        expect_arrays_equal(hML_u.after_fab, to_host_fab(hML_u.before_fab), "hML_u");
    }
    if (hML_v.allocated) {
        expect_arrays_equal(hML_v.after_fab, to_host_fab(hML_v.before_fab), "hML_v");
    }
    if (Kv_u.allocated) {
        expect_arrays_equal(Kv_u.after_fab, to_host_fab(Kv_u.before_fab), "Kv_u");
    }
    if (Kv_v.allocated) {
        expect_arrays_equal(Kv_v.after_fab, to_host_fab(Kv_v.before_fab), "Kv_v");
    }
    if (Kv_gl90_u.allocated) {
        expect_arrays_equal(Kv_gl90_u.after_fab, to_host_fab(Kv_gl90_u.before_fab), "Kv_gl90_u");
    }
    if (Kv_gl90_v.allocated) {
        expect_arrays_equal(Kv_gl90_v.after_fab, to_host_fab(Kv_gl90_v.before_fab), "Kv_gl90_v");
    }
}

// -------------------------------------------------------------------------
// vertvisc_TR
// -------------------------------------------------------------------------
// Field mapping: MOM6 src/parameterizations/vertical/MOM_vert_friction.F90:3243-3311 (vertvisc_TR capture block).
TEST(VertviscTr, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "vertvisc_tr");

    // Open boundary condition structure
    if (captured.logical("_OBC_associated")) {
        GTEST_SKIP() << "vertvisc_tr capture has an associated OBC; OBC-inactive only";
    }

    const auto   bxU0       = captured.box("_bxU0");
    const auto   bxV0       = captured.box("_bxV0");
    auto         u          = captured.fab_device("_u_before");
    auto         v          = captured.fab_device("_v_before");
    const auto   u_after    = captured.fab_host("_u_after");
    const auto   v_after    = captured.fab_host("_v_after");
    const auto   h          = captured.fab_device("_h");
    const double dt         = captured.real64("_dt");
    const auto   mask2dCu   = captured.fab_device("_mask2dCu");
    const auto   mask2dCv   = captured.fab_device("_mask2dCv");
    const auto   areaCu     = captured.fab_device("_areaCu");
    const auto   areaCv     = captured.fab_device("_areaCv");
    const auto   areaT      = captured.fab_device("_areaT");
    const auto   IareaT     = captured.fab_device("_IareaT");
    const auto   dx_Cv      = captured.fab_device("_dx_Cv");
    const auto   dy_Cu      = captured.fab_device("_dy_Cu");

    auto taux            = bind_allocatable_input(captured, "_forces_taux_allocated", "_forces_taux_before");
    auto tauy            = bind_allocatable_input(captured, "_forces_tauy_allocated", "_forces_tauy_before");
    auto Ray_u           = bind_allocatable_input(captured, "_visc_Ray_u_allocated", "_visc_Ray_u_before");
    auto Ray_v           = bind_allocatable_input(captured, "_visc_Ray_v_allocated", "_visc_Ray_v_before");
    auto a_u_gl90        = bind_allocatable_input(captured, "_coupling_state_a_u_gl90_allocated", "_coupling_state_a_u_gl90_before");
    auto a_v_gl90        = bind_allocatable_input(captured, "_coupling_state_a_v_gl90_allocated", "_coupling_state_a_v_gl90_before");
    auto a1_shelf_u      = bind_allocatable_input(captured, "_coupling_state_a1_shelf_u_allocated", "_coupling_state_a1_shelf_u_before");
    auto a1_shelf_v      = bind_allocatable_input(captured, "_coupling_state_a1_shelf_v_allocated", "_coupling_state_a1_shelf_v_before");
    auto Us_x            = bind_allocatable_input(captured, "_Us_x_allocated", "_Us_x");
    auto Us_y            = bind_allocatable_input(captured, "_Us_y_allocated", "_Us_y");

    auto taux_shelf      = bind_allocatable_report(captured, "_visc_taux_shelf");
    auto tauy_shelf      = bind_allocatable_report(captured, "_visc_tauy_shelf");
    auto du_dt_visc      = bind_allocatable_report(captured, "_AD_du_dt_visc");
    auto dv_dt_visc      = bind_allocatable_report(captured, "_AD_dv_dt_visc");
    auto du_dt_visc_gl90 = bind_allocatable_report(captured, "_AD_du_dt_visc_gl90");
    auto dv_dt_visc_gl90 = bind_allocatable_report(captured, "_AD_dv_dt_visc_gl90");
    auto du_dt_str       = bind_allocatable_report(captured, "_AD_du_dt_str");
    auto dv_dt_str       = bind_allocatable_report(captured, "_AD_dv_dt_str");
    auto taux_bot        = bind_allocatable_report(captured, "_taux_bot");
    auto tauy_bot        = bind_allocatable_report(captured, "_tauy_bot");
    auto KE_u            = bind_allocatable_report(captured, "_KE_u");
    auto KE_v            = bind_allocatable_report(captured, "_KE_v");
    auto u_old           = bind_allocatable_report(captured, "_report_u_old");
    auto v_old           = bind_allocatable_report(captured, "_report_v_old");
    auto vel_report_u    = bind_allocatable_report(captured, "_report_vel_report_u");
    auto vel_report_v    = bind_allocatable_report(captured, "_report_vel_report_v");
    auto dowrite_u       = bind_allocatable_int_report(captured, "_report_dowrite_u");
    auto dowrite_v       = bind_allocatable_int_report(captured, "_report_dowrite_v");

    // The vertical grid and unit-scaling constants
    const int    ke            = captured.integer("_vgrid_ke");
    const double H_subroundoff = captured.real64("_vgrid_H_subroundoff");
    const double Angstrom_H    = captured.real64("_vgrid_Angstrom_H");
    const double H_to_RZ       = captured.real64("_vgrid_H_to_RZ");
    const double m_s_to_L_T    = captured.real64("_vgrid_m_s_to_L_T");

    // Vertical viscosity control structure
    vertvisc_opts_CS_C CS{};
    CS.Hmix_stress        = captured.real64("_CS_Hmix_stress");
    CS.direct_stress      = captured.logical("_CS_direct_stress");
    CS.StokesMixing       = captured.logical("_CS_StokesMixing");
    CS.vel_underflow      = captured.real64("_CS_vel_underflow");
    CS.CFL_trunc          = captured.real64("_CS_CFL_trunc");
    CS.CFL_report         = captured.real64("_CS_CFL_report");
    CS.do_report_u        = captured.logical("_CS_do_report_u");
    CS.do_report_v        = captured.logical("_CS_do_report_v");
    CS.do_GLwork          = captured.logical("_CS_do_GLwork");
    CS.do_du_dt_visc_gl90 = captured.logical("_CS_do_du_dt_visc_gl90");
    CS.do_dv_dt_visc_gl90 = captured.logical("_CS_do_dv_dt_visc_gl90");

    // The drag coefficients across interfaces [H T-1 ~> m s-1 or Pa s m-1] and the
    // effective layer thicknesses [H ~> m or kg m-2] at u- and v-points
    ASSERT_TRUE(captured.logical("_coupling_state_a_u_allocated"));
    ASSERT_TRUE(captured.logical("_coupling_state_h_u_allocated"));
    ASSERT_TRUE(captured.logical("_coupling_state_a_v_allocated"));
    ASSERT_TRUE(captured.logical("_coupling_state_h_v_allocated"));
    const auto   a_u        = captured.fab_device("_coupling_state_a_u_before");
    const auto   h_u        = captured.fab_device("_coupling_state_h_u_before");
    const auto   a_v        = captured.fab_device("_coupling_state_a_v_before");
    const auto   h_v        = captured.fab_device("_coupling_state_h_v_before");

    int          ntrunc_sum       = captured.integer("_ntrunc_sum_before");
    const int    ntrunc_sum_after = captured.integer("_ntrunc_sum_after");
    const bool   fpmix            = captured.logical("_fpmix");

    MOM::vertvisc_tr(bxU0, bxV0,
                     u.array(), v.array(), h.const_array(),
                     taux.arr, tauy.arr,
                     Ray_u.arr, Ray_v.arr,
                     taux_shelf.arr, tauy_shelf.arr,
                     dt,
                     nullptr,
                     du_dt_visc.arr, dv_dt_visc.arr,
                     du_dt_visc_gl90.arr, dv_dt_visc_gl90.arr,
                     du_dt_str.arr, dv_dt_str.arr,
                     mask2dCu.const_array(), mask2dCv.const_array(),
                     areaCu.const_array(), areaCv.const_array(),
                     areaT.const_array(), IareaT.const_array(),
                     dx_Cv.const_array(), dy_Cu.const_array(),
                     ke, H_subroundoff, Angstrom_H, H_to_RZ, m_s_to_L_T,
                     CS,
                     a_u.const_array(), a_u_gl90.arr, h_u.const_array(),
                     a_v.const_array(), a_v_gl90.arr, h_v.const_array(),
                     a1_shelf_u.arr, a1_shelf_v.arr,
                     taux_bot.arr, tauy_bot.arr,
                     KE_u.arr, KE_v.arr,
                     u_old.arr, v_old.arr,
                     vel_report_u.arr, vel_report_v.arr,
                     dowrite_u.arr, dowrite_v.arr,
                     ntrunc_sum,
                     fpmix,
                     Us_x.arr, Us_y.arr);
    amrex::Gpu::synchronize();

    expect_arrays_equal(u_after, to_host_fab(u), "u");
    expect_arrays_equal(v_after, to_host_fab(v), "v");
    EXPECT_EQ(ntrunc_sum, ntrunc_sum_after) << "ntrunc_sum";
    if (taux_shelf.allocated) {
        expect_arrays_equal(taux_shelf.after_fab, to_host_fab(taux_shelf.before_fab), "taux_shelf");
    }
    if (tauy_shelf.allocated) {
        expect_arrays_equal(tauy_shelf.after_fab, to_host_fab(tauy_shelf.before_fab), "tauy_shelf");
    }
    if (du_dt_visc.allocated) {
        expect_arrays_equal(du_dt_visc.after_fab, to_host_fab(du_dt_visc.before_fab), "du_dt_visc");
    }
    if (dv_dt_visc.allocated) {
        expect_arrays_equal(dv_dt_visc.after_fab, to_host_fab(dv_dt_visc.before_fab), "dv_dt_visc");
    }
    if (du_dt_visc_gl90.allocated) {
        expect_arrays_equal(du_dt_visc_gl90.after_fab, to_host_fab(du_dt_visc_gl90.before_fab), "du_dt_visc_gl90");
    }
    if (dv_dt_visc_gl90.allocated) {
        expect_arrays_equal(dv_dt_visc_gl90.after_fab, to_host_fab(dv_dt_visc_gl90.before_fab), "dv_dt_visc_gl90");
    }
    if (du_dt_str.allocated) {
        expect_arrays_equal(du_dt_str.after_fab, to_host_fab(du_dt_str.before_fab), "du_dt_str");
    }
    if (dv_dt_str.allocated) {
        expect_arrays_equal(dv_dt_str.after_fab, to_host_fab(dv_dt_str.before_fab), "dv_dt_str");
    }
    if (taux_bot.allocated) {
        expect_arrays_equal(taux_bot.after_fab, to_host_fab(taux_bot.before_fab), "taux_bot");
    }
    if (tauy_bot.allocated) {
        expect_arrays_equal(tauy_bot.after_fab, to_host_fab(tauy_bot.before_fab), "tauy_bot");
    }
    if (KE_u.allocated) {
        expect_arrays_equal(KE_u.after_fab, to_host_fab(KE_u.before_fab), "KE_u");
    }
    if (KE_v.allocated) {
        expect_arrays_equal(KE_v.after_fab, to_host_fab(KE_v.before_fab), "KE_v");
    }
    if (u_old.allocated) {
        expect_arrays_equal(u_old.after_fab, to_host_fab(u_old.before_fab), "u_old");
    }
    if (v_old.allocated) {
        expect_arrays_equal(v_old.after_fab, to_host_fab(v_old.before_fab), "v_old");
    }
    if (vel_report_u.allocated) {
        expect_arrays_equal(vel_report_u.after_fab, to_host_fab(vel_report_u.before_fab), "vel_report_u");
    }
    if (vel_report_v.allocated) {
        expect_arrays_equal(vel_report_v.after_fab, to_host_fab(vel_report_v.before_fab), "vel_report_v");
    }
    if (dowrite_u.allocated) {
        expect_int_arrays_equal(dowrite_u.after_fab, to_host_int_fab(dowrite_u.before_fab), "dowrite_u");
    }
    if (dowrite_v.allocated) {
        expect_int_arrays_equal(dowrite_v.after_fab, to_host_int_fab(dowrite_v.before_fab), "dowrite_v");
    }
}
