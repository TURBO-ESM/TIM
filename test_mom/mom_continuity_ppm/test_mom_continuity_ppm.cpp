// Unit tests for MOM::ppm_limit_pos / PPM_reconstruction_x / PPM_reconstruction_y.
//
// Each test loads a captured Fortran (input, expected-output) pair from
// <data-dir>/<name>.{bin,meta}, runs the C++ kernel over equivalent AMReX
// containers, and compares the result against the captured "after" arrays.

// SKILLS: 0.3.1

#include <gtest/gtest.h>

#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>

#include "amrex_assertions.hpp"
#include "captured_io.hpp"
#include "data_dir.hpp"
#include "mom_continuity_ppm.hpp"

using test_mom::expect_arrays_equal;
using test_mom::to_host_fab;

namespace {

// Binds a "<field>_before" / "<field>_after" in/out array pair that the
// Fortran shim only captures when the corresponding container was
// associated at capture time (may be null-encoded or missing entirely --
// CapturedFile::is_associated() checks both). When absent, `arr` is left
// default-constructed (Array4<Real>{}, a null pointer), matching the
// kernel's own "may be absent (.p == nullptr)" parameter convention;
// `present` gates whether the post-call assertion runs at all.
struct OptionalInOutArray {
    amrex::FArrayBox before_fab;
    amrex::FArrayBox after_fab;
    amrex::Array4<amrex::Real> arr{};
    bool present = false;
};

OptionalInOutArray bind_optional_inout(const test_mom::CapturedFile& captured,
                                       const std::string& field) {
    OptionalInOutArray o;
    o.present = captured.is_associated("_" + field + "_before");
    if (o.present) {
        o.before_fab = captured.fab_device("_" + field + "_before");
        o.arr = o.before_fab.array();
        o.after_fab = captured.fab_host("_" + field + "_after");
    }
    return o;
}

} // namespace

// -------------------------------------------------------------------------
// ppm_limit_pos
// -------------------------------------------------------------------------
TEST(PpmLimitPos, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "ppm_limit_pos");

    const auto   bx        = captured.box("_bx");
    const auto   h_in      = captured.fab_device("_h_in");
    auto         h_L       = captured.fab_device("_h_L_before");
    auto         h_R       = captured.fab_device("_h_R_before");
    const auto   h_L_after = captured.fab_host("_h_L_after");
    const auto   h_R_after = captured.fab_host("_h_R_after");
    const double h_min     = captured.real64("_h_min");

    MOM::ppm_limit_pos(bx,
                       h_in.const_array(),
                       h_L.array(),
                       h_R.array(),
                       h_min);
    amrex::Gpu::synchronize();

    expect_arrays_equal(h_L_after, to_host_fab(h_L), "h_L");
    expect_arrays_equal(h_R_after, to_host_fab(h_R), "h_R");
}

// -------------------------------------------------------------------------
// PPM_reconstruction_x
// -------------------------------------------------------------------------
TEST(PpmReconstructionX, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "ppm_reconstruction_x");

    const auto   bxH        = captured.box("_bxH");
    const auto   h_in       = captured.fab_device("_h_in");
    auto         h_W        = captured.fab_device("_h_W_before");
    auto         h_E        = captured.fab_device("_h_E_before");
    const auto   mask2d     = captured.fab_device("_mask2d_t");
    const auto   h_W_after  = captured.fab_host("_h_W_after");
    const auto   h_E_after  = captured.fab_host("_h_E_after");
    const double h_min      = captured.real64("_h_min");
    const bool   monotonic  = captured.logical("_monotonic");
    const bool   simple_2nd = captured.logical("_simple_2nd");

    MOM::PPM_reconstruction_x(bxH,
                              h_in.const_array(),
                              h_W.array(),
                              h_E.array(),
                              mask2d.const_array(),
                              h_min,
                              monotonic,
                              simple_2nd,
                              /*OBC=*/nullptr);
    amrex::Gpu::synchronize();

    expect_arrays_equal(h_W_after, to_host_fab(h_W), "h_W");
    expect_arrays_equal(h_E_after, to_host_fab(h_E), "h_E");
}

// -------------------------------------------------------------------------
// PPM_reconstruction_y
// -------------------------------------------------------------------------
TEST(PpmReconstructionY, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "ppm_reconstruction_y");

    const auto   bxH        = captured.box("_bxH");
    const auto   h_in       = captured.fab_device("_h_in");
    auto         h_S        = captured.fab_device("_h_S_before");
    auto         h_N        = captured.fab_device("_h_N_before");
    const auto   mask2d     = captured.fab_device("_mask2d_t");
    const auto   h_S_after  = captured.fab_host("_h_S_after");
    const auto   h_N_after  = captured.fab_host("_h_N_after");
    const double h_min      = captured.real64("_h_min");
    const bool   monotonic  = captured.logical("_monotonic");
    const bool   simple_2nd = captured.logical("_simple_2nd");

    MOM::PPM_reconstruction_y(bxH,
                              h_in.const_array(),
                              h_S.array(),
                              h_N.array(),
                              mask2d.const_array(),
                              h_min,
                              monotonic,
                              simple_2nd,
                              /*OBC=*/nullptr);
    amrex::Gpu::synchronize();

    expect_arrays_equal(h_S_after, to_host_fab(h_S), "h_S");
    expect_arrays_equal(h_N_after, to_host_fab(h_N), "h_N");
}

// -------------------------------------------------------------------------
// ppm_limit_cw84 -- no capture available yet
// -------------------------------------------------------------------------
TEST(PpmLimitCw84, MatchesFortranCapture) {
    GTEST_SKIP() << "no captured ppm_limit_cw84.{bin,meta} fixture yet";
}

// -------------------------------------------------------------------------
// continuity_zonal_convergence
// -------------------------------------------------------------------------
//
// No capture/continuity_zonal_convergence.{bin,meta} fixture exists yet, so
// this test's field mapping is grounded directly in Fortran source (not a
// .meta file): submodules/MOM6/src/core/MOM_continuity_PPM.F90:1504-1578,
// the continuity_zonal_convergence shim's TIMH_capture arm (rec%add(...)
// calls at lines 1543-1549 and 1555). Cross-checked against the bind(C)
// interface (lines 158-172) -- both agree.
//
// MOM::continuity_zonal_convergence(bxC, h, uh, dt, IareaT, hin, h_min)
// advances h in place by the convergence of the zonal thickness flux. hin
// is recorded unconditionally by the shim (unlike visc_rem_u/v in the
// flux-thickness kernels), but the Fortran source itself allows hin to be
// unassociated -- if the captured call site had it that way, this fixture's
// _hin field will be empty and CapturedFile::fab_device("_hin") will throw
// rather than silently misbehave.
TEST(ContinuityZonalConvergence, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "continuity_zonal_convergence");

    const auto   bxC     = captured.box("_bxC");
    auto         h       = captured.fab_device("_h_before");
    const auto   h_after = captured.fab_host("_h_after");
    const auto   uh      = captured.fab_device("_uh");
    const double dt      = captured.real64("_dt");
    const auto   IareaT  = captured.fab_device("_IareaT");
    const auto   hin     = captured.fab_device("_hin");
    const double h_min   = captured.real64("_h_min");

    MOM::continuity_zonal_convergence(bxC,
                                      h.array(),
                                      uh.const_array(),
                                      dt,
                                      IareaT.const_array(),
                                      hin.const_array(),
                                      h_min);
    amrex::Gpu::synchronize();

    expect_arrays_equal(h_after, to_host_fab(h), "h");
}

// -------------------------------------------------------------------------
// continuity_meridional_convergence
// -------------------------------------------------------------------------
//
// No capture/continuity_meridional_convergence.{bin,meta} fixture exists
// yet, so this test's field mapping is grounded directly in Fortran source
// (not a .meta file): submodules/MOM6/src/core/MOM_continuity_PPM.F90:1626-1700,
// the continuity_meridional_convergence shim's TIMH_capture arm (rec%add(...)
// calls at lines 1665-1671 and 1677). Cross-checked against the bind(C)
// interface (lines 176-193) -- both agree.
//
// MOM::continuity_meridional_convergence(bxC, h, vh, dt, IareaT, hin, h_min)
// advances h in place by the convergence of the meridional thickness flux.
// hin is recorded unconditionally by the shim (unlike visc_rem_u/v in the
// flux-thickness kernels), but the Fortran source itself allows hin to be
// unassociated -- if the captured call site had it that way, this fixture's
// _hin field will be empty and CapturedFile::fab_device("_hin") will throw
// rather than silently misbehave.
TEST(ContinuityMeridionalConvergence, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "continuity_meridional_convergence");

    const auto   bxC     = captured.box("_bxC");
    auto         h       = captured.fab_device("_h_before");
    const auto   h_after = captured.fab_host("_h_after");
    const auto   vh      = captured.fab_device("_vh");
    const double dt      = captured.real64("_dt");
    const auto   IareaT  = captured.fab_device("_IareaT");
    const auto   hin     = captured.fab_device("_hin");
    const double h_min   = captured.real64("_h_min");

    MOM::continuity_meridional_convergence(bxC,
                                           h.array(),
                                           vh.const_array(),
                                           dt,
                                           IareaT.const_array(),
                                           hin.const_array(),
                                           h_min);
    amrex::Gpu::synchronize();

    expect_arrays_equal(h_after, to_host_fab(h), "h");
}
// -------------------------------------------------------------------------
// zonal_flux_thickness
// -------------------------------------------------------------------------
//
// No capture/zonal_flux_thickness.{bin,meta} fixture exists yet, so this
// test's field mapping is grounded directly in Fortran source (not a .meta
// file): submodules/MOM6/src/core/MOM_continuity_PPM.F90:3451-3568, the
// zonal_flux_thickness shim's TIMH_capture arm (rec%add(...) calls at
// lines 3512-3525 and 3533). Cross-checked against the bind(C) interface
// (lines 338-375) -- both agree.
TEST(ZonalFluxThickness, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "zonal_flux_thickness");

    const auto   bxC             = captured.box("_bxC");
    const auto   u               = captured.fab_device("_u");
    const auto   h               = captured.fab_device("_h");
    const auto   h_W             = captured.fab_device("_h_W");
    const auto   h_E             = captured.fab_device("_h_E");
    auto         h_u             = captured.fab_device("_h_u_before");
    const auto   h_u_after       = captured.fab_host("_h_u_after");
    const double dt              = captured.real64("_dt");
    const auto   dy_Cu           = captured.fab_device("_dy_Cu");
    const auto   IareaT          = captured.fab_device("_IareaT");
    const auto   IdxT            = captured.fab_device("_IdxT");
    const bool   vol_CFL         = captured.logical("_vol_CFL");
    const bool   marginal        = captured.logical("_marginal");
    const auto   por_face_areaU  = captured.fab_device("_por_face_areaU");
    // _visc_rem_u is captured only when associated -- may be absent
    // entirely from this fixture.
    amrex::FArrayBox visc_rem_u_fab;
    amrex::Array4<const amrex::Real> visc_rem_u{};
    if (captured.is_associated("_visc_rem_u")) {
        visc_rem_u_fab = captured.fab_device("_visc_rem_u");
        visc_rem_u = visc_rem_u_fab.const_array();
    }

    MOM::zonal_flux_thickness(bxC,
                              u.const_array(),
                              h.const_array(),
                              h_W.const_array(),
                              h_E.const_array(),
                              h_u.array(),
                              dt,
                              dy_Cu.const_array(),
                              IareaT.const_array(),
                              IdxT.const_array(),
                              vol_CFL,
                              marginal,
                              /*OBC=*/nullptr,
                              por_face_areaU.const_array(),
                              visc_rem_u);
    amrex::Gpu::synchronize();

    expect_arrays_equal(h_u_after, to_host_fab(h_u), "h_u");
}

// -------------------------------------------------------------------------
// meridional_flux_thickness
// -------------------------------------------------------------------------
//
// No capture/meridional_flux_thickness.{bin,meta} fixture exists yet, so
// this test's field mapping is grounded directly in Fortran source (not a
// .meta file): submodules/MOM6/src/core/MOM_continuity_PPM.F90:5292-5409,
// the meridional_flux_thickness shim's TIMH_capture arm (rec%add(...) calls
// at lines 5352-5365 and 5373). Cross-checked against the bind(C) interface
// (lines 380-417) -- both agree.
TEST(MeridionalFluxThickness, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "meridional_flux_thickness");

    const auto   bxC             = captured.box("_bxC");
    const auto   v               = captured.fab_device("_v");
    const auto   h               = captured.fab_device("_h");
    const auto   h_S             = captured.fab_device("_h_S");
    const auto   h_N             = captured.fab_device("_h_N");
    auto         h_v             = captured.fab_device("_h_v_before");
    const auto   h_v_after       = captured.fab_host("_h_v_after");
    const double dt              = captured.real64("_dt");
    const auto   dx_Cv           = captured.fab_device("_dx_Cv");
    const auto   IareaT          = captured.fab_device("_IareaT");
    const auto   IdyT            = captured.fab_device("_IdyT");
    const bool   vol_CFL         = captured.logical("_vol_CFL");
    const bool   marginal        = captured.logical("_marginal");
    const auto   por_face_areaV  = captured.fab_device("_por_face_areaV");
    // _visc_rem_v is captured only when associated -- may be absent
    // entirely from this fixture.
    amrex::FArrayBox visc_rem_v_fab;
    amrex::Array4<const amrex::Real> visc_rem_v{};
    if (captured.is_associated("_visc_rem_v")) {
        visc_rem_v_fab = captured.fab_device("_visc_rem_v");
        visc_rem_v = visc_rem_v_fab.const_array();
    }

    MOM::meridional_flux_thickness(bxC,
                                   v.const_array(),
                                   h.const_array(),
                                   h_S.const_array(),
                                   h_N.const_array(),
                                   h_v.array(),
                                   dt,
                                   dx_Cv.const_array(),
                                   IareaT.const_array(),
                                   IdyT.const_array(),
                                   vol_CFL,
                                   marginal,
                                   /*OBC=*/nullptr,
                                   por_face_areaV.const_array(),
                                   visc_rem_v);
    amrex::Gpu::synchronize();

    expect_arrays_equal(h_v_after, to_host_fab(h_v), "h_v");
}
// -------------------------------------------------------------------------
// set_merid_BT_cont
// -------------------------------------------------------------------------
//
// No capture/set_merid_bt_cont.{bin,meta} fixture exists yet, so this
// test's field mapping is grounded directly in Fortran source (not a .meta
// file): submodules/MOM6/src/core/MOM_continuity_PPM.F90:5993-6153, the
// set_merid_BT_cont shim's TIMH_capture arm (rec%add(...) calls at lines
// 6066-6091 and 6100-6105). Cross-checked against the bind(C) interface
// (lines 585-638) -- both agree.
//
// MOM::set_merid_BT_cont computes the effective open face areas and
// barotropic-velocity corrections at meridional faces as a function of
// barotropic flow, for use by the barotropic solver's transport-adjustment
// iteration. do_I is captured as a LogicalArray_t, read here via
// int_fab_device() into an amrex::IArrayBox. Only CS.vol_CFL is captured --
// the kernel reads no other transport_adjust_CS_C field, so the rest of CS
// is left default-initialized.
TEST(SetMeridBtCont, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "set_merid_bt_cont");

    const auto   bxC             = captured.box("_bxC");
    const auto   v                = captured.fab_device("_v");
    const auto   h_in             = captured.fab_device("_h_in");
    const auto   h_S              = captured.fab_device("_h_S");
    const auto   h_N              = captured.fab_device("_h_N");
    auto         FA_v_S0          = captured.fab_device("_FA_v_S0_before");
    auto         FA_v_N0          = captured.fab_device("_FA_v_N0_before");
    auto         FA_v_SS          = captured.fab_device("_FA_v_SS_before");
    auto         FA_v_NN          = captured.fab_device("_FA_v_NN_before");
    auto         vBT_SS           = captured.fab_device("_vBT_SS_before");
    auto         vBT_NN           = captured.fab_device("_vBT_NN_before");
    const auto   FA_v_S0_after    = captured.fab_host("_FA_v_S0_after");
    const auto   FA_v_N0_after    = captured.fab_host("_FA_v_N0_after");
    const auto   FA_v_SS_after    = captured.fab_host("_FA_v_SS_after");
    const auto   FA_v_NN_after    = captured.fab_host("_FA_v_NN_after");
    const auto   vBT_SS_after     = captured.fab_host("_vBT_SS_after");
    const auto   vBT_NN_after     = captured.fab_host("_vBT_NN_after");
    const auto   dv0              = captured.fab_device("_dv0");
    const auto   vh_tot_0         = captured.fab_device("_vh_tot_0");
    const auto   dvhdv_tot_0      = captured.fab_device("_dvhdv_tot_0");
    const auto   dv_max_CFL       = captured.fab_device("_dv_max_CFL");
    const auto   dv_min_CFL       = captured.fab_device("_dv_min_CFL");
    const double dt               = captured.real64("_dt");
    const auto   dyCv             = captured.fab_device("_dyCv");
    const auto   dx_Cv            = captured.fab_device("_dx_Cv");
    const auto   IareaT           = captured.fab_device("_IareaT");
    const auto   IdyT             = captured.fab_device("_IdyT");
    transport_adjust_CS_C CS{};
    CS.vol_CFL                    = captured.logical("_vol_CFL");
    const auto   visc_rem         = captured.fab_device("_visc_rem");
    const auto   visc_rem_max     = captured.fab_device("_visc_rem_max");
    const auto   do_I             = captured.int_fab_device("_do_I");
    const auto   por_face_areaV   = captured.fab_device("_por_face_areaV");

    MOM::set_merid_BT_cont(bxC,
                           v.const_array(),
                           h_in.const_array(),
                           h_S.const_array(),
                           h_N.const_array(),
                           FA_v_S0.array(),
                           FA_v_N0.array(),
                           FA_v_SS.array(),
                           FA_v_NN.array(),
                           vBT_SS.array(),
                           vBT_NN.array(),
                           dv0.const_array(),
                           vh_tot_0.const_array(),
                           dvhdv_tot_0.const_array(),
                           dv_max_CFL.const_array(),
                           dv_min_CFL.const_array(),
                           dt,
                           dyCv.const_array(),
                           dx_Cv.const_array(),
                           IareaT.const_array(),
                           IdyT.const_array(),
                           CS,
                           visc_rem.const_array(),
                           visc_rem_max.const_array(),
                           do_I.const_array(),
                           por_face_areaV.const_array());
    amrex::Gpu::synchronize();

    expect_arrays_equal(FA_v_S0_after, to_host_fab(FA_v_S0), "FA_v_S0");
    expect_arrays_equal(FA_v_N0_after, to_host_fab(FA_v_N0), "FA_v_N0");
    expect_arrays_equal(FA_v_SS_after, to_host_fab(FA_v_SS), "FA_v_SS");
    expect_arrays_equal(FA_v_NN_after, to_host_fab(FA_v_NN), "FA_v_NN");
    expect_arrays_equal(vBT_SS_after,  to_host_fab(vBT_SS),  "vBT_SS");
    expect_arrays_equal(vBT_NN_after,  to_host_fab(vBT_NN),  "vBT_NN");
}

// -------------------------------------------------------------------------
// set_zonal_BT_cont
// -------------------------------------------------------------------------
//
// No capture/set_zonal_bt_cont.{bin,meta} fixture exists yet, so this
// test's field mapping is grounded directly in Fortran source (not a .meta
// file): submodules/MOM6/src/core/MOM_continuity_PPM.F90:4152-4311, the
// set_zonal_BT_cont shim's TIMH_capture arm (rec%add(...) calls at lines
// 4224-4249 and 4258-4263). Cross-checked against the bind(C) interface
// (lines 527-580) -- both agree.
//
// MOM::set_zonal_BT_cont computes the effective open face areas and
// barotropic-velocity corrections at zonal faces as a function of
// barotropic flow, for use by the barotropic solver's transport-adjustment
// iteration. do_I is captured as a LogicalArray_t, read here via
// int_fab_device() into an amrex::IArrayBox. Only CS.vol_CFL is captured --
// the kernel reads no other transport_adjust_CS_C field, so the rest of CS
// is left default-initialized.
TEST(SetZonalBtCont, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "set_zonal_bt_cont");

    const auto   bxC             = captured.box("_bxC");
    const auto   u                = captured.fab_device("_u");
    const auto   h_in             = captured.fab_device("_h_in");
    const auto   h_W              = captured.fab_device("_h_W");
    const auto   h_E              = captured.fab_device("_h_E");
    auto         FA_u_W0          = captured.fab_device("_FA_u_W0_before");
    auto         FA_u_E0          = captured.fab_device("_FA_u_E0_before");
    auto         FA_u_WW          = captured.fab_device("_FA_u_WW_before");
    auto         FA_u_EE          = captured.fab_device("_FA_u_EE_before");
    auto         uBT_WW           = captured.fab_device("_uBT_WW_before");
    auto         uBT_EE           = captured.fab_device("_uBT_EE_before");
    const auto   FA_u_W0_after    = captured.fab_host("_FA_u_W0_after");
    const auto   FA_u_E0_after    = captured.fab_host("_FA_u_E0_after");
    const auto   FA_u_WW_after    = captured.fab_host("_FA_u_WW_after");
    const auto   FA_u_EE_after    = captured.fab_host("_FA_u_EE_after");
    const auto   uBT_WW_after     = captured.fab_host("_uBT_WW_after");
    const auto   uBT_EE_after     = captured.fab_host("_uBT_EE_after");
    const auto   du0              = captured.fab_device("_du0");
    const auto   uh_tot_0         = captured.fab_device("_uh_tot_0");
    const auto   duhdu_tot_0      = captured.fab_device("_duhdu_tot_0");
    const auto   du_max_CFL       = captured.fab_device("_du_max_CFL");
    const auto   du_min_CFL       = captured.fab_device("_du_min_CFL");
    const double dt               = captured.real64("_dt");
    const auto   dxCu             = captured.fab_device("_dxCu");
    const auto   dy_Cu            = captured.fab_device("_dy_Cu");
    const auto   IareaT           = captured.fab_device("_IareaT");
    const auto   IdxT             = captured.fab_device("_IdxT");
    transport_adjust_CS_C CS{};
    CS.vol_CFL                    = captured.logical("_vol_CFL");
    const auto   visc_rem         = captured.fab_device("_visc_rem");
    const auto   visc_rem_max     = captured.fab_device("_visc_rem_max");
    const auto   do_I             = captured.int_fab_device("_do_I");
    const auto   por_face_areaU   = captured.fab_device("_por_face_areaU");

    MOM::set_zonal_BT_cont(bxC,
                           u.const_array(),
                           h_in.const_array(),
                           h_W.const_array(),
                           h_E.const_array(),
                           FA_u_W0.array(),
                           FA_u_E0.array(),
                           FA_u_WW.array(),
                           FA_u_EE.array(),
                           uBT_WW.array(),
                           uBT_EE.array(),
                           du0.const_array(),
                           uh_tot_0.const_array(),
                           duhdu_tot_0.const_array(),
                           du_max_CFL.const_array(),
                           du_min_CFL.const_array(),
                           dt,
                           dxCu.const_array(),
                           dy_Cu.const_array(),
                           IareaT.const_array(),
                           IdxT.const_array(),
                           CS,
                           visc_rem.const_array(),
                           visc_rem_max.const_array(),
                           do_I.const_array(),
                           por_face_areaU.const_array());
    amrex::Gpu::synchronize();

    expect_arrays_equal(FA_u_W0_after, to_host_fab(FA_u_W0), "FA_u_W0");
    expect_arrays_equal(FA_u_E0_after, to_host_fab(FA_u_E0), "FA_u_E0");
    expect_arrays_equal(FA_u_WW_after, to_host_fab(FA_u_WW), "FA_u_WW");
    expect_arrays_equal(FA_u_EE_after, to_host_fab(FA_u_EE), "FA_u_EE");
    expect_arrays_equal(uBT_WW_after,  to_host_fab(uBT_WW),  "uBT_WW");
    expect_arrays_equal(uBT_EE_after,  to_host_fab(uBT_EE),  "uBT_EE");
}

// -------------------------------------------------------------------------
// zonal_flux_adjust
// -------------------------------------------------------------------------
//
// No capture/zonal_flux_adjust.{bin,meta} fixture exists yet, so this
// test's field mapping is grounded directly in Fortran source (not a .meta
// file): submodules/MOM6/src/core/MOM_continuity_PPM.F90:3799-3954, the
// zonal_flux_adjust shim's TIMH_capture arm (rec%add(...) calls at lines
// 3868-3894 and 3898-3900). Cross-checked against the bind(C) interface
// (lines 422-472) -- both agree. _uhbt/_uh_3d_before are captured only when
// unassociated; is_associated()/bind_optional_inout() handle both cases.
//
// MOM::zonal_flux_adjust Newton-iterates a barotropic velocity correction
// per zonal face so that the vertically-summed zonal mass/volume transport
// matches the target barotropic transport. do_I_in is captured as a
// LogicalArray_t, read here via int_fab_device() into an amrex::IArrayBox.
// obc is passed as nullptr -- OceanOBC is not yet implemented on the C++
// side; the kernel aborts if given a non-null obc.
TEST(ZonalFluxAdjust, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "zonal_flux_adjust");

    const auto   bxC             = captured.box("_bxC");
    const auto   u                = captured.fab_device("_u");
    const auto   h_in             = captured.fab_device("_h_in");
    const auto   h_W              = captured.fab_device("_h_W");
    const auto   h_E              = captured.fab_device("_h_E");
    const auto   uh_tot_0         = captured.fab_device("_uh_tot_0");
    const auto   duhdu_tot_0      = captured.fab_device("_duhdu_tot_0");
    auto         du               = captured.fab_device("_du_before");
    const auto   du_after         = captured.fab_host("_du_after");
    const auto   du_max_CFL       = captured.fab_device("_du_max_CFL");
    const auto   du_min_CFL       = captured.fab_device("_du_min_CFL");
    const double dt               = captured.real64("_dt");
    const auto   dy_Cu            = captured.fab_device("_dy_Cu");
    const auto   IareaT           = captured.fab_device("_IareaT");
    const auto   IdxT             = captured.fab_device("_IdxT");
    transport_adjust_CS_C CS{};
    CS.tol_eta                    = captured.real64("_tol_eta");
    CS.tol_vel                    = captured.real64("_tol_vel");
    CS.better_iter                = captured.logical("_better_iter");
    CS.vol_CFL                    = captured.logical("_vol_CFL");
    const auto   visc_rem         = captured.fab_device("_visc_rem");
    const auto   do_I_in          = captured.int_fab_device("_do_I_in");
    const auto   por_face_areaU   = captured.fab_device("_por_face_areaU");
    // _uhbt is captured unconditionally but may still be null-encoded
    // (Fortran container unassociated at capture time); _uh_3d is an
    // optional in/out output captured only when associated.
    amrex::FArrayBox uhbt_fab;
    amrex::Array4<const amrex::Real> uhbt{};
    if (captured.is_associated("_uhbt")) {
        uhbt_fab = captured.fab_device("_uhbt");
        uhbt = uhbt_fab.const_array();
    }
    auto uh_3d = bind_optional_inout(captured, "uh_3d");

    MOM::zonal_flux_adjust(bxC,
                           u.const_array(),
                           h_in.const_array(),
                           h_W.const_array(),
                           h_E.const_array(),
                           uh_tot_0.const_array(),
                           duhdu_tot_0.const_array(),
                           du.array(),
                           du_max_CFL.const_array(),
                           du_min_CFL.const_array(),
                           dt,
                           dy_Cu.const_array(),
                           IareaT.const_array(),
                           IdxT.const_array(),
                           CS,
                           visc_rem.const_array(),
                           do_I_in.const_array(),
                           por_face_areaU.const_array(),
                           uhbt,
                           uh_3d.arr,
                           /*obc=*/nullptr);
    amrex::Gpu::synchronize();

    expect_arrays_equal(du_after, to_host_fab(du), "du");
    if (uh_3d.present) expect_arrays_equal(uh_3d.after_fab, to_host_fab(uh_3d.before_fab), "uh_3d");
}

// -------------------------------------------------------------------------
// meridional_flux_adjust
// -------------------------------------------------------------------------
//
// No capture/meridional_flux_adjust.{bin,meta} fixture exists yet, so this
// test's field mapping is grounded directly in Fortran source (not a .meta
// file): submodules/MOM6/src/core/MOM_continuity_PPM.F90:5638-5791, the
// meridional_flux_adjust shim's TIMH_capture arm (rec%add(...) calls at
// lines 5712-5738 and 5742-5744). Cross-checked against the bind(C)
// interface (lines 477-522) -- both agree. _vhbt/_vh_3d_before are captured
// only when unassociated; is_associated()/bind_optional_inout() handle both
// cases.
//
// MOM::meridional_flux_adjust Newton-iterates a barotropic velocity
// correction per meridional face so that the vertically-summed meridional
// mass/volume transport matches the target barotropic transport. do_I_in
// is captured as a LogicalArray_t, read here via int_fab_device() into an
// amrex::IArrayBox. obc is passed as nullptr -- OceanOBC is not yet
// implemented on the C++ side; the kernel aborts if given a non-null obc.
TEST(MeridionalFluxAdjust, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "meridional_flux_adjust");

    const auto   bxC             = captured.box("_bxC");
    const auto   v                = captured.fab_device("_v");
    const auto   h_in             = captured.fab_device("_h_in");
    const auto   h_S              = captured.fab_device("_h_S");
    const auto   h_N              = captured.fab_device("_h_N");
    const auto   vh_tot_0         = captured.fab_device("_vh_tot_0");
    const auto   dvhdv_tot_0      = captured.fab_device("_dvhdv_tot_0");
    auto         dv               = captured.fab_device("_dv_before");
    const auto   dv_after         = captured.fab_host("_dv_after");
    const auto   dv_max_CFL       = captured.fab_device("_dv_max_CFL");
    const auto   dv_min_CFL       = captured.fab_device("_dv_min_CFL");
    const double dt               = captured.real64("_dt");
    const auto   dx_Cv            = captured.fab_device("_dx_Cv");
    const auto   IareaT           = captured.fab_device("_IareaT");
    const auto   IdyT             = captured.fab_device("_IdyT");
    transport_adjust_CS_C CS{};
    CS.tol_eta                    = captured.real64("_tol_eta");
    CS.tol_vel                    = captured.real64("_tol_vel");
    CS.better_iter                = captured.logical("_better_iter");
    CS.vol_CFL                    = captured.logical("_vol_CFL");
    const auto   visc_rem         = captured.fab_device("_visc_rem");
    const auto   do_I_in          = captured.int_fab_device("_do_I_in");
    const auto   por_face_areaV   = captured.fab_device("_por_face_areaV");
    // _vhbt is captured unconditionally but may still be null-encoded
    // (Fortran container unassociated at capture time); _vh_3d is an
    // optional in/out output captured only when associated.
    amrex::FArrayBox vhbt_fab;
    amrex::Array4<const amrex::Real> vhbt{};
    if (captured.is_associated("_vhbt")) {
        vhbt_fab = captured.fab_device("_vhbt");
        vhbt = vhbt_fab.const_array();
    }
    auto vh_3d = bind_optional_inout(captured, "vh_3d");

    MOM::meridional_flux_adjust(bxC,
                                v.const_array(),
                                h_in.const_array(),
                                h_S.const_array(),
                                h_N.const_array(),
                                vh_tot_0.const_array(),
                                dvhdv_tot_0.const_array(),
                                dv.array(),
                                dv_max_CFL.const_array(),
                                dv_min_CFL.const_array(),
                                dt,
                                dx_Cv.const_array(),
                                IareaT.const_array(),
                                IdyT.const_array(),
                                CS,
                                visc_rem.const_array(),
                                do_I_in.const_array(),
                                por_face_areaV.const_array(),
                                vhbt,
                                vh_3d.arr,
                                /*obc=*/nullptr);
    amrex::Gpu::synchronize();

    expect_arrays_equal(dv_after, to_host_fab(dv), "dv");
    if (vh_3d.present) expect_arrays_equal(vh_3d.after_fab, to_host_fab(vh_3d.before_fab), "vh_3d");
}
}

// -------------------------------------------------------------------------
// zonal_BT_mass_flux
// -------------------------------------------------------------------------
//
// No capture/zonal_bt_mass_flux.{bin,meta} fixture exists yet, so this
// test's field mapping is grounded directly in Fortran source (not a .meta
// file): submodules/MOM6/src/core/MOM_continuity_PPM.F90:3095-3205, the
// zonal_BT_mass_flux shim's TIMH_capture arm (rec%add(...) calls at lines
// 3149-3160 and 3168). Cross-checked against the bind(C) interface
// (lines 270-299) -- both agree.
//
// MOM::zonal_BT_mass_flux accumulates the vertically-summed zonal
// barotropic mass/volume transport across the water column, for use as the
// barotropic solver's target transport. Only CS.vol_CFL is captured -- the
// kernel reads no other transport_adjust_CS_C field, so the rest of CS is
// left default-initialized. obc is passed as nullptr -- OceanOBC is not
// yet implemented on the C++ side; the kernel aborts if given a non-null
// obc.
TEST(ZonalBtMassFlux, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "zonal_bt_mass_flux");

    const auto   bxC             = captured.box("_bxC");
    const auto   u                = captured.fab_device("_u");
    const auto   h_in             = captured.fab_device("_h_in");
    const auto   h_W              = captured.fab_device("_h_W");
    const auto   h_E              = captured.fab_device("_h_E");
    auto         uhbt             = captured.fab_device("_uhbt_before");
    const auto   uhbt_after       = captured.fab_host("_uhbt_after");
    const double dt               = captured.real64("_dt");
    const auto   dy_Cu            = captured.fab_device("_dy_Cu");
    const auto   IareaT           = captured.fab_device("_IareaT");
    const auto   IdxT             = captured.fab_device("_IdxT");
    transport_adjust_CS_C CS{};
    CS.vol_CFL                    = captured.logical("_vol_CFL");
    const auto   por_face_areaU   = captured.fab_device("_por_face_areaU");

    MOM::zonal_BT_mass_flux(bxC,
                            u.const_array(),
                            h_in.const_array(),
                            h_W.const_array(),
                            h_E.const_array(),
                            uhbt.array(),
                            dt,
                            dy_Cu.const_array(),
                            IareaT.const_array(),
                            IdxT.const_array(),
                            CS,
                            /*obc=*/nullptr,
                            por_face_areaU.const_array());
    amrex::Gpu::synchronize();

    expect_arrays_equal(uhbt_after, to_host_fab(uhbt), "uhbt");
}

// -------------------------------------------------------------------------
// meridional_BT_mass_flux
// -------------------------------------------------------------------------
//
// No capture/meridional_bt_mass_flux.{bin,meta} fixture exists yet, so this
// test's field mapping is grounded directly in Fortran source (not a .meta
// file): submodules/MOM6/src/core/MOM_continuity_PPM.F90:5030-5143, the
// meridional_BT_mass_flux shim's TIMH_capture arm (rec%add(...) calls at
// lines 5087-5098 and 5106). Cross-checked against the bind(C) interface
// (lines 304-333) -- both agree.
//
// MOM::meridional_BT_mass_flux accumulates the vertically-summed
// meridional barotropic mass/volume transport across the water column, for
// use as the barotropic solver's target transport. Only CS.vol_CFL is
// captured -- the kernel reads no other transport_adjust_CS_C field, so
// the rest of CS is left default-initialized. obc is passed as nullptr --
// OceanOBC is not yet implemented on the C++ side; the kernel aborts if
// given a non-null obc.
TEST(MeridionalBtMassFlux, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "meridional_bt_mass_flux");

    const auto   bxC             = captured.box("_bxC");
    const auto   v                = captured.fab_device("_v");
    const auto   h_in             = captured.fab_device("_h_in");
    const auto   h_S              = captured.fab_device("_h_S");
    const auto   h_N              = captured.fab_device("_h_N");
    auto         vhbt             = captured.fab_device("_vhbt_before");
    const auto   vhbt_after       = captured.fab_host("_vhbt_after");
    const double dt               = captured.real64("_dt");
    const auto   dx_Cv            = captured.fab_device("_dx_Cv");
    const auto   IareaT           = captured.fab_device("_IareaT");
    const auto   IdyT             = captured.fab_device("_IdyT");
    transport_adjust_CS_C CS{};
    CS.vol_CFL                    = captured.logical("_vol_CFL");
    const auto   por_face_areaV   = captured.fab_device("_por_face_areaV");

    MOM::meridional_BT_mass_flux(bxC,
                                 v.const_array(),
                                 h_in.const_array(),
                                 h_S.const_array(),
                                 h_N.const_array(),
                                 vhbt.array(),
                                 dt,
                                 dx_Cv.const_array(),
                                 IareaT.const_array(),
                                 IdyT.const_array(),
                                 CS,
                                 /*obc=*/nullptr,
                                 por_face_areaV.const_array());
    amrex::Gpu::synchronize();

    expect_arrays_equal(vhbt_after, to_host_fab(vhbt), "vhbt");
}
