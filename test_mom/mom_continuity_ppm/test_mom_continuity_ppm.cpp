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
