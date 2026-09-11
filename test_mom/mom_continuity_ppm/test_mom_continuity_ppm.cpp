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
