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
