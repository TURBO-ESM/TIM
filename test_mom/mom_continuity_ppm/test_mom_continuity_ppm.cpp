// Unit tests for MOM::ppm_limit_pos / PPM_reconstruction_x / PPM_reconstruction_y.
//
// Each test loads a captured Fortran (input, expected-output) pair from
// <data-dir>/<name>.{bin,meta}, runs the C++ kernel over equivalent AMReX
// containers, and compares the result against the captured "after" arrays.

#include <gtest/gtest.h>

#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>

#include "amrex_assertions.hpp"
#include "captured_io.hpp"
#include "data_dir.hpp"
#include "mom_continuity_ppm.hpp"

using test_mom::expect_arrays_equal;

// -------------------------------------------------------------------------
// ppm_limit_pos
// -------------------------------------------------------------------------
TEST(PpmLimitPos, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "ppm_limit_pos");

    const auto   bx        = captured.box("_bx");
    const auto   h_in      = captured.fab("_h_in");
    auto         h_L       = captured.fab("_h_L_before");
    auto         h_R       = captured.fab("_h_R_before");
    const auto   h_L_after = captured.fab("_h_L_after");
    const auto   h_R_after = captured.fab("_h_R_after");
    const double h_min     = captured.real64("_h_min");

    MOM::ppm_limit_pos(bx,
                       h_in.const_array(),
                       h_L.array(),
                       h_R.array(),
                       h_min);
    amrex::Gpu::synchronize();

    expect_arrays_equal(h_L_after, h_L, "h_L");
    expect_arrays_equal(h_R_after, h_R, "h_R");
}

// -------------------------------------------------------------------------
// PPM_reconstruction_x
// -------------------------------------------------------------------------
TEST(PpmReconstructionX, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "ppm_reconstruction_x");

    const auto   bxH        = captured.box("_bxH");
    const auto   h_in       = captured.fab("_h_in");
    auto         h_W        = captured.fab("_h_W_before");
    auto         h_E        = captured.fab("_h_E_before");
    const auto   mask2d     = captured.fab("_mask2d_t");
    const auto   h_W_after  = captured.fab("_h_W_after");
    const auto   h_E_after  = captured.fab("_h_E_after");
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

    expect_arrays_equal(h_W_after, h_W, "h_W");
    expect_arrays_equal(h_E_after, h_E, "h_E");
}

// -------------------------------------------------------------------------
// PPM_reconstruction_y
// -------------------------------------------------------------------------
TEST(PpmReconstructionY, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "ppm_reconstruction_y");

    const auto   bxH        = captured.box("_bxH");
    const auto   h_in       = captured.fab("_h_in");
    auto         h_S        = captured.fab("_h_S_before");
    auto         h_N        = captured.fab("_h_N_before");
    const auto   mask2d     = captured.fab("_mask2d_t");
    const auto   h_S_after  = captured.fab("_h_S_after");
    const auto   h_N_after  = captured.fab("_h_N_after");
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

    expect_arrays_equal(h_S_after, h_S, "h_S");
    expect_arrays_equal(h_N_after, h_N, "h_N");
}

// -------------------------------------------------------------------------
// ppm_limit_cw84 -- no capture available yet
// -------------------------------------------------------------------------
TEST(PpmLimitCw84, MatchesFortranCapture) {
    GTEST_SKIP() << "no captured ppm_limit_cw84.{bin,meta} fixture yet";
}
