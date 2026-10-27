// Unit tests for MOM::smooth_x9_h / smooth_x9_uv / smooth_gme_h_pass / smooth_gme_q_pass /
// hor_visc_backscatter_h / hor_visc_backscatter_q / hor_visc_gme_bt / hor_visc_leith_grad_a /
// hor_visc_leith_grad_b / hor_visc_leithy_m / hor_visc_leithy_ah / hor_visc_leithy_ah_h /
// hor_visc_strain / hor_visc_h_coeffs_a / hor_visc_h_coeffs_b / hor_visc_stress_div /
// hor_visc_q_coeffs / horizontal_viscosity_tr.
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
#include "mom_hor_visc.hpp"

using test_mom::expect_arrays_equal;
using test_mom::to_host_fab;

// -------------------------------------------------------------------------
// Field mapping: MOM6 src/parameterizations/lateral/MOM_hor_visc.F90:10079-10096 (smooth_x9_h capture block).
TEST(SmoothX9H, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "smooth_x9_h");

    const auto bxH0          = captured.box("_bxH0");
    const auto mask2dT       = captured.fab_device("_mask2dT");
    auto       field_h       = captured.fab_device("_field_h_before");
    const auto field_h_after = captured.fab_host("_field_h_after");
    const bool zero_land     = captured.logical("_zero_land");

    MOM::smooth_x9_h(bxH0,
                     mask2dT.const_array(),
                     field_h.array(),
                     zero_land);
    amrex::Gpu::synchronize();

    expect_arrays_equal(field_h_after, to_host_fab(field_h), "field_h");
}

// -------------------------------------------------------------------------
// Field mapping: MOM6 src/parameterizations/lateral/MOM_hor_visc.F90:10228-10249 (smooth_x9_uv capture block).
TEST(SmoothX9Uv, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "smooth_x9_uv");

    const auto bxU0          = captured.box("_bxU0");
    const auto bxV0          = captured.box("_bxV0");
    const auto mask2dCu      = captured.fab_device("_mask2dCu");
    const auto mask2dCv      = captured.fab_device("_mask2dCv");
    auto       field_u       = captured.fab_device("_field_u_before");
    auto       field_v       = captured.fab_device("_field_v_before");
    const auto field_u_after = captured.fab_host("_field_u_after");
    const auto field_v_after = captured.fab_host("_field_v_after");
    const bool zero_land     = captured.logical("_zero_land");

    MOM::smooth_x9_uv(bxU0, bxV0,
                      mask2dCu.const_array(),
                      mask2dCv.const_array(),
                      field_u.array(),
                      field_v.array(),
                      zero_land);
    amrex::Gpu::synchronize();

    expect_arrays_equal(field_u_after, to_host_fab(field_u), "field_u");
    expect_arrays_equal(field_v_after, to_host_fab(field_v), "field_v");
}

// -------------------------------------------------------------------------
// Field mapping: MOM6 src/parameterizations/lateral/MOM_hor_visc.F90:9788-9804 (smooth_GME_h_pass capture block).
TEST(SmoothGmeHPass, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "smooth_gme_h_pass");

    const auto bxPass           = captured.box("_bxPass");
    const auto mask2dT          = captured.fab_device("_mask2dT");
    auto       GME_flux_h       = captured.fab_device("_GME_flux_h_before");
    const auto GME_flux_h_after = captured.fab_host("_GME_flux_h_after");

    MOM::smooth_gme_h_pass(bxPass,
                           mask2dT.const_array(),
                           GME_flux_h.array());
    amrex::Gpu::synchronize();

    expect_arrays_equal(GME_flux_h_after, to_host_fab(GME_flux_h), "GME_flux_h");
}

// -------------------------------------------------------------------------
// Field mapping: MOM6 src/parameterizations/lateral/MOM_hor_visc.F90:9887-9903 (smooth_GME_q_pass capture block).
TEST(SmoothGmeQPass, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "smooth_gme_q_pass");

    const auto bxPass           = captured.box("_bxPass");
    const auto mask2dBu         = captured.fab_device("_mask2dBu");
    auto       GME_flux_q       = captured.fab_device("_GME_flux_q_before");
    const auto GME_flux_q_after = captured.fab_host("_GME_flux_q_after");

    MOM::smooth_gme_q_pass(bxPass,
                           mask2dBu.const_array(),
                           GME_flux_q.array());
    amrex::Gpu::synchronize();

    expect_arrays_equal(GME_flux_q_after, to_host_fab(GME_flux_q), "GME_flux_q");
}

// -------------------------------------------------------------------------
// Field mapping: MOM6 src/parameterizations/lateral/MOM_hor_visc.F90:7358-7386 (hor_visc_backscatter_h capture block).
TEST(HorViscBackscatterH, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "hor_visc_backscatter_h");

    const auto bxHstr = captured.box("_bxHstr");

    hv_opts_CS_C CS{};
    CS.Laplacian              = captured.logical("_CS_Laplacian");
    CS.biharmonic             = captured.logical("_CS_biharmonic");
    CS.debug                  = captured.logical("_CS_debug");
    CS.no_slip                = captured.logical("_CS_no_slip");
    CS.bound_Kh               = captured.logical("_CS_bound_Kh");
    CS.EY24_EBT_BS            = captured.logical("_CS_EY24_EBT_BS");
    CS.bound_Ah               = captured.logical("_CS_bound_Ah");
    CS.Re_Ah                  = captured.real64("_CS_Re_Ah");
    CS.KS_coef                = captured.real64("_CS_KS_coef");
    CS.backscatter_underbound = captured.logical("_CS_backscatter_underbound");
    CS.Smagorinsky_Kh         = captured.logical("_CS_Smagorinsky_Kh");
    CS.Smagorinsky_Ah         = captured.logical("_CS_Smagorinsky_Ah");
    CS.Leith_Kh               = captured.logical("_CS_Leith_Kh");
    CS.Modified_Leith         = captured.logical("_CS_Modified_Leith");
    CS.use_beta_in_Leith      = captured.logical("_CS_use_beta_in_Leith");
    CS.Leith_Ah               = captured.logical("_CS_Leith_Ah");
    CS.use_Leithy             = captured.logical("_CS_use_Leithy");
    CS.c_K                    = captured.real64("_CS_c_K");
    CS.smooth_Ah              = captured.logical("_CS_smooth_Ah");
    CS.use_QG_Leith_visc      = captured.logical("_CS_use_QG_Leith_visc");
    CS.bound_Coriolis         = captured.logical("_CS_bound_Coriolis");
    CS.Kh_bg_min              = captured.real64("_CS_Kh_bg_min");
    CS.FrictWork_bug          = captured.logical("_CS_FrictWork_bug");
    CS.OBC_strain_bug         = captured.logical("_CS_OBC_strain_bug");
    CS.use_land_mask          = captured.logical("_CS_use_land_mask");
    CS.anisotropic            = captured.logical("_CS_anisotropic");
    CS.add_LES_viscosity      = captured.logical("_CS_add_LES_viscosity");
    CS.Kh_aniso               = captured.real64("_CS_Kh_aniso");
    CS.res_scale_MEKE         = captured.logical("_CS_res_scale_MEKE");
    CS.use_GME                = captured.logical("_CS_use_GME");
    CS.answer_date            = captured.integer("_CS_answer_date");
    CS.GME_h0                 = captured.real64("_CS_GME_h0");
    CS.GME_efficiency         = captured.real64("_CS_GME_efficiency");
    CS.GME_limiter            = captured.real64("_CS_GME_limiter");
    CS.min_grid_Kh            = captured.real64("_CS_min_grid_Kh");
    CS.min_grid_Ah            = captured.real64("_CS_min_grid_Ah");
    CS.use_cont_thick         = captured.logical("_CS_use_cont_thick");
    CS.use_ZB2020             = captured.logical("_CS_use_ZB2020");
    CS.use_circulation        = captured.logical("_CS_use_circulation");
    CS.do_grid_Re_Ah          = captured.logical("_CS_do_grid_Re_Ah");
    CS.do_grid_Re_Kh          = captured.logical("_CS_do_grid_Re_Kh");
    CS.do_Ah_h                = captured.logical("_CS_do_Ah_h");
    CS.do_Ah_q                = captured.logical("_CS_do_Ah_q");
    CS.do_Kh_h                = captured.logical("_CS_do_Kh_h");
    CS.do_Kh_q                = captured.logical("_CS_do_Kh_q");
    CS.do_GME_coeff_h         = captured.logical("_CS_do_GME_coeff_h");
    CS.do_GME_coeff_q         = captured.logical("_CS_do_GME_coeff_q");
    CS.do_vort_xy_q           = captured.logical("_CS_do_vort_xy_q");
    CS.do_div_xx_h            = captured.logical("_CS_do_div_xx_h");
    CS.do_sh_xy_q             = captured.logical("_CS_do_sh_xy_q");
    CS.do_sh_xx_h             = captured.logical("_CS_do_sh_xx_h");
    CS.do_FrictWork           = captured.logical("_CS_do_FrictWork");
    CS.do_FrictWorkIntz       = captured.logical("_CS_do_FrictWorkIntz");
    CS.do_FrictWork_bh        = captured.logical("_CS_do_FrictWork_bh");
    CS.do_FrictWorkIntz_bh    = captured.logical("_CS_do_FrictWorkIntz_bh");
    CS.do_normstress          = captured.logical("_CS_do_normstress");
    CS.do_shearstress         = captured.logical("_CS_do_shearstress");
    CS.do_BS_coeff_h          = captured.logical("_CS_do_BS_coeff_h");
    CS.do_BS_coeff_q          = captured.logical("_CS_do_BS_coeff_q");

    // The MEKE-derived lateral viscosity coefficient [L2 T-1 ~> m2 s-1]; may be unallocated
    amrex::FArrayBox MEKE_Ku_fab;
    amrex::Array4<const amrex::Real> MEKE_Ku{};
    if (captured.logical("_MEKE_Ku_allocated")) {
        MEKE_Ku_fab = captured.fab_device("_MEKE_Ku_before");
        MEKE_Ku     = MEKE_Ku_fab.const_array();
    }

    // Vertical structure function used in backscatter [nondim]; may be unallocated
    amrex::FArrayBox BS_struct_fab;
    amrex::Array4<const amrex::Real> BS_struct{};
    if (captured.logical("_VarMix_BS_struct_allocated")) {
        BS_struct_fab = captured.fab_device("_VarMix_BS_struct_before");
        BS_struct     = BS_struct_fab.const_array();
    }

    const bool use_kh_struct     = captured.logical("_use_kh_struct");
    const int  kstart            = captured.integer("_kstart");
    const auto visc_limit_h_flag = captured.fab_device("_visc_limit_h_flag");
    const auto sh_xx             = captured.fab_device("_sh_xx");
    auto       str_xx            = captured.fab_device("_str_xx_before");
    auto       BS_coeff_h        = captured.fab_device("_BS_coeff_h_before");
    const auto str_xx_after      = captured.fab_host("_str_xx_after");
    const auto BS_coeff_h_after  = captured.fab_host("_BS_coeff_h_after");

    MOM::hor_visc_backscatter_h(bxHstr, CS,
                                MEKE_Ku, BS_struct,
                                use_kh_struct, kstart,
                                visc_limit_h_flag.const_array(),
                                sh_xx.const_array(),
                                str_xx.array(),
                                BS_coeff_h.array());
    amrex::Gpu::synchronize();

    expect_arrays_equal(str_xx_after,     to_host_fab(str_xx),     "str_xx");
    expect_arrays_equal(BS_coeff_h_after, to_host_fab(BS_coeff_h), "BS_coeff_h");
}

// -------------------------------------------------------------------------
// Field mapping: MOM6 src/parameterizations/lateral/MOM_hor_visc.F90:7536-7564 (hor_visc_backscatter_q capture block).
TEST(HorViscBackscatterQ, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "hor_visc_backscatter_q");

    const auto bxQstr = captured.box("_bxQstr");

    hv_opts_CS_C CS{};
    CS.Laplacian              = captured.logical("_CS_Laplacian");
    CS.biharmonic             = captured.logical("_CS_biharmonic");
    CS.debug                  = captured.logical("_CS_debug");
    CS.no_slip                = captured.logical("_CS_no_slip");
    CS.bound_Kh               = captured.logical("_CS_bound_Kh");
    CS.EY24_EBT_BS            = captured.logical("_CS_EY24_EBT_BS");
    CS.bound_Ah               = captured.logical("_CS_bound_Ah");
    CS.Re_Ah                  = captured.real64("_CS_Re_Ah");
    CS.KS_coef                = captured.real64("_CS_KS_coef");
    CS.backscatter_underbound = captured.logical("_CS_backscatter_underbound");
    CS.Smagorinsky_Kh         = captured.logical("_CS_Smagorinsky_Kh");
    CS.Smagorinsky_Ah         = captured.logical("_CS_Smagorinsky_Ah");
    CS.Leith_Kh               = captured.logical("_CS_Leith_Kh");
    CS.Modified_Leith         = captured.logical("_CS_Modified_Leith");
    CS.use_beta_in_Leith      = captured.logical("_CS_use_beta_in_Leith");
    CS.Leith_Ah               = captured.logical("_CS_Leith_Ah");
    CS.use_Leithy             = captured.logical("_CS_use_Leithy");
    CS.c_K                    = captured.real64("_CS_c_K");
    CS.smooth_Ah              = captured.logical("_CS_smooth_Ah");
    CS.use_QG_Leith_visc      = captured.logical("_CS_use_QG_Leith_visc");
    CS.bound_Coriolis         = captured.logical("_CS_bound_Coriolis");
    CS.Kh_bg_min              = captured.real64("_CS_Kh_bg_min");
    CS.FrictWork_bug          = captured.logical("_CS_FrictWork_bug");
    CS.OBC_strain_bug         = captured.logical("_CS_OBC_strain_bug");
    CS.use_land_mask          = captured.logical("_CS_use_land_mask");
    CS.anisotropic            = captured.logical("_CS_anisotropic");
    CS.add_LES_viscosity      = captured.logical("_CS_add_LES_viscosity");
    CS.Kh_aniso               = captured.real64("_CS_Kh_aniso");
    CS.res_scale_MEKE         = captured.logical("_CS_res_scale_MEKE");
    CS.use_GME                = captured.logical("_CS_use_GME");
    CS.answer_date            = captured.integer("_CS_answer_date");
    CS.GME_h0                 = captured.real64("_CS_GME_h0");
    CS.GME_efficiency         = captured.real64("_CS_GME_efficiency");
    CS.GME_limiter            = captured.real64("_CS_GME_limiter");
    CS.min_grid_Kh            = captured.real64("_CS_min_grid_Kh");
    CS.min_grid_Ah            = captured.real64("_CS_min_grid_Ah");
    CS.use_cont_thick         = captured.logical("_CS_use_cont_thick");
    CS.use_ZB2020             = captured.logical("_CS_use_ZB2020");
    CS.use_circulation        = captured.logical("_CS_use_circulation");
    CS.do_grid_Re_Ah          = captured.logical("_CS_do_grid_Re_Ah");
    CS.do_grid_Re_Kh          = captured.logical("_CS_do_grid_Re_Kh");
    CS.do_Ah_h                = captured.logical("_CS_do_Ah_h");
    CS.do_Ah_q                = captured.logical("_CS_do_Ah_q");
    CS.do_Kh_h                = captured.logical("_CS_do_Kh_h");
    CS.do_Kh_q                = captured.logical("_CS_do_Kh_q");
    CS.do_GME_coeff_h         = captured.logical("_CS_do_GME_coeff_h");
    CS.do_GME_coeff_q         = captured.logical("_CS_do_GME_coeff_q");
    CS.do_vort_xy_q           = captured.logical("_CS_do_vort_xy_q");
    CS.do_div_xx_h            = captured.logical("_CS_do_div_xx_h");
    CS.do_sh_xy_q             = captured.logical("_CS_do_sh_xy_q");
    CS.do_sh_xx_h             = captured.logical("_CS_do_sh_xx_h");
    CS.do_FrictWork           = captured.logical("_CS_do_FrictWork");
    CS.do_FrictWorkIntz       = captured.logical("_CS_do_FrictWorkIntz");
    CS.do_FrictWork_bh        = captured.logical("_CS_do_FrictWork_bh");
    CS.do_FrictWorkIntz_bh    = captured.logical("_CS_do_FrictWorkIntz_bh");
    CS.do_normstress          = captured.logical("_CS_do_normstress");
    CS.do_shearstress         = captured.logical("_CS_do_shearstress");
    CS.do_BS_coeff_h          = captured.logical("_CS_do_BS_coeff_h");
    CS.do_BS_coeff_q          = captured.logical("_CS_do_BS_coeff_q");

    // The MEKE-derived lateral viscosity coefficient [L2 T-1 ~> m2 s-1]; may be unallocated
    amrex::FArrayBox MEKE_Ku_fab;
    amrex::Array4<const amrex::Real> MEKE_Ku{};
    if (captured.logical("_MEKE_Ku_allocated")) {
        MEKE_Ku_fab = captured.fab_device("_MEKE_Ku_before");
        MEKE_Ku     = MEKE_Ku_fab.const_array();
    }

    // Vertical structure function used in backscatter [nondim]; may be unallocated
    amrex::FArrayBox BS_struct_fab;
    amrex::Array4<const amrex::Real> BS_struct{};
    if (captured.logical("_VarMix_BS_struct_allocated")) {
        BS_struct_fab = captured.fab_device("_VarMix_BS_struct_before");
        BS_struct     = BS_struct_fab.const_array();
    }

    const bool use_kh_struct     = captured.logical("_use_kh_struct");
    const int  kstart            = captured.integer("_kstart");
    const auto visc_limit_q_flag = captured.fab_device("_visc_limit_q_flag");
    const auto sh_xy             = captured.fab_device("_sh_xy");
    auto       str_xy            = captured.fab_device("_str_xy_before");
    auto       BS_coeff_q        = captured.fab_device("_BS_coeff_q_before");
    const auto str_xy_after      = captured.fab_host("_str_xy_after");
    const auto BS_coeff_q_after  = captured.fab_host("_BS_coeff_q_after");

    MOM::hor_visc_backscatter_q(bxQstr, CS,
                                MEKE_Ku, BS_struct,
                                use_kh_struct, kstart,
                                visc_limit_q_flag.const_array(),
                                sh_xy.const_array(),
                                str_xy.array(),
                                BS_coeff_q.array());
    amrex::Gpu::synchronize();

    expect_arrays_equal(str_xy_after,     to_host_fab(str_xy),     "str_xy");
    expect_arrays_equal(BS_coeff_q_after, to_host_fab(BS_coeff_q), "BS_coeff_q");
}

// -------------------------------------------------------------------------
// Field mapping: MOM6 src/parameterizations/lateral/MOM_hor_visc.F90:8426-8472 (hor_visc_GME_bt capture block).
TEST(HorViscGmeBt, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "hor_visc_gme_bt");

    const auto bxH0 = captured.box("_bxH0");
    const auto bxQ0 = captured.box("_bxQ0");

    hv_opts_CS_C CS{};
    CS.Laplacian              = captured.logical("_CS_Laplacian");
    CS.biharmonic             = captured.logical("_CS_biharmonic");
    CS.debug                  = captured.logical("_CS_debug");
    CS.no_slip                = captured.logical("_CS_no_slip");
    CS.bound_Kh               = captured.logical("_CS_bound_Kh");
    CS.EY24_EBT_BS            = captured.logical("_CS_EY24_EBT_BS");
    CS.bound_Ah               = captured.logical("_CS_bound_Ah");
    CS.Re_Ah                  = captured.real64("_CS_Re_Ah");
    CS.KS_coef                = captured.real64("_CS_KS_coef");
    CS.backscatter_underbound = captured.logical("_CS_backscatter_underbound");
    CS.Smagorinsky_Kh         = captured.logical("_CS_Smagorinsky_Kh");
    CS.Smagorinsky_Ah         = captured.logical("_CS_Smagorinsky_Ah");
    CS.Leith_Kh               = captured.logical("_CS_Leith_Kh");
    CS.Modified_Leith         = captured.logical("_CS_Modified_Leith");
    CS.use_beta_in_Leith      = captured.logical("_CS_use_beta_in_Leith");
    CS.Leith_Ah               = captured.logical("_CS_Leith_Ah");
    CS.use_Leithy             = captured.logical("_CS_use_Leithy");
    CS.c_K                    = captured.real64("_CS_c_K");
    CS.smooth_Ah              = captured.logical("_CS_smooth_Ah");
    CS.use_QG_Leith_visc      = captured.logical("_CS_use_QG_Leith_visc");
    CS.bound_Coriolis         = captured.logical("_CS_bound_Coriolis");
    CS.Kh_bg_min              = captured.real64("_CS_Kh_bg_min");
    CS.FrictWork_bug          = captured.logical("_CS_FrictWork_bug");
    CS.OBC_strain_bug         = captured.logical("_CS_OBC_strain_bug");
    CS.use_land_mask          = captured.logical("_CS_use_land_mask");
    CS.anisotropic            = captured.logical("_CS_anisotropic");
    CS.add_LES_viscosity      = captured.logical("_CS_add_LES_viscosity");
    CS.Kh_aniso               = captured.real64("_CS_Kh_aniso");
    CS.res_scale_MEKE         = captured.logical("_CS_res_scale_MEKE");
    CS.use_GME                = captured.logical("_CS_use_GME");
    CS.answer_date            = captured.integer("_CS_answer_date");
    CS.GME_h0                 = captured.real64("_CS_GME_h0");
    CS.GME_efficiency         = captured.real64("_CS_GME_efficiency");
    CS.GME_limiter            = captured.real64("_CS_GME_limiter");
    CS.min_grid_Kh            = captured.real64("_CS_min_grid_Kh");
    CS.min_grid_Ah            = captured.real64("_CS_min_grid_Ah");
    CS.use_cont_thick         = captured.logical("_CS_use_cont_thick");
    CS.use_ZB2020             = captured.logical("_CS_use_ZB2020");
    CS.use_circulation        = captured.logical("_CS_use_circulation");
    CS.do_grid_Re_Ah          = captured.logical("_CS_do_grid_Re_Ah");
    CS.do_grid_Re_Kh          = captured.logical("_CS_do_grid_Re_Kh");
    CS.do_Ah_h                = captured.logical("_CS_do_Ah_h");
    CS.do_Ah_q                = captured.logical("_CS_do_Ah_q");
    CS.do_Kh_h                = captured.logical("_CS_do_Kh_h");
    CS.do_Kh_q                = captured.logical("_CS_do_Kh_q");
    CS.do_GME_coeff_h         = captured.logical("_CS_do_GME_coeff_h");
    CS.do_GME_coeff_q         = captured.logical("_CS_do_GME_coeff_q");
    CS.do_vort_xy_q           = captured.logical("_CS_do_vort_xy_q");
    CS.do_div_xx_h            = captured.logical("_CS_do_div_xx_h");
    CS.do_sh_xy_q             = captured.logical("_CS_do_sh_xy_q");
    CS.do_sh_xx_h             = captured.logical("_CS_do_sh_xx_h");
    CS.do_FrictWork           = captured.logical("_CS_do_FrictWork");
    CS.do_FrictWorkIntz       = captured.logical("_CS_do_FrictWorkIntz");
    CS.do_FrictWork_bh        = captured.logical("_CS_do_FrictWork_bh");
    CS.do_FrictWorkIntz_bh    = captured.logical("_CS_do_FrictWorkIntz_bh");
    CS.do_normstress          = captured.logical("_CS_do_normstress");
    CS.do_shearstress         = captured.logical("_CS_do_shearstress");
    CS.do_BS_coeff_h          = captured.logical("_CS_do_BS_coeff_h");
    CS.do_BS_coeff_q          = captured.logical("_CS_do_BS_coeff_q");

    const auto IdyCu    = captured.fab_device("_IdyCu");
    const auto IdxCv    = captured.fab_device("_IdxCv");
    const auto IdyCv    = captured.fab_device("_IdyCv");
    const auto IdxCu    = captured.fab_device("_IdxCu");
    const auto mask2dBu = captured.fab_device("_mask2dBu");
    const auto mask2dCu = captured.fab_device("_mask2dCu");
    const auto mask2dCv = captured.fab_device("_mask2dCv");
    const auto mask2dT  = captured.fab_device("_mask2dT");

    // Pre-calculated dx/dy at h points [nondim]; may be unallocated
    amrex::FArrayBox dx_dyT_fab;
    amrex::Array4<const amrex::Real> dx_dyT{};
    if (captured.logical("_metrics_dx_dyT_allocated")) {
        dx_dyT_fab = captured.fab_device("_metrics_dx_dyT_before");
        dx_dyT = dx_dyT_fab.const_array();
    }
    // Pre-calculated dy/dx at h points [nondim]; may be unallocated
    amrex::FArrayBox dy_dxT_fab;
    amrex::Array4<const amrex::Real> dy_dxT{};
    if (captured.logical("_metrics_dy_dxT_allocated")) {
        dy_dxT_fab = captured.fab_device("_metrics_dy_dxT_before");
        dy_dxT = dy_dxT_fab.const_array();
    }
    // Pre-calculated dx/dy at q points [nondim]; may be unallocated
    amrex::FArrayBox dx_dyBu_fab;
    amrex::Array4<const amrex::Real> dx_dyBu{};
    if (captured.logical("_metrics_dx_dyBu_allocated")) {
        dx_dyBu_fab = captured.fab_device("_metrics_dx_dyBu_before");
        dx_dyBu = dx_dyBu_fab.const_array();
    }
    // Pre-calculated dy/dx at q points [nondim]; may be unallocated
    amrex::FArrayBox dy_dxBu_fab;
    amrex::Array4<const amrex::Real> dy_dxBu{};
    if (captured.logical("_metrics_dy_dxBu_allocated")) {
        dy_dxBu_fab = captured.fab_device("_metrics_dy_dxBu_before");
        dy_dxBu = dy_dxBu_fab.const_array();
    }

    const auto h     = captured.fab_device("_h_before");
    const auto ubtav = captured.fab_device("_ubtav");
    const auto vbtav = captured.fab_device("_vbtav");
    auto       dudx_bt           = captured.fab_device("_dudx_bt_before");
    auto       dvdy_bt           = captured.fab_device("_dvdy_bt_before");
    auto       dvdx_bt           = captured.fab_device("_dvdx_bt_before");
    auto       dudy_bt           = captured.fab_device("_dudy_bt_before");
    auto       sh_xx_bt          = captured.fab_device("_sh_xx_bt_before");
    auto       sh_xy_bt          = captured.fab_device("_sh_xy_bt_before");
    auto       GME_effic_h       = captured.fab_device("_GME_effic_h_before");
    auto       GME_effic_q       = captured.fab_device("_GME_effic_q_before");
    const auto dudx_bt_after     = captured.fab_host("_dudx_bt_after");
    const auto dvdy_bt_after     = captured.fab_host("_dvdy_bt_after");
    const auto dvdx_bt_after     = captured.fab_host("_dvdx_bt_after");
    const auto dudy_bt_after     = captured.fab_host("_dudy_bt_after");
    const auto sh_xx_bt_after    = captured.fab_host("_sh_xx_bt_after");
    const auto sh_xy_bt_after    = captured.fab_host("_sh_xy_bt_after");
    const auto GME_effic_h_after = captured.fab_host("_GME_effic_h_after");
    const auto GME_effic_q_after = captured.fab_host("_GME_effic_q_after");

    MOM::hor_visc_gme_bt(bxH0, bxQ0, CS,
                         IdyCu.const_array(), IdxCv.const_array(), IdyCv.const_array(), IdxCu.const_array(),
                         mask2dBu.const_array(), mask2dCu.const_array(), mask2dCv.const_array(), mask2dT.const_array(),
                         dx_dyT, dy_dxT, dx_dyBu, dy_dxBu,
                         h.const_array(), ubtav.const_array(), vbtav.const_array(),
                         dudx_bt.array(), dvdy_bt.array(), dvdx_bt.array(), dudy_bt.array(),
                         sh_xx_bt.array(), sh_xy_bt.array(), GME_effic_h.array(), GME_effic_q.array());
    amrex::Gpu::synchronize();

    expect_arrays_equal(dudx_bt_after,     to_host_fab(dudx_bt),     "dudx_bt");
    expect_arrays_equal(dvdy_bt_after,     to_host_fab(dvdy_bt),     "dvdy_bt");
    expect_arrays_equal(dvdx_bt_after,     to_host_fab(dvdx_bt),     "dvdx_bt");
    expect_arrays_equal(dudy_bt_after,     to_host_fab(dudy_bt),     "dudy_bt");
    expect_arrays_equal(sh_xx_bt_after,    to_host_fab(sh_xx_bt),    "sh_xx_bt");
    expect_arrays_equal(sh_xy_bt_after,    to_host_fab(sh_xy_bt),    "sh_xy_bt");
    expect_arrays_equal(GME_effic_h_after, to_host_fab(GME_effic_h), "GME_effic_h");
    expect_arrays_equal(GME_effic_q_after, to_host_fab(GME_effic_q), "GME_effic_q");
}

// -------------------------------------------------------------------------
// Field mapping: MOM6 src/parameterizations/lateral/MOM_hor_visc.F90:6918-6976 (hor_visc_Leith_grad_a capture block).
TEST(HorViscLeithGradA, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "hor_visc_leith_grad_a");

    const auto bxH        = captured.box("_bxH");
    const auto bxKh       = captured.box("_bxKh");
    const auto bxKh_i     = captured.box("_bxKh_i");
    const auto bxKh_j     = captured.box("_bxKh_j");
    const auto bxQstr     = captured.box("_bxQstr");
    const auto bxQvsmooth = captured.box("_bxQvsmooth");

    hv_opts_CS_C CS{};
    CS.Laplacian              = captured.logical("_CS_Laplacian");
    CS.biharmonic             = captured.logical("_CS_biharmonic");
    CS.debug                  = captured.logical("_CS_debug");
    CS.no_slip                = captured.logical("_CS_no_slip");
    CS.bound_Kh               = captured.logical("_CS_bound_Kh");
    CS.EY24_EBT_BS            = captured.logical("_CS_EY24_EBT_BS");
    CS.bound_Ah               = captured.logical("_CS_bound_Ah");
    CS.Re_Ah                  = captured.real64("_CS_Re_Ah");
    CS.KS_coef                = captured.real64("_CS_KS_coef");
    CS.backscatter_underbound = captured.logical("_CS_backscatter_underbound");
    CS.Smagorinsky_Kh         = captured.logical("_CS_Smagorinsky_Kh");
    CS.Smagorinsky_Ah         = captured.logical("_CS_Smagorinsky_Ah");
    CS.Leith_Kh               = captured.logical("_CS_Leith_Kh");
    CS.Modified_Leith         = captured.logical("_CS_Modified_Leith");
    CS.use_beta_in_Leith      = captured.logical("_CS_use_beta_in_Leith");
    CS.Leith_Ah               = captured.logical("_CS_Leith_Ah");
    CS.use_Leithy             = captured.logical("_CS_use_Leithy");
    CS.c_K                    = captured.real64("_CS_c_K");
    CS.smooth_Ah              = captured.logical("_CS_smooth_Ah");
    CS.use_QG_Leith_visc      = captured.logical("_CS_use_QG_Leith_visc");
    CS.bound_Coriolis         = captured.logical("_CS_bound_Coriolis");
    CS.Kh_bg_min              = captured.real64("_CS_Kh_bg_min");
    CS.FrictWork_bug          = captured.logical("_CS_FrictWork_bug");
    CS.OBC_strain_bug         = captured.logical("_CS_OBC_strain_bug");
    CS.use_land_mask          = captured.logical("_CS_use_land_mask");
    CS.anisotropic            = captured.logical("_CS_anisotropic");
    CS.add_LES_viscosity      = captured.logical("_CS_add_LES_viscosity");
    CS.Kh_aniso               = captured.real64("_CS_Kh_aniso");
    CS.res_scale_MEKE         = captured.logical("_CS_res_scale_MEKE");
    CS.use_GME                = captured.logical("_CS_use_GME");
    CS.answer_date            = captured.integer("_CS_answer_date");
    CS.GME_h0                 = captured.real64("_CS_GME_h0");
    CS.GME_efficiency         = captured.real64("_CS_GME_efficiency");
    CS.GME_limiter            = captured.real64("_CS_GME_limiter");
    CS.min_grid_Kh            = captured.real64("_CS_min_grid_Kh");
    CS.min_grid_Ah            = captured.real64("_CS_min_grid_Ah");
    CS.use_cont_thick         = captured.logical("_CS_use_cont_thick");
    CS.use_ZB2020             = captured.logical("_CS_use_ZB2020");
    CS.use_circulation        = captured.logical("_CS_use_circulation");
    CS.do_grid_Re_Ah          = captured.logical("_CS_do_grid_Re_Ah");
    CS.do_grid_Re_Kh          = captured.logical("_CS_do_grid_Re_Kh");
    CS.do_Ah_h                = captured.logical("_CS_do_Ah_h");
    CS.do_Ah_q                = captured.logical("_CS_do_Ah_q");
    CS.do_Kh_h                = captured.logical("_CS_do_Kh_h");
    CS.do_Kh_q                = captured.logical("_CS_do_Kh_q");
    CS.do_GME_coeff_h         = captured.logical("_CS_do_GME_coeff_h");
    CS.do_GME_coeff_q         = captured.logical("_CS_do_GME_coeff_q");
    CS.do_vort_xy_q           = captured.logical("_CS_do_vort_xy_q");
    CS.do_div_xx_h            = captured.logical("_CS_do_div_xx_h");
    CS.do_sh_xy_q             = captured.logical("_CS_do_sh_xy_q");
    CS.do_sh_xx_h             = captured.logical("_CS_do_sh_xx_h");
    CS.do_FrictWork           = captured.logical("_CS_do_FrictWork");
    CS.do_FrictWorkIntz       = captured.logical("_CS_do_FrictWorkIntz");
    CS.do_FrictWork_bh        = captured.logical("_CS_do_FrictWork_bh");
    CS.do_FrictWorkIntz_bh    = captured.logical("_CS_do_FrictWorkIntz_bh");
    CS.do_normstress          = captured.logical("_CS_do_normstress");
    CS.do_shearstress         = captured.logical("_CS_do_shearstress");
    CS.do_BS_coeff_h          = captured.logical("_CS_do_BS_coeff_h");
    CS.do_BS_coeff_q          = captured.logical("_CS_do_BS_coeff_q");

    const auto dyBu           = captured.fab_device("_dyBu");
    const auto IdxBu          = captured.fab_device("_IdxBu");
    const auto IdyCu          = captured.fab_device("_IdyCu");
    const auto dxBu           = captured.fab_device("_dxBu");
    const auto IdyBu          = captured.fab_device("_IdyBu");
    const auto IdxCv          = captured.fab_device("_IdxCv");
    const auto IdyCv          = captured.fab_device("_IdyCv");
    const auto IdxCu          = captured.fab_device("_IdxCu");
    const auto dF_dx          = captured.fab_device("_dF_dx");
    const auto dF_dy          = captured.fab_device("_dF_dy");
    const auto dudx           = captured.fab_device("_dudx");
    const auto dvdy           = captured.fab_device("_dvdy");
    const auto vort_xy        = captured.fab_device("_vort_xy");
    const auto vort_xy_smooth = captured.fab_device("_vort_xy_smooth");
    auto       grad_vort_mag_h_2d       = captured.fab_device("_grad_vort_mag_h_2d_before");
    auto       grad_div_mag_h           = captured.fab_device("_grad_div_mag_h_before");
    auto       grad_vort_mag_q_2d       = captured.fab_device("_grad_vort_mag_q_2d_before");
    auto       grad_div_mag_q           = captured.fab_device("_grad_div_mag_q_before");
    auto       Del2vort_q               = captured.fab_device("_Del2vort_q_before");
    auto       vort_xy_dx               = captured.fab_device("_vort_xy_dx_before");
    auto       vort_xy_dx_smooth        = captured.fab_device("_vort_xy_dx_smooth_before");
    auto       div_xx_dy                = captured.fab_device("_div_xx_dy_before");
    auto       vort_xy_dy               = captured.fab_device("_vort_xy_dy_before");
    auto       vort_xy_dy_smooth        = captured.fab_device("_vort_xy_dy_smooth_before");
    auto       div_xx_dx                = captured.fab_device("_div_xx_dx_before");
    const auto grad_vort_mag_h_2d_after = captured.fab_host("_grad_vort_mag_h_2d_after");
    const auto grad_div_mag_h_after     = captured.fab_host("_grad_div_mag_h_after");
    const auto grad_vort_mag_q_2d_after = captured.fab_host("_grad_vort_mag_q_2d_after");
    const auto grad_div_mag_q_after     = captured.fab_host("_grad_div_mag_q_after");
    const auto Del2vort_q_after         = captured.fab_host("_Del2vort_q_after");
    const auto vort_xy_dx_after         = captured.fab_host("_vort_xy_dx_after");
    const auto vort_xy_dx_smooth_after  = captured.fab_host("_vort_xy_dx_smooth_after");
    const auto div_xx_dy_after          = captured.fab_host("_div_xx_dy_after");
    const auto vort_xy_dy_after         = captured.fab_host("_vort_xy_dy_after");
    const auto vort_xy_dy_smooth_after  = captured.fab_host("_vort_xy_dy_smooth_after");
    const auto div_xx_dx_after          = captured.fab_host("_div_xx_dx_after");

    MOM::hor_visc_leith_grad_a(bxH, bxKh, bxKh_i, bxKh_j, bxQstr, bxQvsmooth, CS,
                               dyBu.const_array(), IdxBu.const_array(), IdyCu.const_array(), dxBu.const_array(),
                               IdyBu.const_array(),
                               IdxCv.const_array(), IdyCv.const_array(), IdxCu.const_array(), dF_dx.const_array(),
                               dF_dy.const_array(),
                               dudx.const_array(), dvdy.const_array(), vort_xy.const_array(),
                               vort_xy_smooth.const_array(),
                               grad_vort_mag_h_2d.array(), grad_div_mag_h.array(), grad_vort_mag_q_2d.array(),
                               grad_div_mag_q.array(),
                               Del2vort_q.array(), vort_xy_dx.array(), vort_xy_dx_smooth.array(), div_xx_dy.array(),
                               vort_xy_dy.array(), vort_xy_dy_smooth.array(), div_xx_dx.array());
    amrex::Gpu::synchronize();

    expect_arrays_equal(grad_vort_mag_h_2d_after, to_host_fab(grad_vort_mag_h_2d), "grad_vort_mag_h_2d");
    expect_arrays_equal(grad_div_mag_h_after,     to_host_fab(grad_div_mag_h),     "grad_div_mag_h");
    expect_arrays_equal(grad_vort_mag_q_2d_after, to_host_fab(grad_vort_mag_q_2d), "grad_vort_mag_q_2d");
    expect_arrays_equal(grad_div_mag_q_after,     to_host_fab(grad_div_mag_q),     "grad_div_mag_q");
    expect_arrays_equal(Del2vort_q_after,         to_host_fab(Del2vort_q),         "Del2vort_q");
    expect_arrays_equal(vort_xy_dx_after,         to_host_fab(vort_xy_dx),         "vort_xy_dx");
    expect_arrays_equal(vort_xy_dx_smooth_after,  to_host_fab(vort_xy_dx_smooth),  "vort_xy_dx_smooth");
    expect_arrays_equal(div_xx_dy_after,          to_host_fab(div_xx_dy),          "div_xx_dy");
    expect_arrays_equal(vort_xy_dy_after,         to_host_fab(vort_xy_dy),         "vort_xy_dy");
    expect_arrays_equal(vort_xy_dy_smooth_after,  to_host_fab(vort_xy_dy_smooth),  "vort_xy_dy_smooth");
    expect_arrays_equal(div_xx_dx_after,          to_host_fab(div_xx_dx),          "div_xx_dx");
}

// -------------------------------------------------------------------------
// Field mapping: MOM6 src/parameterizations/lateral/MOM_hor_visc.F90:7130-7162 (hor_visc_Leith_grad_b capture block).
TEST(HorViscLeithGradB, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "hor_visc_leith_grad_b");

    const auto bxKh   = captured.box("_bxKh");
    const auto bxQstr = captured.box("_bxQstr");

    hv_opts_CS_C CS{};
    CS.Laplacian              = captured.logical("_CS_Laplacian");
    CS.biharmonic             = captured.logical("_CS_biharmonic");
    CS.debug                  = captured.logical("_CS_debug");
    CS.no_slip                = captured.logical("_CS_no_slip");
    CS.bound_Kh               = captured.logical("_CS_bound_Kh");
    CS.EY24_EBT_BS            = captured.logical("_CS_EY24_EBT_BS");
    CS.bound_Ah               = captured.logical("_CS_bound_Ah");
    CS.Re_Ah                  = captured.real64("_CS_Re_Ah");
    CS.KS_coef                = captured.real64("_CS_KS_coef");
    CS.backscatter_underbound = captured.logical("_CS_backscatter_underbound");
    CS.Smagorinsky_Kh         = captured.logical("_CS_Smagorinsky_Kh");
    CS.Smagorinsky_Ah         = captured.logical("_CS_Smagorinsky_Ah");
    CS.Leith_Kh               = captured.logical("_CS_Leith_Kh");
    CS.Modified_Leith         = captured.logical("_CS_Modified_Leith");
    CS.use_beta_in_Leith      = captured.logical("_CS_use_beta_in_Leith");
    CS.Leith_Ah               = captured.logical("_CS_Leith_Ah");
    CS.use_Leithy             = captured.logical("_CS_use_Leithy");
    CS.c_K                    = captured.real64("_CS_c_K");
    CS.smooth_Ah              = captured.logical("_CS_smooth_Ah");
    CS.use_QG_Leith_visc      = captured.logical("_CS_use_QG_Leith_visc");
    CS.bound_Coriolis         = captured.logical("_CS_bound_Coriolis");
    CS.Kh_bg_min              = captured.real64("_CS_Kh_bg_min");
    CS.FrictWork_bug          = captured.logical("_CS_FrictWork_bug");
    CS.OBC_strain_bug         = captured.logical("_CS_OBC_strain_bug");
    CS.use_land_mask          = captured.logical("_CS_use_land_mask");
    CS.anisotropic            = captured.logical("_CS_anisotropic");
    CS.add_LES_viscosity      = captured.logical("_CS_add_LES_viscosity");
    CS.Kh_aniso               = captured.real64("_CS_Kh_aniso");
    CS.res_scale_MEKE         = captured.logical("_CS_res_scale_MEKE");
    CS.use_GME                = captured.logical("_CS_use_GME");
    CS.answer_date            = captured.integer("_CS_answer_date");
    CS.GME_h0                 = captured.real64("_CS_GME_h0");
    CS.GME_efficiency         = captured.real64("_CS_GME_efficiency");
    CS.GME_limiter            = captured.real64("_CS_GME_limiter");
    CS.min_grid_Kh            = captured.real64("_CS_min_grid_Kh");
    CS.min_grid_Ah            = captured.real64("_CS_min_grid_Ah");
    CS.use_cont_thick         = captured.logical("_CS_use_cont_thick");
    CS.use_ZB2020             = captured.logical("_CS_use_ZB2020");
    CS.use_circulation        = captured.logical("_CS_use_circulation");
    CS.do_grid_Re_Ah          = captured.logical("_CS_do_grid_Re_Ah");
    CS.do_grid_Re_Kh          = captured.logical("_CS_do_grid_Re_Kh");
    CS.do_Ah_h                = captured.logical("_CS_do_Ah_h");
    CS.do_Ah_q                = captured.logical("_CS_do_Ah_q");
    CS.do_Kh_h                = captured.logical("_CS_do_Kh_h");
    CS.do_Kh_q                = captured.logical("_CS_do_Kh_q");
    CS.do_GME_coeff_h         = captured.logical("_CS_do_GME_coeff_h");
    CS.do_GME_coeff_q         = captured.logical("_CS_do_GME_coeff_q");
    CS.do_vort_xy_q           = captured.logical("_CS_do_vort_xy_q");
    CS.do_div_xx_h            = captured.logical("_CS_do_div_xx_h");
    CS.do_sh_xy_q             = captured.logical("_CS_do_sh_xy_q");
    CS.do_sh_xx_h             = captured.logical("_CS_do_sh_xx_h");
    CS.do_FrictWork           = captured.logical("_CS_do_FrictWork");
    CS.do_FrictWorkIntz       = captured.logical("_CS_do_FrictWorkIntz");
    CS.do_FrictWork_bh        = captured.logical("_CS_do_FrictWork_bh");
    CS.do_FrictWorkIntz_bh    = captured.logical("_CS_do_FrictWorkIntz_bh");
    CS.do_normstress          = captured.logical("_CS_do_normstress");
    CS.do_shearstress         = captured.logical("_CS_do_shearstress");
    CS.do_BS_coeff_h          = captured.logical("_CS_do_BS_coeff_h");
    CS.do_BS_coeff_q          = captured.logical("_CS_do_BS_coeff_q");

    auto       grad_vort_mag_h            = captured.fab_device("_grad_vort_mag_h_before");
    auto       grad_vort_mag_q            = captured.fab_device("_grad_vort_mag_q_before");
    auto       vert_vort_mag_smooth       = captured.fab_device("_vert_vort_mag_smooth_before");
    const auto vort_xy_dx                 = captured.fab_device("_vort_xy_dx_before");
    const auto vort_xy_dx_smooth          = captured.fab_device("_vort_xy_dx_smooth_before");
    const auto vort_xy_dy                 = captured.fab_device("_vort_xy_dy_before");
    const auto vort_xy_dy_smooth          = captured.fab_device("_vort_xy_dy_smooth_before");
    const auto grad_vort_mag_h_after      = captured.fab_host("_grad_vort_mag_h_after");
    const auto grad_vort_mag_q_after      = captured.fab_host("_grad_vort_mag_q_after");
    const auto vert_vort_mag_smooth_after = captured.fab_host("_vert_vort_mag_smooth_after");

    MOM::hor_visc_leith_grad_b(bxKh, bxQstr, CS,
                               grad_vort_mag_h.array(), grad_vort_mag_q.array(), vert_vort_mag_smooth.array(),
                               vort_xy_dx.const_array(), vort_xy_dx_smooth.const_array(),
                               vort_xy_dy.const_array(), vort_xy_dy_smooth.const_array());
    amrex::Gpu::synchronize();

    expect_arrays_equal(grad_vort_mag_h_after,      to_host_fab(grad_vort_mag_h),      "grad_vort_mag_h");
    expect_arrays_equal(grad_vort_mag_q_after,      to_host_fab(grad_vort_mag_q),      "grad_vort_mag_q");
    expect_arrays_equal(vert_vort_mag_smooth_after, to_host_fab(vert_vort_mag_smooth), "vert_vort_mag_smooth");
}

// -------------------------------------------------------------------------
// Field mapping: MOM6 src/parameterizations/lateral/MOM_hor_visc.F90:7693-7718 (hor_visc_Leithy_m capture block).
TEST(HorViscLeithyM, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "hor_visc_leithy_m");

    const auto bxKh = captured.box("_bxKh");

    hv_opts_CS_C CS{};
    CS.Laplacian              = captured.logical("_CS_Laplacian");
    CS.biharmonic             = captured.logical("_CS_biharmonic");
    CS.debug                  = captured.logical("_CS_debug");
    CS.no_slip                = captured.logical("_CS_no_slip");
    CS.bound_Kh               = captured.logical("_CS_bound_Kh");
    CS.EY24_EBT_BS            = captured.logical("_CS_EY24_EBT_BS");
    CS.bound_Ah               = captured.logical("_CS_bound_Ah");
    CS.Re_Ah                  = captured.real64("_CS_Re_Ah");
    CS.KS_coef                = captured.real64("_CS_KS_coef");
    CS.backscatter_underbound = captured.logical("_CS_backscatter_underbound");
    CS.Smagorinsky_Kh         = captured.logical("_CS_Smagorinsky_Kh");
    CS.Smagorinsky_Ah         = captured.logical("_CS_Smagorinsky_Ah");
    CS.Leith_Kh               = captured.logical("_CS_Leith_Kh");
    CS.Modified_Leith         = captured.logical("_CS_Modified_Leith");
    CS.use_beta_in_Leith      = captured.logical("_CS_use_beta_in_Leith");
    CS.Leith_Ah               = captured.logical("_CS_Leith_Ah");
    CS.use_Leithy             = captured.logical("_CS_use_Leithy");
    CS.c_K                    = captured.real64("_CS_c_K");
    CS.smooth_Ah              = captured.logical("_CS_smooth_Ah");
    CS.use_QG_Leith_visc      = captured.logical("_CS_use_QG_Leith_visc");
    CS.bound_Coriolis         = captured.logical("_CS_bound_Coriolis");
    CS.Kh_bg_min              = captured.real64("_CS_Kh_bg_min");
    CS.FrictWork_bug          = captured.logical("_CS_FrictWork_bug");
    CS.OBC_strain_bug         = captured.logical("_CS_OBC_strain_bug");
    CS.use_land_mask          = captured.logical("_CS_use_land_mask");
    CS.anisotropic            = captured.logical("_CS_anisotropic");
    CS.add_LES_viscosity      = captured.logical("_CS_add_LES_viscosity");
    CS.Kh_aniso               = captured.real64("_CS_Kh_aniso");
    CS.res_scale_MEKE         = captured.logical("_CS_res_scale_MEKE");
    CS.use_GME                = captured.logical("_CS_use_GME");
    CS.answer_date            = captured.integer("_CS_answer_date");
    CS.GME_h0                 = captured.real64("_CS_GME_h0");
    CS.GME_efficiency         = captured.real64("_CS_GME_efficiency");
    CS.GME_limiter            = captured.real64("_CS_GME_limiter");
    CS.min_grid_Kh            = captured.real64("_CS_min_grid_Kh");
    CS.min_grid_Ah            = captured.real64("_CS_min_grid_Ah");
    CS.use_cont_thick         = captured.logical("_CS_use_cont_thick");
    CS.use_ZB2020             = captured.logical("_CS_use_ZB2020");
    CS.use_circulation        = captured.logical("_CS_use_circulation");
    CS.do_grid_Re_Ah          = captured.logical("_CS_do_grid_Re_Ah");
    CS.do_grid_Re_Kh          = captured.logical("_CS_do_grid_Re_Kh");
    CS.do_Ah_h                = captured.logical("_CS_do_Ah_h");
    CS.do_Ah_q                = captured.logical("_CS_do_Ah_q");
    CS.do_Kh_h                = captured.logical("_CS_do_Kh_h");
    CS.do_Kh_q                = captured.logical("_CS_do_Kh_q");
    CS.do_GME_coeff_h         = captured.logical("_CS_do_GME_coeff_h");
    CS.do_GME_coeff_q         = captured.logical("_CS_do_GME_coeff_q");
    CS.do_vort_xy_q           = captured.logical("_CS_do_vort_xy_q");
    CS.do_div_xx_h            = captured.logical("_CS_do_div_xx_h");
    CS.do_sh_xy_q             = captured.logical("_CS_do_sh_xy_q");
    CS.do_sh_xx_h             = captured.logical("_CS_do_sh_xx_h");
    CS.do_FrictWork           = captured.logical("_CS_do_FrictWork");
    CS.do_FrictWorkIntz       = captured.logical("_CS_do_FrictWorkIntz");
    CS.do_FrictWork_bh        = captured.logical("_CS_do_FrictWork_bh");
    CS.do_FrictWorkIntz_bh    = captured.logical("_CS_do_FrictWorkIntz_bh");
    CS.do_normstress          = captured.logical("_CS_do_normstress");
    CS.do_shearstress         = captured.logical("_CS_do_shearstress");
    CS.do_BS_coeff_h          = captured.logical("_CS_do_BS_coeff_h");
    CS.do_BS_coeff_q          = captured.logical("_CS_do_BS_coeff_q");

    const auto mask2dBu = captured.fab_device("_mask2dBu");

    // The background biharmonic viscosity at h points [L4 T-1 ~> m4 s-1]; may be unallocated
    amrex::FArrayBox Ah_bg_xx_fab;
    amrex::Array4<const amrex::Real> Ah_bg_xx{};
    if (captured.logical("_metrics_Ah_bg_xx_allocated")) {
        Ah_bg_xx_fab = captured.fab_device("_metrics_Ah_bg_xx_before");
        Ah_bg_xx = Ah_bg_xx_fab.const_array();
    }
    // Pre-calculated .5*sqrt(c_K)*max{dx,dy} [L ~> m]; may be unallocated
    amrex::FArrayBox m_const_leithy_fab;
    amrex::Array4<const amrex::Real> m_const_leithy{};
    if (captured.logical("_metrics_m_const_leithy_allocated")) {
        m_const_leithy_fab = captured.fab_device("_metrics_m_const_leithy_before");
        m_const_leithy = m_const_leithy_fab.const_array();
    }
    // Pre-calculated 4./max(dx,dy)^2 at h points [L-2 ~> m-2]; may be unallocated
    amrex::FArrayBox m_leithy_max_fab;
    amrex::Array4<const amrex::Real> m_leithy_max{};
    if (captured.logical("_metrics_m_leithy_max_allocated")) {
        m_leithy_max_fab = captured.fab_device("_metrics_m_leithy_max_before");
        m_leithy_max = m_leithy_max_fab.const_array();
    }
    // Biharmonic metric-dependent constants [L6 ~> m6]; may be unallocated
    amrex::FArrayBox Biharm6_const_xx_fab;
    amrex::Array4<const amrex::Real> Biharm6_const_xx{};
    if (captured.logical("_metrics_Biharm6_const_xx_allocated")) {
        Biharm6_const_xx_fab = captured.fab_device("_metrics_Biharm6_const_xx_before");
        Biharm6_const_xx = Biharm6_const_xx_fab.const_array();
    }

    const double inv_PI6        = captured.real64("_inv_PI6");
    const auto   Del2vort_q     = captured.fab_device("_Del2vort_q");
    const auto   vert_vort_mag  = captured.fab_device("_vert_vort_mag");
    const auto   vort_xy_smooth = captured.fab_device("_vort_xy_smooth");
    auto         m_leithy       = captured.fab_device("_m_leithy_before");
    const auto   m_leithy_after = captured.fab_host("_m_leithy_after");

    MOM::hor_visc_leithy_m(bxKh, CS, mask2dBu.const_array(),
                           Ah_bg_xx, m_const_leithy, m_leithy_max, Biharm6_const_xx,
                           inv_PI6,
                           Del2vort_q.const_array(), vert_vort_mag.const_array(), vort_xy_smooth.const_array(),
                           m_leithy.array());
    amrex::Gpu::synchronize();

    expect_arrays_equal(m_leithy_after, to_host_fab(m_leithy), "m_leithy");
}

// -------------------------------------------------------------------------
// Field mapping: MOM6 src/parameterizations/lateral/MOM_hor_visc.F90:7838-7864 (hor_visc_Leithy_Ah capture block).
TEST(HorViscLeithyAh, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "hor_visc_leithy_ah");

    const auto bxKh = captured.box("_bxKh");

    hv_opts_CS_C CS{};
    CS.Laplacian              = captured.logical("_CS_Laplacian");
    CS.biharmonic             = captured.logical("_CS_biharmonic");
    CS.debug                  = captured.logical("_CS_debug");
    CS.no_slip                = captured.logical("_CS_no_slip");
    CS.bound_Kh               = captured.logical("_CS_bound_Kh");
    CS.EY24_EBT_BS            = captured.logical("_CS_EY24_EBT_BS");
    CS.bound_Ah               = captured.logical("_CS_bound_Ah");
    CS.Re_Ah                  = captured.real64("_CS_Re_Ah");
    CS.KS_coef                = captured.real64("_CS_KS_coef");
    CS.backscatter_underbound = captured.logical("_CS_backscatter_underbound");
    CS.Smagorinsky_Kh         = captured.logical("_CS_Smagorinsky_Kh");
    CS.Smagorinsky_Ah         = captured.logical("_CS_Smagorinsky_Ah");
    CS.Leith_Kh               = captured.logical("_CS_Leith_Kh");
    CS.Modified_Leith         = captured.logical("_CS_Modified_Leith");
    CS.use_beta_in_Leith      = captured.logical("_CS_use_beta_in_Leith");
    CS.Leith_Ah               = captured.logical("_CS_Leith_Ah");
    CS.use_Leithy             = captured.logical("_CS_use_Leithy");
    CS.c_K                    = captured.real64("_CS_c_K");
    CS.smooth_Ah              = captured.logical("_CS_smooth_Ah");
    CS.use_QG_Leith_visc      = captured.logical("_CS_use_QG_Leith_visc");
    CS.bound_Coriolis         = captured.logical("_CS_bound_Coriolis");
    CS.Kh_bg_min              = captured.real64("_CS_Kh_bg_min");
    CS.FrictWork_bug          = captured.logical("_CS_FrictWork_bug");
    CS.OBC_strain_bug         = captured.logical("_CS_OBC_strain_bug");
    CS.use_land_mask          = captured.logical("_CS_use_land_mask");
    CS.anisotropic            = captured.logical("_CS_anisotropic");
    CS.add_LES_viscosity      = captured.logical("_CS_add_LES_viscosity");
    CS.Kh_aniso               = captured.real64("_CS_Kh_aniso");
    CS.res_scale_MEKE         = captured.logical("_CS_res_scale_MEKE");
    CS.use_GME                = captured.logical("_CS_use_GME");
    CS.answer_date            = captured.integer("_CS_answer_date");
    CS.GME_h0                 = captured.real64("_CS_GME_h0");
    CS.GME_efficiency         = captured.real64("_CS_GME_efficiency");
    CS.GME_limiter            = captured.real64("_CS_GME_limiter");
    CS.min_grid_Kh            = captured.real64("_CS_min_grid_Kh");
    CS.min_grid_Ah            = captured.real64("_CS_min_grid_Ah");
    CS.use_cont_thick         = captured.logical("_CS_use_cont_thick");
    CS.use_ZB2020             = captured.logical("_CS_use_ZB2020");
    CS.use_circulation        = captured.logical("_CS_use_circulation");
    CS.do_grid_Re_Ah          = captured.logical("_CS_do_grid_Re_Ah");
    CS.do_grid_Re_Kh          = captured.logical("_CS_do_grid_Re_Kh");
    CS.do_Ah_h                = captured.logical("_CS_do_Ah_h");
    CS.do_Ah_q                = captured.logical("_CS_do_Ah_q");
    CS.do_Kh_h                = captured.logical("_CS_do_Kh_h");
    CS.do_Kh_q                = captured.logical("_CS_do_Kh_q");
    CS.do_GME_coeff_h         = captured.logical("_CS_do_GME_coeff_h");
    CS.do_GME_coeff_q         = captured.logical("_CS_do_GME_coeff_q");
    CS.do_vort_xy_q           = captured.logical("_CS_do_vort_xy_q");
    CS.do_div_xx_h            = captured.logical("_CS_do_div_xx_h");
    CS.do_sh_xy_q             = captured.logical("_CS_do_sh_xy_q");
    CS.do_sh_xx_h             = captured.logical("_CS_do_sh_xx_h");
    CS.do_FrictWork           = captured.logical("_CS_do_FrictWork");
    CS.do_FrictWorkIntz       = captured.logical("_CS_do_FrictWorkIntz");
    CS.do_FrictWork_bh        = captured.logical("_CS_do_FrictWork_bh");
    CS.do_FrictWorkIntz_bh    = captured.logical("_CS_do_FrictWorkIntz_bh");
    CS.do_normstress          = captured.logical("_CS_do_normstress");
    CS.do_shearstress         = captured.logical("_CS_do_shearstress");
    CS.do_BS_coeff_h          = captured.logical("_CS_do_BS_coeff_h");
    CS.do_BS_coeff_q          = captured.logical("_CS_do_BS_coeff_q");

    // The background biharmonic viscosity at h points [L4 T-1 ~> m4 s-1]; may be unallocated
    amrex::FArrayBox Ah_bg_xx_fab;
    amrex::Array4<const amrex::Real> Ah_bg_xx{};
    if (captured.logical("_metrics_Ah_bg_xx_allocated")) {
        Ah_bg_xx_fab = captured.fab_device("_metrics_Ah_bg_xx_before");
        Ah_bg_xx = Ah_bg_xx_fab.const_array();
    }
    // Biharmonic metric-dependent constants [L6 ~> m6]; may be unallocated
    amrex::FArrayBox Biharm6_const_xx_fab;
    amrex::Array4<const amrex::Real> Biharm6_const_xx{};
    if (captured.logical("_metrics_Biharm6_const_xx_allocated")) {
        Biharm6_const_xx_fab = captured.fab_device("_metrics_Biharm6_const_xx_before");
        Biharm6_const_xx = Biharm6_const_xx_fab.const_array();
    }

    const double inv_PI6              = captured.real64("_inv_PI6");
    const auto   Del2vort_q           = captured.fab_device("_Del2vort_q");
    const auto   vert_vort_mag_smooth = captured.fab_device("_vert_vort_mag_smooth");
    auto         Ah                   = captured.fab_device("_Ah_before");
    const auto   m_leithy             = captured.fab_device("_m_leithy");
    auto         Ah_sq                = captured.fab_device("_Ah_sq_before");
    const auto   Ah_after             = captured.fab_host("_Ah_after");
    const auto   Ah_sq_after          = captured.fab_host("_Ah_sq_after");

    MOM::hor_visc_leithy_ah(bxKh, CS,
                            Ah_bg_xx, Biharm6_const_xx,
                            inv_PI6,
                            Del2vort_q.const_array(), vert_vort_mag_smooth.const_array(),
                            Ah.array(), m_leithy.const_array(), Ah_sq.array());
    amrex::Gpu::synchronize();

    expect_arrays_equal(Ah_after,    to_host_fab(Ah),    "Ah");
    expect_arrays_equal(Ah_sq_after, to_host_fab(Ah_sq), "Ah_sq");
}

// -------------------------------------------------------------------------
// Field mapping: MOM6 src/parameterizations/lateral/MOM_hor_visc.F90:7964-7990 (hor_visc_Leithy_Ah_h capture block).
TEST(HorViscLeithyAhH, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "hor_visc_leithy_ah_h");

    const auto bxKh = captured.box("_bxKh");

    hv_opts_CS_C CS{};
    CS.Laplacian              = captured.logical("_CS_Laplacian");
    CS.biharmonic             = captured.logical("_CS_biharmonic");
    CS.debug                  = captured.logical("_CS_debug");
    CS.no_slip                = captured.logical("_CS_no_slip");
    CS.bound_Kh               = captured.logical("_CS_bound_Kh");
    CS.EY24_EBT_BS            = captured.logical("_CS_EY24_EBT_BS");
    CS.bound_Ah               = captured.logical("_CS_bound_Ah");
    CS.Re_Ah                  = captured.real64("_CS_Re_Ah");
    CS.KS_coef                = captured.real64("_CS_KS_coef");
    CS.backscatter_underbound = captured.logical("_CS_backscatter_underbound");
    CS.Smagorinsky_Kh         = captured.logical("_CS_Smagorinsky_Kh");
    CS.Smagorinsky_Ah         = captured.logical("_CS_Smagorinsky_Ah");
    CS.Leith_Kh               = captured.logical("_CS_Leith_Kh");
    CS.Modified_Leith         = captured.logical("_CS_Modified_Leith");
    CS.use_beta_in_Leith      = captured.logical("_CS_use_beta_in_Leith");
    CS.Leith_Ah               = captured.logical("_CS_Leith_Ah");
    CS.use_Leithy             = captured.logical("_CS_use_Leithy");
    CS.c_K                    = captured.real64("_CS_c_K");
    CS.smooth_Ah              = captured.logical("_CS_smooth_Ah");
    CS.use_QG_Leith_visc      = captured.logical("_CS_use_QG_Leith_visc");
    CS.bound_Coriolis         = captured.logical("_CS_bound_Coriolis");
    CS.Kh_bg_min              = captured.real64("_CS_Kh_bg_min");
    CS.FrictWork_bug          = captured.logical("_CS_FrictWork_bug");
    CS.OBC_strain_bug         = captured.logical("_CS_OBC_strain_bug");
    CS.use_land_mask          = captured.logical("_CS_use_land_mask");
    CS.anisotropic            = captured.logical("_CS_anisotropic");
    CS.add_LES_viscosity      = captured.logical("_CS_add_LES_viscosity");
    CS.Kh_aniso               = captured.real64("_CS_Kh_aniso");
    CS.res_scale_MEKE         = captured.logical("_CS_res_scale_MEKE");
    CS.use_GME                = captured.logical("_CS_use_GME");
    CS.answer_date            = captured.integer("_CS_answer_date");
    CS.GME_h0                 = captured.real64("_CS_GME_h0");
    CS.GME_efficiency         = captured.real64("_CS_GME_efficiency");
    CS.GME_limiter            = captured.real64("_CS_GME_limiter");
    CS.min_grid_Kh            = captured.real64("_CS_min_grid_Kh");
    CS.min_grid_Ah            = captured.real64("_CS_min_grid_Ah");
    CS.use_cont_thick         = captured.logical("_CS_use_cont_thick");
    CS.use_ZB2020             = captured.logical("_CS_use_ZB2020");
    CS.use_circulation        = captured.logical("_CS_use_circulation");
    CS.do_grid_Re_Ah          = captured.logical("_CS_do_grid_Re_Ah");
    CS.do_grid_Re_Kh          = captured.logical("_CS_do_grid_Re_Kh");
    CS.do_Ah_h                = captured.logical("_CS_do_Ah_h");
    CS.do_Ah_q                = captured.logical("_CS_do_Ah_q");
    CS.do_Kh_h                = captured.logical("_CS_do_Kh_h");
    CS.do_Kh_q                = captured.logical("_CS_do_Kh_q");
    CS.do_GME_coeff_h         = captured.logical("_CS_do_GME_coeff_h");
    CS.do_GME_coeff_q         = captured.logical("_CS_do_GME_coeff_q");
    CS.do_vort_xy_q           = captured.logical("_CS_do_vort_xy_q");
    CS.do_div_xx_h            = captured.logical("_CS_do_div_xx_h");
    CS.do_sh_xy_q             = captured.logical("_CS_do_sh_xy_q");
    CS.do_sh_xx_h             = captured.logical("_CS_do_sh_xx_h");
    CS.do_FrictWork           = captured.logical("_CS_do_FrictWork");
    CS.do_FrictWorkIntz       = captured.logical("_CS_do_FrictWorkIntz");
    CS.do_FrictWork_bh        = captured.logical("_CS_do_FrictWork_bh");
    CS.do_FrictWorkIntz_bh    = captured.logical("_CS_do_FrictWorkIntz_bh");
    CS.do_normstress          = captured.logical("_CS_do_normstress");
    CS.do_shearstress         = captured.logical("_CS_do_shearstress");
    CS.do_BS_coeff_h          = captured.logical("_CS_do_BS_coeff_h");
    CS.do_BS_coeff_q          = captured.logical("_CS_do_BS_coeff_q");

    // The background biharmonic viscosity at h points [L4 T-1 ~> m4 s-1]; may be unallocated
    amrex::FArrayBox Ah_bg_xx_fab;
    amrex::Array4<const amrex::Real> Ah_bg_xx{};
    if (captured.logical("_metrics_Ah_bg_xx_allocated")) {
        Ah_bg_xx_fab = captured.fab_device("_metrics_Ah_bg_xx_before");
        Ah_bg_xx = Ah_bg_xx_fab.const_array();
    }

    const int  kstart     = captured.integer("_kstart");
    auto       Ah         = captured.fab_device("_Ah_before");
    auto       Ah_h       = captured.fab_device("_Ah_h_before");
    const auto Ah_sq      = captured.fab_device("_Ah_sq");
    const auto Ah_after   = captured.fab_host("_Ah_after");
    const auto Ah_h_after = captured.fab_host("_Ah_h_after");

    MOM::hor_visc_leithy_ah_h(bxKh, CS,
                              Ah_bg_xx,
                              kstart,
                              Ah.array(), Ah_h.array(), Ah_sq.const_array());
    amrex::Gpu::synchronize();

    expect_arrays_equal(Ah_after,   to_host_fab(Ah),   "Ah");
    expect_arrays_equal(Ah_h_after, to_host_fab(Ah_h), "Ah_h");
}

// -------------------------------------------------------------------------
// Field mapping: MOM6 src/parameterizations/lateral/MOM_hor_visc.F90:3529-3626 (hor_visc_strain capture block).
TEST(HorViscStrain, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "hor_visc_strain");

    const auto bxH        = captured.box("_bxH");
    const auto bxQ        = captured.box("_bxQ");
    const auto bxU        = captured.box("_bxU");
    const auto bxV        = captured.box("_bxV");
    const auto bxKh       = captured.box("_bxKh");
    const auto bxHstr     = captured.box("_bxHstr");
    const auto bxQstr     = captured.box("_bxQstr");
    const auto bxHgrad    = captured.box("_bxHgrad");
    const auto bxQshear   = captured.box("_bxQshear");
    const auto bxQvsmooth = captured.box("_bxQvsmooth");

    hv_opts_CS_C CS{};
    CS.Laplacian              = captured.logical("_CS_Laplacian");
    CS.biharmonic             = captured.logical("_CS_biharmonic");
    CS.debug                  = captured.logical("_CS_debug");
    CS.no_slip                = captured.logical("_CS_no_slip");
    CS.bound_Kh               = captured.logical("_CS_bound_Kh");
    CS.EY24_EBT_BS            = captured.logical("_CS_EY24_EBT_BS");
    CS.bound_Ah               = captured.logical("_CS_bound_Ah");
    CS.Re_Ah                  = captured.real64("_CS_Re_Ah");
    CS.KS_coef                = captured.real64("_CS_KS_coef");
    CS.backscatter_underbound = captured.logical("_CS_backscatter_underbound");
    CS.Smagorinsky_Kh         = captured.logical("_CS_Smagorinsky_Kh");
    CS.Smagorinsky_Ah         = captured.logical("_CS_Smagorinsky_Ah");
    CS.Leith_Kh               = captured.logical("_CS_Leith_Kh");
    CS.Modified_Leith         = captured.logical("_CS_Modified_Leith");
    CS.use_beta_in_Leith      = captured.logical("_CS_use_beta_in_Leith");
    CS.Leith_Ah               = captured.logical("_CS_Leith_Ah");
    CS.use_Leithy             = captured.logical("_CS_use_Leithy");
    CS.c_K                    = captured.real64("_CS_c_K");
    CS.smooth_Ah              = captured.logical("_CS_smooth_Ah");
    CS.use_QG_Leith_visc      = captured.logical("_CS_use_QG_Leith_visc");
    CS.bound_Coriolis         = captured.logical("_CS_bound_Coriolis");
    CS.Kh_bg_min              = captured.real64("_CS_Kh_bg_min");
    CS.FrictWork_bug          = captured.logical("_CS_FrictWork_bug");
    CS.OBC_strain_bug         = captured.logical("_CS_OBC_strain_bug");
    CS.use_land_mask          = captured.logical("_CS_use_land_mask");
    CS.anisotropic            = captured.logical("_CS_anisotropic");
    CS.add_LES_viscosity      = captured.logical("_CS_add_LES_viscosity");
    CS.Kh_aniso               = captured.real64("_CS_Kh_aniso");
    CS.res_scale_MEKE         = captured.logical("_CS_res_scale_MEKE");
    CS.use_GME                = captured.logical("_CS_use_GME");
    CS.answer_date            = captured.integer("_CS_answer_date");
    CS.GME_h0                 = captured.real64("_CS_GME_h0");
    CS.GME_efficiency         = captured.real64("_CS_GME_efficiency");
    CS.GME_limiter            = captured.real64("_CS_GME_limiter");
    CS.min_grid_Kh            = captured.real64("_CS_min_grid_Kh");
    CS.min_grid_Ah            = captured.real64("_CS_min_grid_Ah");
    CS.use_cont_thick         = captured.logical("_CS_use_cont_thick");
    CS.use_ZB2020             = captured.logical("_CS_use_ZB2020");
    CS.use_circulation        = captured.logical("_CS_use_circulation");
    CS.do_grid_Re_Ah          = captured.logical("_CS_do_grid_Re_Ah");
    CS.do_grid_Re_Kh          = captured.logical("_CS_do_grid_Re_Kh");
    CS.do_Ah_h                = captured.logical("_CS_do_Ah_h");
    CS.do_Ah_q                = captured.logical("_CS_do_Ah_q");
    CS.do_Kh_h                = captured.logical("_CS_do_Kh_h");
    CS.do_Kh_q                = captured.logical("_CS_do_Kh_q");
    CS.do_GME_coeff_h         = captured.logical("_CS_do_GME_coeff_h");
    CS.do_GME_coeff_q         = captured.logical("_CS_do_GME_coeff_q");
    CS.do_vort_xy_q           = captured.logical("_CS_do_vort_xy_q");
    CS.do_div_xx_h            = captured.logical("_CS_do_div_xx_h");
    CS.do_sh_xy_q             = captured.logical("_CS_do_sh_xy_q");
    CS.do_sh_xx_h             = captured.logical("_CS_do_sh_xx_h");
    CS.do_FrictWork           = captured.logical("_CS_do_FrictWork");
    CS.do_FrictWorkIntz       = captured.logical("_CS_do_FrictWorkIntz");
    CS.do_FrictWork_bh        = captured.logical("_CS_do_FrictWork_bh");
    CS.do_FrictWorkIntz_bh    = captured.logical("_CS_do_FrictWorkIntz_bh");
    CS.do_normstress          = captured.logical("_CS_do_normstress");
    CS.do_shearstress         = captured.logical("_CS_do_shearstress");
    CS.do_BS_coeff_h          = captured.logical("_CS_do_BS_coeff_h");
    CS.do_BS_coeff_q          = captured.logical("_CS_do_BS_coeff_q");

    const auto IdyCu    = captured.fab_device("_IdyCu");
    const auto IdxCv    = captured.fab_device("_IdxCv");
    const auto IdyCv    = captured.fab_device("_IdyCv");
    const auto IdxCu    = captured.fab_device("_IdxCu");
    const auto mask2dT  = captured.fab_device("_mask2dT");
    const auto dxBu     = captured.fab_device("_dxBu");
    const auto mask2dBu = captured.fab_device("_mask2dBu");
    const auto IareaBu  = captured.fab_device("_IareaBu");
    const auto dyCv     = captured.fab_device("_dyCv");
    const auto dxCu     = captured.fab_device("_dxCu");

    // Pre-calculated dx^2 at h points [L2 ~> m2]; may be unallocated
    amrex::FArrayBox dx2h_fab;
    amrex::Array4<const amrex::Real> dx2h{};
    if (captured.logical("_metrics_dx2h_allocated")) {
        dx2h_fab = captured.fab_device("_metrics_dx2h_before");
        dx2h = dx2h_fab.const_array();
    }

    // Pre-calculated dy^2 at h points [L2 ~> m2]; may be unallocated
    amrex::FArrayBox dy2h_fab;
    amrex::Array4<const amrex::Real> dy2h{};
    if (captured.logical("_metrics_dy2h_allocated")) {
        dy2h_fab = captured.fab_device("_metrics_dy2h_before");
        dy2h = dy2h_fab.const_array();
    }

    // Pre-calculated dx/dy at h points [nondim]; may be unallocated
    amrex::FArrayBox dx_dyT_fab;
    amrex::Array4<const amrex::Real> dx_dyT{};
    if (captured.logical("_metrics_dx_dyT_allocated")) {
        dx_dyT_fab = captured.fab_device("_metrics_dx_dyT_before");
        dx_dyT = dx_dyT_fab.const_array();
    }

    // Pre-calculated dy/dx at h points [nondim]; may be unallocated
    amrex::FArrayBox dy_dxT_fab;
    amrex::Array4<const amrex::Real> dy_dxT{};
    if (captured.logical("_metrics_dy_dxT_allocated")) {
        dy_dxT_fab = captured.fab_device("_metrics_dy_dxT_before");
        dy_dxT = dy_dxT_fab.const_array();
    }

    // Pre-calculated dx^2 at q points [L2 ~> m2]; may be unallocated
    amrex::FArrayBox dx2q_fab;
    amrex::Array4<const amrex::Real> dx2q{};
    if (captured.logical("_metrics_dx2q_allocated")) {
        dx2q_fab = captured.fab_device("_metrics_dx2q_before");
        dx2q = dx2q_fab.const_array();
    }

    // Pre-calculated dy^2 at q points [L2 ~> m2]; may be unallocated
    amrex::FArrayBox dy2q_fab;
    amrex::Array4<const amrex::Real> dy2q{};
    if (captured.logical("_metrics_dy2q_allocated")) {
        dy2q_fab = captured.fab_device("_metrics_dy2q_before");
        dy2q = dy2q_fab.const_array();
    }

    // Pre-calculated dx/dy at q points [nondim]; may be unallocated
    amrex::FArrayBox dx_dyBu_fab;
    amrex::Array4<const amrex::Real> dx_dyBu{};
    if (captured.logical("_metrics_dx_dyBu_allocated")) {
        dx_dyBu_fab = captured.fab_device("_metrics_dx_dyBu_before");
        dx_dyBu = dx_dyBu_fab.const_array();
    }

    // Pre-calculated dy/dx at q points [nondim]; may be unallocated
    amrex::FArrayBox dy_dxBu_fab;
    amrex::Array4<const amrex::Real> dy_dxBu{};
    if (captured.logical("_metrics_dy_dxBu_allocated")) {
        dy_dxBu_fab = captured.fab_device("_metrics_dy_dxBu_before");
        dy_dxBu = dy_dxBu_fab.const_array();
    }

    // 1/(dx^2 dy) at u points [L-3 ~> m-3]; may be unallocated
    amrex::FArrayBox Idx2dyCu_fab;
    amrex::Array4<const amrex::Real> Idx2dyCu{};
    if (captured.logical("_metrics_Idx2dyCu_allocated")) {
        Idx2dyCu_fab = captured.fab_device("_metrics_Idx2dyCu_before");
        Idx2dyCu = Idx2dyCu_fab.const_array();
    }

    // 1/(dx dy^2) at u points [L-3 ~> m-3]; may be unallocated
    amrex::FArrayBox Idxdy2u_fab;
    amrex::Array4<const amrex::Real> Idxdy2u{};
    if (captured.logical("_metrics_Idxdy2u_allocated")) {
        Idxdy2u_fab = captured.fab_device("_metrics_Idxdy2u_before");
        Idxdy2u = Idxdy2u_fab.const_array();
    }

    // 1/(dx^2 dy) at v points [L-3 ~> m-3]; may be unallocated
    amrex::FArrayBox Idx2dyCv_fab;
    amrex::Array4<const amrex::Real> Idx2dyCv{};
    if (captured.logical("_metrics_Idx2dyCv_allocated")) {
        Idx2dyCv_fab = captured.fab_device("_metrics_Idx2dyCv_before");
        Idx2dyCv = Idx2dyCv_fab.const_array();
    }

    // 1/(dx dy^2) at v points [L-3 ~> m-3]; may be unallocated
    amrex::FArrayBox Idxdy2v_fab;
    amrex::Array4<const amrex::Real> Idxdy2v{};
    if (captured.logical("_metrics_Idxdy2v_allocated")) {
        Idxdy2v_fab = captured.fab_device("_metrics_Idxdy2v_before");
        Idxdy2v = Idxdy2v_fab.const_array();
    }

    const int  kstart           = captured.integer("_kstart");
    const bool use_vort_xy      = captured.logical("_use_vort_xy");
    const bool apply_OBC        = captured.logical("_apply_OBC");
    const bool apply_OBC_strain = captured.logical("_apply_OBC_strain");
    const bool use_cont_huv     = captured.logical("_use_cont_huv");
    // Loop ranges for vorticity terms

    const auto u        = captured.fab_device("_u");
    const auto v        = captured.fab_device("_v");
    const auto h        = captured.fab_device("_h");

    // Layer thickness at u-points [H ~> m or kg m-2]; may be unallocated
    amrex::FArrayBox hu_cont_fab;
    amrex::Array4<const amrex::Real> hu_cont{};
    if (captured.logical("_hu_cont_allocated")) {
        hu_cont_fab = captured.fab_device("_hu_cont_before");
        hu_cont = hu_cont_fab.const_array();
    }

    // Layer thickness at v-points [H ~> m or kg m-2]; may be unallocated
    amrex::FArrayBox hv_cont_fab;
    amrex::Array4<const amrex::Real> hv_cont{};
    if (captured.logical("_hv_cont_allocated")) {
        hv_cont_fab = captured.fab_device("_hv_cont_before");
        hv_cont = hv_cont_fab.const_array();
    }

    const auto u_smooth = captured.fab_device("_u_smooth");
    const auto v_smooth = captured.fab_device("_v_smooth");

    auto       NoSt                 = captured.fab_device("_NoSt_before");
    auto       ShSt                 = captured.fab_device("_ShSt_before");
    auto       Del2u                = captured.fab_device("_Del2u_before");
    auto       h_u                  = captured.fab_device("_h_u_before");
    auto       Del2v                = captured.fab_device("_Del2v_before");
    auto       h_v                  = captured.fab_device("_h_v_before");
    auto       sh_xx                = captured.fab_device("_sh_xx_before");
    auto       sh_xx_smooth         = captured.fab_device("_sh_xx_smooth_before");
    auto       dudx                 = captured.fab_device("_dudx_before");
    auto       dvdy                 = captured.fab_device("_dvdy_before");
    auto       dvdx                 = captured.fab_device("_dvdx_before");
    auto       dudy                 = captured.fab_device("_dudy_before");
    auto       dvdx_smooth          = captured.fab_device("_dvdx_smooth_before");
    auto       dudy_smooth          = captured.fab_device("_dudy_smooth_before");
    auto       sh_xy                = captured.fab_device("_sh_xy_before");
    auto       sh_xy_smooth         = captured.fab_device("_sh_xy_smooth_before");
    auto       vort_xy              = captured.fab_device("_vort_xy_before");
    auto       vort_xy_smooth       = captured.fab_device("_vort_xy_smooth_before");
    const auto NoSt_after           = captured.fab_host("_NoSt_after");
    const auto ShSt_after           = captured.fab_host("_ShSt_after");
    const auto Del2u_after          = captured.fab_host("_Del2u_after");
    const auto h_u_after            = captured.fab_host("_h_u_after");
    const auto Del2v_after          = captured.fab_host("_Del2v_after");
    const auto h_v_after            = captured.fab_host("_h_v_after");
    const auto sh_xx_after          = captured.fab_host("_sh_xx_after");
    const auto sh_xx_smooth_after   = captured.fab_host("_sh_xx_smooth_after");
    const auto dudx_after           = captured.fab_host("_dudx_after");
    const auto dvdy_after           = captured.fab_host("_dvdy_after");
    const auto dvdx_after           = captured.fab_host("_dvdx_after");
    const auto dudy_after           = captured.fab_host("_dudy_after");
    const auto dvdx_smooth_after    = captured.fab_host("_dvdx_smooth_after");
    const auto dudy_smooth_after    = captured.fab_host("_dudy_smooth_after");
    const auto sh_xy_after          = captured.fab_host("_sh_xy_after");
    const auto sh_xy_smooth_after   = captured.fab_host("_sh_xy_smooth_after");
    const auto vort_xy_after        = captured.fab_host("_vort_xy_after");
    const auto vort_xy_smooth_after = captured.fab_host("_vort_xy_smooth_after");

    MOM::hor_visc_strain(bxH, bxQ, bxU, bxV, bxKh, bxHstr, bxQstr, bxHgrad, bxQshear, bxQvsmooth,
                         CS,
                         IdyCu.const_array(), IdxCv.const_array(), IdyCv.const_array(), IdxCu.const_array(),
                         mask2dT.const_array(), dxBu.const_array(), mask2dBu.const_array(), IareaBu.const_array(),
                         dyCv.const_array(), dxCu.const_array(),
                         dx2h, dy2h, dx_dyT, dy_dxT, dx2q, dy2q, dx_dyBu, dy_dxBu,
                         Idx2dyCu, Idxdy2u, Idx2dyCv, Idxdy2v,
                         nullptr,
                         kstart,
                         use_vort_xy, apply_OBC, apply_OBC_strain, use_cont_huv,
                         u.const_array(), v.const_array(), h.const_array(),
                         hu_cont, hv_cont,
                         u_smooth.const_array(), v_smooth.const_array(),
                         NoSt.array(), ShSt.array(), Del2u.array(), h_u.array(),
                         Del2v.array(), h_v.array(), sh_xx.array(), sh_xx_smooth.array(),
                         dudx.array(), dvdy.array(), dvdx.array(), dudy.array(),
                         dvdx_smooth.array(), dudy_smooth.array(), sh_xy.array(), sh_xy_smooth.array(),
                         vort_xy.array(), vort_xy_smooth.array());
    amrex::Gpu::synchronize();

    expect_arrays_equal(NoSt_after,           to_host_fab(NoSt),           "NoSt");
    expect_arrays_equal(ShSt_after,           to_host_fab(ShSt),           "ShSt");
    expect_arrays_equal(Del2u_after,          to_host_fab(Del2u),          "Del2u");
    expect_arrays_equal(h_u_after,            to_host_fab(h_u),            "h_u");
    expect_arrays_equal(Del2v_after,          to_host_fab(Del2v),          "Del2v");
    expect_arrays_equal(h_v_after,            to_host_fab(h_v),            "h_v");
    expect_arrays_equal(sh_xx_after,          to_host_fab(sh_xx),          "sh_xx");
    expect_arrays_equal(sh_xx_smooth_after,   to_host_fab(sh_xx_smooth),   "sh_xx_smooth");
    expect_arrays_equal(dudx_after,           to_host_fab(dudx),           "dudx");
    expect_arrays_equal(dvdy_after,           to_host_fab(dvdy),           "dvdy");
    expect_arrays_equal(dvdx_after,           to_host_fab(dvdx),           "dvdx");
    expect_arrays_equal(dudy_after,           to_host_fab(dudy),           "dudy");
    expect_arrays_equal(dvdx_smooth_after,    to_host_fab(dvdx_smooth),    "dvdx_smooth");
    expect_arrays_equal(dudy_smooth_after,    to_host_fab(dudy_smooth),    "dudy_smooth");
    expect_arrays_equal(sh_xy_after,          to_host_fab(sh_xy),          "sh_xy");
    expect_arrays_equal(sh_xy_smooth_after,   to_host_fab(sh_xy_smooth),   "sh_xy_smooth");
    expect_arrays_equal(vort_xy_after,        to_host_fab(vort_xy),        "vort_xy");
    expect_arrays_equal(vort_xy_smooth_after, to_host_fab(vort_xy_smooth), "vort_xy_smooth");
}

// -------------------------------------------------------------------------
// Field mapping: MOM6 src/parameterizations/lateral/MOM_hor_visc.F90:4164-4248 (hor_visc_h_coeffs_a capture block).
TEST(HorViscHCoeffsA, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "hor_visc_h_coeffs_a");

    const auto bxH    = captured.box("_bxH");
    const auto bxKh   = captured.box("_bxKh");
    const auto bxHstr = captured.box("_bxHstr");

    hv_opts_CS_C CS{};
    CS.Laplacian              = captured.logical("_CS_Laplacian");
    CS.biharmonic             = captured.logical("_CS_biharmonic");
    CS.debug                  = captured.logical("_CS_debug");
    CS.no_slip                = captured.logical("_CS_no_slip");
    CS.bound_Kh               = captured.logical("_CS_bound_Kh");
    CS.EY24_EBT_BS            = captured.logical("_CS_EY24_EBT_BS");
    CS.bound_Ah               = captured.logical("_CS_bound_Ah");
    CS.Re_Ah                  = captured.real64("_CS_Re_Ah");
    CS.KS_coef                = captured.real64("_CS_KS_coef");
    CS.backscatter_underbound = captured.logical("_CS_backscatter_underbound");
    CS.Smagorinsky_Kh         = captured.logical("_CS_Smagorinsky_Kh");
    CS.Smagorinsky_Ah         = captured.logical("_CS_Smagorinsky_Ah");
    CS.Leith_Kh               = captured.logical("_CS_Leith_Kh");
    CS.Modified_Leith         = captured.logical("_CS_Modified_Leith");
    CS.use_beta_in_Leith      = captured.logical("_CS_use_beta_in_Leith");
    CS.Leith_Ah               = captured.logical("_CS_Leith_Ah");
    CS.use_Leithy             = captured.logical("_CS_use_Leithy");
    CS.c_K                    = captured.real64("_CS_c_K");
    CS.smooth_Ah              = captured.logical("_CS_smooth_Ah");
    CS.use_QG_Leith_visc      = captured.logical("_CS_use_QG_Leith_visc");
    CS.bound_Coriolis         = captured.logical("_CS_bound_Coriolis");
    CS.Kh_bg_min              = captured.real64("_CS_Kh_bg_min");
    CS.FrictWork_bug          = captured.logical("_CS_FrictWork_bug");
    CS.OBC_strain_bug         = captured.logical("_CS_OBC_strain_bug");
    CS.use_land_mask          = captured.logical("_CS_use_land_mask");
    CS.anisotropic            = captured.logical("_CS_anisotropic");
    CS.add_LES_viscosity      = captured.logical("_CS_add_LES_viscosity");
    CS.Kh_aniso               = captured.real64("_CS_Kh_aniso");
    CS.res_scale_MEKE         = captured.logical("_CS_res_scale_MEKE");
    CS.use_GME                = captured.logical("_CS_use_GME");
    CS.answer_date            = captured.integer("_CS_answer_date");
    CS.GME_h0                 = captured.real64("_CS_GME_h0");
    CS.GME_efficiency         = captured.real64("_CS_GME_efficiency");
    CS.GME_limiter            = captured.real64("_CS_GME_limiter");
    CS.min_grid_Kh            = captured.real64("_CS_min_grid_Kh");
    CS.min_grid_Ah            = captured.real64("_CS_min_grid_Ah");
    CS.use_cont_thick         = captured.logical("_CS_use_cont_thick");
    CS.use_ZB2020             = captured.logical("_CS_use_ZB2020");
    CS.use_circulation        = captured.logical("_CS_use_circulation");
    CS.do_grid_Re_Ah          = captured.logical("_CS_do_grid_Re_Ah");
    CS.do_grid_Re_Kh          = captured.logical("_CS_do_grid_Re_Kh");
    CS.do_Ah_h                = captured.logical("_CS_do_Ah_h");
    CS.do_Ah_q                = captured.logical("_CS_do_Ah_q");
    CS.do_Kh_h                = captured.logical("_CS_do_Kh_h");
    CS.do_Kh_q                = captured.logical("_CS_do_Kh_q");
    CS.do_GME_coeff_h         = captured.logical("_CS_do_GME_coeff_h");
    CS.do_GME_coeff_q         = captured.logical("_CS_do_GME_coeff_q");
    CS.do_vort_xy_q           = captured.logical("_CS_do_vort_xy_q");
    CS.do_div_xx_h            = captured.logical("_CS_do_div_xx_h");
    CS.do_sh_xy_q             = captured.logical("_CS_do_sh_xy_q");
    CS.do_sh_xx_h             = captured.logical("_CS_do_sh_xx_h");
    CS.do_FrictWork           = captured.logical("_CS_do_FrictWork");
    CS.do_FrictWorkIntz       = captured.logical("_CS_do_FrictWorkIntz");
    CS.do_FrictWork_bh        = captured.logical("_CS_do_FrictWork_bh");
    CS.do_FrictWorkIntz_bh    = captured.logical("_CS_do_FrictWorkIntz_bh");
    CS.do_normstress          = captured.logical("_CS_do_normstress");
    CS.do_shearstress         = captured.logical("_CS_do_shearstress");
    CS.do_BS_coeff_h          = captured.logical("_CS_do_BS_coeff_h");
    CS.do_BS_coeff_q          = captured.logical("_CS_do_BS_coeff_q");

    // The MEKE-derived lateral viscosity coefficient [L2 T-1 ~> m2 s-1]; may be unallocated
    amrex::FArrayBox MEKE_Ku_fab;
    amrex::Array4<const amrex::Real> MEKE_Ku{};
    if (captured.logical("_MEKE_Ku_allocated")) {
        MEKE_Ku_fab = captured.fab_device("_MEKE_Ku_before");
        MEKE_Ku = MEKE_Ku_fab.const_array();
    }

    // Non-dimensional function of the ratio the first baroclinic deformation radius to the grid spacing
    // at h points [nondim]; may be unallocated
    amrex::FArrayBox Res_fn_h_fab;
    amrex::Array4<const amrex::Real> Res_fn_h{};
    if (captured.logical("_VarMix_Res_fn_h_allocated")) {
        Res_fn_h_fab = captured.fab_device("_VarMix_Res_fn_h_before");
        Res_fn_h = Res_fn_h_fab.const_array();
    }

    // Vertical structure function used in backscatter [nondim]; may be unallocated
    amrex::FArrayBox BS_struct_fab;
    amrex::Array4<const amrex::Real> BS_struct{};
    if (captured.logical("_VarMix_BS_struct_allocated")) {
        BS_struct_fab = captured.fab_device("_VarMix_BS_struct_before");
        BS_struct = BS_struct_fab.const_array();
    }

    // The background Laplacian viscosity at h points [L2 T-1 ~> m2 s-1]; may be unallocated
    amrex::FArrayBox Kh_bg_xx_fab;
    amrex::Array4<const amrex::Real> Kh_bg_xx{};
    if (captured.logical("_metrics_Kh_bg_xx_allocated")) {
        Kh_bg_xx_fab = captured.fab_device("_metrics_Kh_bg_xx_before");
        Kh_bg_xx = Kh_bg_xx_fab.const_array();
    }

    // The background biharmonic viscosity at h points [L4 T-1 ~> m4 s-1]; may be unallocated
    amrex::FArrayBox Ah_bg_xx_fab;
    amrex::Array4<const amrex::Real> Ah_bg_xx{};
    if (captured.logical("_metrics_Ah_bg_xx_allocated")) {
        Ah_bg_xx_fab = captured.fab_device("_metrics_Ah_bg_xx_before");
        Ah_bg_xx = Ah_bg_xx_fab.const_array();
    }

    // The maximum permitted Laplacian viscosity [L2 T-1 ~> m2 s-1]; may be unallocated
    amrex::FArrayBox Kh_Max_xx_fab;
    amrex::Array4<const amrex::Real> Kh_Max_xx{};
    if (captured.logical("_metrics_Kh_Max_xx_allocated")) {
        Kh_Max_xx_fab = captured.fab_device("_metrics_Kh_Max_xx_before");
        Kh_Max_xx = Kh_Max_xx_fab.const_array();
    }

    // Factor n1*n2 in the anisotropic direction tensor at h-points [nondim]; may be unallocated
    amrex::FArrayBox n1n2_h_fab;
    amrex::Array4<const amrex::Real> n1n2_h{};
    if (captured.logical("_metrics_n1n2_h_allocated")) {
        n1n2_h_fab = captured.fab_device("_metrics_n1n2_h_before");
        n1n2_h = n1n2_h_fab.const_array();
    }

    // Factor n1**2-n2**2 in the anisotropic direction tensor at h-points [nondim]; may be unallocated
    amrex::FArrayBox n1n1_m_n2n2_h_fab;
    amrex::Array4<const amrex::Real> n1n1_m_n2n2_h{};
    if (captured.logical("_metrics_n1n1_m_n2n2_h_allocated")) {
        n1n1_m_n2n2_h_fab = captured.fab_device("_metrics_n1n1_m_n2n2_h_before");
        n1n1_m_n2n2_h = n1n1_m_n2n2_h_fab.const_array();
    }

    // Harmonic mean of the squares of the grid [L2 ~> m2]; may be unallocated
    amrex::FArrayBox grid_sp_h2_fab;
    amrex::Array4<const amrex::Real> grid_sp_h2{};
    if (captured.logical("_metrics_grid_sp_h2_allocated")) {
        grid_sp_h2_fab = captured.fab_device("_metrics_grid_sp_h2_before");
        grid_sp_h2 = grid_sp_h2_fab.const_array();
    }

    // Laplacian metric-dependent constants [L2 ~> m2]; may be unallocated
    amrex::FArrayBox Laplac2_const_xx_fab;
    amrex::Array4<const amrex::Real> Laplac2_const_xx{};
    if (captured.logical("_metrics_Laplac2_const_xx_allocated")) {
        Laplac2_const_xx_fab = captured.fab_device("_metrics_Laplac2_const_xx_before");
        Laplac2_const_xx = Laplac2_const_xx_fab.const_array();
    }

    // Biharmonic metric-dependent constants [L6 ~> m6]; may be unallocated
    amrex::FArrayBox Biharm6_const_xx_fab;
    amrex::Array4<const amrex::Real> Biharm6_const_xx{};
    if (captured.logical("_metrics_Biharm6_const_xx_allocated")) {
        Biharm6_const_xx_fab = captured.fab_device("_metrics_Biharm6_const_xx_before");
        Biharm6_const_xx = Biharm6_const_xx_fab.const_array();
    }

    // Laplacian metric-dependent constants [L3 ~> m3]; may be unallocated
    amrex::FArrayBox Laplac3_const_xx_fab;
    amrex::Array4<const amrex::Real> Laplac3_const_xx{};
    if (captured.logical("_metrics_Laplac3_const_xx_allocated")) {
        Laplac3_const_xx_fab = captured.fab_device("_metrics_Laplac3_const_xx_before");
        Laplac3_const_xx = Laplac3_const_xx_fab.const_array();
    }

    // Biharmonic metric-dependent constants [L4 ~> m4]; may be unallocated
    amrex::FArrayBox Biharm_const_xx_fab;
    amrex::Array4<const amrex::Real> Biharm_const_xx{};
    if (captured.logical("_metrics_Biharm_const_xx_allocated")) {
        Biharm_const_xx_fab = captured.fab_device("_metrics_Biharm_const_xx_before");
        Biharm_const_xx = Biharm_const_xx_fab.const_array();
    }

    // Biharmonic metric-dependent constants [T L4 ~> s m4]; may be unallocated
    amrex::FArrayBox Biharm_const2_xx_fab;
    amrex::Array4<const amrex::Real> Biharm_const2_xx{};
    if (captured.logical("_metrics_Biharm_const2_xx_allocated")) {
        Biharm_const2_xx_fab = captured.fab_device("_metrics_Biharm_const2_xx_before");
        Biharm_const2_xx = Biharm_const2_xx_fab.const_array();
    }

    const int    kstart        = captured.integer("_kstart");
    const double h_neglect     = captured.real64("_h_neglect");
    const double inv_PI3       = captured.real64("_inv_PI3");
    const double inv_PI6       = captured.real64("_inv_PI6");
    const bool   use_kh_struct = captured.logical("_use_kh_struct");
    const bool   use_MEKE_Ku   = captured.logical("_use_MEKE_Ku");
    const bool   rescale_Kh    = captured.logical("_rescale_Kh");

    const auto u                        = captured.fab_device("_u");
    const auto v                        = captured.fab_device("_v");
    const auto h                        = captured.fab_device("_h");
    const auto h_u                      = captured.fab_device("_h_u_before");
    const auto h_v                      = captured.fab_device("_h_v_before");
    const auto sh_xx                    = captured.fab_device("_sh_xx_before");
    const auto grad_vort_mag_h          = captured.fab_device("_grad_vort_mag_h_before");
    const auto grad_vort_mag_h_2d       = captured.fab_device("_grad_vort_mag_h_2d_before");
    const auto grad_div_mag_h           = captured.fab_device("_grad_div_mag_h_before");
    const auto dudx                     = captured.fab_device("_dudx_before");
    const auto dvdy                     = captured.fab_device("_dvdy_before");
    const auto sh_xy                    = captured.fab_device("_sh_xy_before");
    const auto Del2vort_q               = captured.fab_device("_Del2vort_q_before");
    auto       Kh_h                     = captured.fab_device("_Kh_h_before");
    auto       grid_Re_Kh               = captured.fab_device("_grid_Re_Kh_before");
    auto       div_xx_h                 = captured.fab_device("_div_xx_h_before");
    auto       sh_xx_h                  = captured.fab_device("_sh_xx_h_before");
    auto       str_xx                   = captured.fab_device("_str_xx_before");
    auto       Ah                       = captured.fab_device("_Ah_before");
    auto       Kh                       = captured.fab_device("_Kh_before");
    auto       Shear_mag                = captured.fab_device("_Shear_mag_before");
    auto       vert_vort_mag            = captured.fab_device("_vert_vort_mag_before");
    auto       hrat_min                 = captured.fab_device("_hrat_min_before");
    auto       visc_bound_rem           = captured.fab_device("_visc_bound_rem_before");
    const auto Kh_h_after               = captured.fab_host("_Kh_h_after");
    const auto grid_Re_Kh_after         = captured.fab_host("_grid_Re_Kh_after");
    const auto div_xx_h_after           = captured.fab_host("_div_xx_h_after");
    const auto sh_xx_h_after            = captured.fab_host("_sh_xx_h_after");
    const auto str_xx_after             = captured.fab_host("_str_xx_after");
    const auto Ah_after                 = captured.fab_host("_Ah_after");
    const auto Kh_after                 = captured.fab_host("_Kh_after");
    const auto Shear_mag_after          = captured.fab_host("_Shear_mag_after");
    const auto vert_vort_mag_after      = captured.fab_host("_vert_vort_mag_after");
    const auto hrat_min_after           = captured.fab_host("_hrat_min_after");
    const auto visc_bound_rem_after     = captured.fab_host("_visc_bound_rem_after");

    MOM::hor_visc_h_coeffs_a(bxH, bxKh, bxHstr, CS, MEKE_Ku, Res_fn_h, BS_struct, Kh_bg_xx, Ah_bg_xx, Kh_Max_xx, n1n2_h,
                             n1n1_m_n2n2_h, grid_sp_h2, Laplac2_const_xx, Biharm6_const_xx, Laplac3_const_xx,
                             Biharm_const_xx, Biharm_const2_xx, kstart, h_neglect, inv_PI3, inv_PI6, use_kh_struct,
                             use_MEKE_Ku, rescale_Kh, u.const_array(), v.const_array(), h.const_array(),
                             Kh_h.array(), grid_Re_Kh.array(), div_xx_h.array(), sh_xx_h.array(), h_u.const_array(),
                             h_v.const_array(), sh_xx.const_array(), str_xx.array(), grad_vort_mag_h.const_array(),
                             grad_vort_mag_h_2d.const_array(), grad_div_mag_h.const_array(), dudx.const_array(),
                             dvdy.const_array(), sh_xy.const_array(), Del2vort_q.const_array(), Ah.array(),
                             Kh.array(), Shear_mag.array(), vert_vort_mag.array(), hrat_min.array(),
                             visc_bound_rem.array());
    amrex::Gpu::synchronize();

    expect_arrays_equal(Kh_h_after,               to_host_fab(Kh_h),               "Kh_h");
    expect_arrays_equal(grid_Re_Kh_after,         to_host_fab(grid_Re_Kh),         "grid_Re_Kh");
    expect_arrays_equal(div_xx_h_after,           to_host_fab(div_xx_h),           "div_xx_h");
    expect_arrays_equal(sh_xx_h_after,            to_host_fab(sh_xx_h),            "sh_xx_h");
    expect_arrays_equal(str_xx_after,             to_host_fab(str_xx),             "str_xx");
    expect_arrays_equal(Ah_after,                 to_host_fab(Ah),                 "Ah");
    expect_arrays_equal(Kh_after,                 to_host_fab(Kh),                 "Kh");
    expect_arrays_equal(Shear_mag_after,          to_host_fab(Shear_mag),          "Shear_mag");
    expect_arrays_equal(vert_vort_mag_after,      to_host_fab(vert_vort_mag),      "vert_vort_mag");
    expect_arrays_equal(hrat_min_after,           to_host_fab(hrat_min),           "hrat_min");
    expect_arrays_equal(visc_bound_rem_after,     to_host_fab(visc_bound_rem),     "visc_bound_rem");
}

// -------------------------------------------------------------------------
// Field mapping: MOM6 src/parameterizations/lateral/MOM_hor_visc.F90:4767-4866 (hor_visc_h_coeffs_b capture block).
TEST(HorViscHCoeffsB, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "hor_visc_h_coeffs_b");

    const auto bxH    = captured.box("_bxH");
    const auto bxQ    = captured.box("_bxQ");
    const auto bxKh   = captured.box("_bxKh");
    const auto bxHstr = captured.box("_bxHstr");
    const auto bxQstr = captured.box("_bxQstr");

    hv_opts_CS_C CS{};
    CS.Laplacian              = captured.logical("_CS_Laplacian");
    CS.biharmonic             = captured.logical("_CS_biharmonic");
    CS.debug                  = captured.logical("_CS_debug");
    CS.no_slip                = captured.logical("_CS_no_slip");
    CS.bound_Kh               = captured.logical("_CS_bound_Kh");
    CS.EY24_EBT_BS            = captured.logical("_CS_EY24_EBT_BS");
    CS.bound_Ah               = captured.logical("_CS_bound_Ah");
    CS.Re_Ah                  = captured.real64("_CS_Re_Ah");
    CS.KS_coef                = captured.real64("_CS_KS_coef");
    CS.backscatter_underbound = captured.logical("_CS_backscatter_underbound");
    CS.Smagorinsky_Kh         = captured.logical("_CS_Smagorinsky_Kh");
    CS.Smagorinsky_Ah         = captured.logical("_CS_Smagorinsky_Ah");
    CS.Leith_Kh               = captured.logical("_CS_Leith_Kh");
    CS.Modified_Leith         = captured.logical("_CS_Modified_Leith");
    CS.use_beta_in_Leith      = captured.logical("_CS_use_beta_in_Leith");
    CS.Leith_Ah               = captured.logical("_CS_Leith_Ah");
    CS.use_Leithy             = captured.logical("_CS_use_Leithy");
    CS.c_K                    = captured.real64("_CS_c_K");
    CS.smooth_Ah              = captured.logical("_CS_smooth_Ah");
    CS.use_QG_Leith_visc      = captured.logical("_CS_use_QG_Leith_visc");
    CS.bound_Coriolis         = captured.logical("_CS_bound_Coriolis");
    CS.Kh_bg_min              = captured.real64("_CS_Kh_bg_min");
    CS.FrictWork_bug          = captured.logical("_CS_FrictWork_bug");
    CS.OBC_strain_bug         = captured.logical("_CS_OBC_strain_bug");
    CS.use_land_mask          = captured.logical("_CS_use_land_mask");
    CS.anisotropic            = captured.logical("_CS_anisotropic");
    CS.add_LES_viscosity      = captured.logical("_CS_add_LES_viscosity");
    CS.Kh_aniso               = captured.real64("_CS_Kh_aniso");
    CS.res_scale_MEKE         = captured.logical("_CS_res_scale_MEKE");
    CS.use_GME                = captured.logical("_CS_use_GME");
    CS.answer_date            = captured.integer("_CS_answer_date");
    CS.GME_h0                 = captured.real64("_CS_GME_h0");
    CS.GME_efficiency         = captured.real64("_CS_GME_efficiency");
    CS.GME_limiter            = captured.real64("_CS_GME_limiter");
    CS.min_grid_Kh            = captured.real64("_CS_min_grid_Kh");
    CS.min_grid_Ah            = captured.real64("_CS_min_grid_Ah");
    CS.use_cont_thick         = captured.logical("_CS_use_cont_thick");
    CS.use_ZB2020             = captured.logical("_CS_use_ZB2020");
    CS.use_circulation        = captured.logical("_CS_use_circulation");
    CS.do_grid_Re_Ah          = captured.logical("_CS_do_grid_Re_Ah");
    CS.do_grid_Re_Kh          = captured.logical("_CS_do_grid_Re_Kh");
    CS.do_Ah_h                = captured.logical("_CS_do_Ah_h");
    CS.do_Ah_q                = captured.logical("_CS_do_Ah_q");
    CS.do_Kh_h                = captured.logical("_CS_do_Kh_h");
    CS.do_Kh_q                = captured.logical("_CS_do_Kh_q");
    CS.do_GME_coeff_h         = captured.logical("_CS_do_GME_coeff_h");
    CS.do_GME_coeff_q         = captured.logical("_CS_do_GME_coeff_q");
    CS.do_vort_xy_q           = captured.logical("_CS_do_vort_xy_q");
    CS.do_div_xx_h            = captured.logical("_CS_do_div_xx_h");
    CS.do_sh_xy_q             = captured.logical("_CS_do_sh_xy_q");
    CS.do_sh_xx_h             = captured.logical("_CS_do_sh_xx_h");
    CS.do_FrictWork           = captured.logical("_CS_do_FrictWork");
    CS.do_FrictWorkIntz       = captured.logical("_CS_do_FrictWorkIntz");
    CS.do_FrictWork_bh        = captured.logical("_CS_do_FrictWork_bh");
    CS.do_FrictWorkIntz_bh    = captured.logical("_CS_do_FrictWorkIntz_bh");
    CS.do_normstress          = captured.logical("_CS_do_normstress");
    CS.do_shearstress         = captured.logical("_CS_do_shearstress");
    CS.do_BS_coeff_h          = captured.logical("_CS_do_BS_coeff_h");
    CS.do_BS_coeff_q          = captured.logical("_CS_do_BS_coeff_q");

    // The MEKE-derived lateral biharmonic viscosity coefficient [L4 T-1 ~> m4 s-1]; may be unallocated
    amrex::FArrayBox MEKE_Au_fab;
    amrex::Array4<const amrex::Real> MEKE_Au{};
    if (captured.logical("_MEKE_Au_allocated")) {
        MEKE_Au_fab = captured.fab_device("_MEKE_Au_before");
        MEKE_Au = MEKE_Au_fab.const_array();
    }

    // The MEKE-derived lateral viscosity coefficient [L2 T-1 ~> m2 s-1]; may be unallocated
    amrex::FArrayBox MEKE_Ku_fab;
    amrex::Array4<const amrex::Real> MEKE_Ku{};
    if (captured.logical("_MEKE_Ku_allocated")) {
        MEKE_Ku_fab = captured.fab_device("_MEKE_Ku_before");
        MEKE_Ku = MEKE_Ku_fab.const_array();
    }

    // Vertical structure function used in backscatter [nondim]; may be unallocated
    amrex::FArrayBox BS_struct_fab;
    amrex::Array4<const amrex::Real> BS_struct{};
    if (captured.logical("_VarMix_BS_struct_allocated")) {
        BS_struct_fab = captured.fab_device("_VarMix_BS_struct_before");
        BS_struct = BS_struct_fab.const_array();
    }

    // The amount by which stresses through h points are reduced due to partial barriers [nondim]; may be unallocated
    amrex::FArrayBox reduction_xx_fab;
    amrex::Array4<const amrex::Real> reduction_xx{};
    if (captured.logical("_metrics_reduction_xx_allocated")) {
        reduction_xx_fab = captured.fab_device("_metrics_reduction_xx_before");
        reduction_xx = reduction_xx_fab.const_array();
    }

    // The maximum permitted biharmonic viscosity [L4 T-1 ~> m4 s-1]; may be unallocated
    amrex::FArrayBox Ah_Max_xx_fab;
    amrex::Array4<const amrex::Real> Ah_Max_xx{};
    if (captured.logical("_metrics_Ah_Max_xx_allocated")) {
        Ah_Max_xx_fab = captured.fab_device("_metrics_Ah_Max_xx_before");
        Ah_Max_xx = Ah_Max_xx_fab.const_array();
    }

    // The maximum permitted biharmonic viscosity for the kill switch [L4 T-1 ~> m4 s-1]; may be unallocated
    amrex::FArrayBox Ah_Max_xx_KS_fab;
    amrex::Array4<const amrex::Real> Ah_Max_xx_KS{};
    if (captured.logical("_metrics_Ah_Max_xx_KS_allocated")) {
        Ah_Max_xx_KS_fab = captured.fab_device("_metrics_Ah_Max_xx_KS_before");
        Ah_Max_xx_KS = Ah_Max_xx_KS_fab.const_array();
    }

    // Harmonic mean of the squares of the grid^(3/2) [L3 ~> m3]; may be unallocated
    amrex::FArrayBox grid_sp_h3_fab;
    amrex::Array4<const amrex::Real> grid_sp_h3{};
    if (captured.logical("_metrics_grid_sp_h3_allocated")) {
        grid_sp_h3_fab = captured.fab_device("_metrics_grid_sp_h3_before");
        grid_sp_h3 = grid_sp_h3_fab.const_array();
    }

    // Pre-calculated dx/dy at h points [nondim]; may be unallocated
    amrex::FArrayBox dx_dyT_fab;
    amrex::Array4<const amrex::Real> dx_dyT{};
    if (captured.logical("_metrics_dx_dyT_allocated")) {
        dx_dyT_fab = captured.fab_device("_metrics_dx_dyT_before");
        dx_dyT = dx_dyT_fab.const_array();
    }

    // Pre-calculated dy/dx at h points [nondim]; may be unallocated
    amrex::FArrayBox dy_dxT_fab;
    amrex::Array4<const amrex::Real> dy_dxT{};
    if (captured.logical("_metrics_dy_dxT_allocated")) {
        dy_dxT_fab = captured.fab_device("_metrics_dy_dxT_before");
        dy_dxT = dy_dxT_fab.const_array();
    }

    // Pre-calculated dx/dy at q points [nondim]; may be unallocated
    amrex::FArrayBox dx_dyBu_fab;
    amrex::Array4<const amrex::Real> dx_dyBu{};
    if (captured.logical("_metrics_dx_dyBu_allocated")) {
        dx_dyBu_fab = captured.fab_device("_metrics_dx_dyBu_before");
        dx_dyBu = dx_dyBu_fab.const_array();
    }

    // Pre-calculated dy/dx at q points [nondim]; may be unallocated
    amrex::FArrayBox dy_dxBu_fab;
    amrex::Array4<const amrex::Real> dy_dxBu{};
    if (captured.logical("_metrics_dy_dxBu_allocated")) {
        dy_dxBu_fab = captured.fab_device("_metrics_dy_dxBu_before");
        dy_dxBu = dy_dxBu_fab.const_array();
    }

    // Biharmonic metric-dependent constants [L3 ~> m3]; may be unallocated
    amrex::FArrayBox Re_Ah_const_xx_fab;
    amrex::Array4<const amrex::Real> Re_Ah_const_xx{};
    if (captured.logical("_metrics_Re_Ah_const_xx_allocated")) {
        Re_Ah_const_xx_fab = captured.fab_device("_metrics_Re_Ah_const_xx_before");
        Re_Ah_const_xx = Re_Ah_const_xx_fab.const_array();
    }

    const auto IdyCu    = captured.fab_device("_IdyCu");
    const auto IdxCv    = captured.fab_device("_IdxCv");
    const auto IdyCv    = captured.fab_device("_IdyCv");
    const auto IdxCu    = captured.fab_device("_IdxCu");
    const auto mask2dBu = captured.fab_device("_mask2dBu");
    const auto mask2dCu = captured.fab_device("_mask2dCu");
    const auto mask2dCv = captured.fab_device("_mask2dCv");

    const int    kstart        = captured.integer("_kstart");
    const double h_neglect     = captured.real64("_h_neglect");
    const double h_neglect3    = captured.real64("_h_neglect3");
    const bool   apply_OBC     = captured.logical("_apply_OBC");
    const bool   use_kh_struct = captured.logical("_use_kh_struct");
    const bool   use_MEKE_Au   = captured.logical("_use_MEKE_Au");

    const auto u                       = captured.fab_device("_u");
    const auto v                       = captured.fab_device("_v");
    const auto h                       = captured.fab_device("_h");
    const auto Del2u                   = captured.fab_device("_Del2u_before");
    const auto h_u                     = captured.fab_device("_h_u_before");
    const auto Del2v                   = captured.fab_device("_Del2v_before");
    const auto h_v                     = captured.fab_device("_h_v_before");
    const auto sh_xx                   = captured.fab_device("_sh_xx_before");
    const auto sh_xx_smooth            = captured.fab_device("_sh_xx_smooth_before");
    const auto m_leithy                = captured.fab_device("_m_leithy_before");
    const auto sh_xy                   = captured.fab_device("_sh_xy_before");
    const auto visc_bound_rem          = captured.fab_device("_visc_bound_rem_before");
    auto       Ah_h                    = captured.fab_device("_Ah_h_before");
    auto       Kh_h                    = captured.fab_device("_Kh_h_before");
    auto       grid_Re_Ah              = captured.fab_device("_grid_Re_Ah_before");
    auto       visc_limit_h            = captured.fab_device("_visc_limit_h_before");
    auto       visc_limit_h_flag       = captured.fab_device("_visc_limit_h_flag_before");
    auto       visc_limit_h_frac       = captured.fab_device("_visc_limit_h_frac_before");
    auto       BS_coeff_h              = captured.fab_device("_BS_coeff_h_before");
    auto       str_xx                  = captured.fab_device("_str_xx_before");
    auto       bhstr_xx                = captured.fab_device("_bhstr_xx_before");
    auto       dDel2vdx                = captured.fab_device("_dDel2vdx_before");
    auto       dDel2udy                = captured.fab_device("_dDel2udy_before");
    auto       hq                      = captured.fab_device("_hq_before");
    auto       Ah                      = captured.fab_device("_Ah_before");
    auto       Kh                      = captured.fab_device("_Kh_before");
    auto       Shear_mag               = captured.fab_device("_Shear_mag_before");
    auto       hrat_min                = captured.fab_device("_hrat_min_before");
    const auto Ah_h_after              = captured.fab_host("_Ah_h_after");
    const auto Kh_h_after              = captured.fab_host("_Kh_h_after");
    const auto grid_Re_Ah_after        = captured.fab_host("_grid_Re_Ah_after");
    const auto visc_limit_h_after      = captured.fab_host("_visc_limit_h_after");
    const auto visc_limit_h_flag_after = captured.fab_host("_visc_limit_h_flag_after");
    const auto visc_limit_h_frac_after = captured.fab_host("_visc_limit_h_frac_after");
    const auto BS_coeff_h_after        = captured.fab_host("_BS_coeff_h_after");
    const auto str_xx_after            = captured.fab_host("_str_xx_after");
    const auto bhstr_xx_after          = captured.fab_host("_bhstr_xx_after");
    const auto dDel2vdx_after          = captured.fab_host("_dDel2vdx_after");
    const auto dDel2udy_after          = captured.fab_host("_dDel2udy_after");
    const auto hq_after                = captured.fab_host("_hq_after");
    const auto Ah_after                = captured.fab_host("_Ah_after");
    const auto Kh_after                = captured.fab_host("_Kh_after");
    const auto Shear_mag_after         = captured.fab_host("_Shear_mag_after");
    const auto hrat_min_after          = captured.fab_host("_hrat_min_after");

    MOM::hor_visc_h_coeffs_b(bxH, bxQ, bxKh, bxHstr, bxQstr, CS, MEKE_Au, MEKE_Ku, BS_struct, IdyCu.const_array(),
                             IdxCv.const_array(), IdyCv.const_array(), IdxCu.const_array(), mask2dBu.const_array(),
                             mask2dCu.const_array(), mask2dCv.const_array(), reduction_xx, Ah_Max_xx, Ah_Max_xx_KS,
                             grid_sp_h3, dx_dyT, dy_dxT, dx_dyBu, dy_dxBu, Re_Ah_const_xx, nullptr, kstart,
                             h_neglect, h_neglect3, apply_OBC, use_kh_struct, use_MEKE_Au, u.const_array(),
                             v.const_array(), h.const_array(), Ah_h.array(), Kh_h.array(), grid_Re_Ah.array(),
                             visc_limit_h.array(), visc_limit_h_flag.array(), visc_limit_h_frac.array(),
                             BS_coeff_h.array(), Del2u.const_array(), h_u.const_array(), Del2v.const_array(),
                             h_v.const_array(), sh_xx.const_array(), sh_xx_smooth.const_array(), str_xx.array(),
                             bhstr_xx.array(), m_leithy.const_array(), dDel2vdx.array(), dDel2udy.array(),
                             sh_xy.const_array(), hq.array(), Ah.array(), Kh.array(), Shear_mag.array(),
                             hrat_min.array(), visc_bound_rem.const_array());
    amrex::Gpu::synchronize();

    expect_arrays_equal(Ah_h_after,              to_host_fab(Ah_h),              "Ah_h");
    expect_arrays_equal(Kh_h_after,              to_host_fab(Kh_h),              "Kh_h");
    expect_arrays_equal(grid_Re_Ah_after,        to_host_fab(grid_Re_Ah),        "grid_Re_Ah");
    expect_arrays_equal(visc_limit_h_after,      to_host_fab(visc_limit_h),      "visc_limit_h");
    expect_arrays_equal(visc_limit_h_flag_after, to_host_fab(visc_limit_h_flag), "visc_limit_h_flag");
    expect_arrays_equal(visc_limit_h_frac_after, to_host_fab(visc_limit_h_frac), "visc_limit_h_frac");
    expect_arrays_equal(BS_coeff_h_after,        to_host_fab(BS_coeff_h),        "BS_coeff_h");
    expect_arrays_equal(str_xx_after,            to_host_fab(str_xx),            "str_xx");
    expect_arrays_equal(bhstr_xx_after,          to_host_fab(bhstr_xx),          "bhstr_xx");
    expect_arrays_equal(dDel2vdx_after,          to_host_fab(dDel2vdx),          "dDel2vdx");
    expect_arrays_equal(dDel2udy_after,          to_host_fab(dDel2udy),          "dDel2udy");
    expect_arrays_equal(hq_after,                to_host_fab(hq),                "hq");
    expect_arrays_equal(Ah_after,                to_host_fab(Ah),                "Ah");
    expect_arrays_equal(Kh_after,                to_host_fab(Kh),                "Kh");
    expect_arrays_equal(Shear_mag_after,         to_host_fab(Shear_mag),         "Shear_mag");
    expect_arrays_equal(hrat_min_after,          to_host_fab(hrat_min),          "hrat_min");
}

// -------------------------------------------------------------------------
// Field mapping: MOM6 src/parameterizations/lateral/MOM_hor_visc.F90:6452-6549 (hor_visc_stress_div capture block).
TEST(HorViscStressDiv, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "hor_visc_stress_div");

    const auto bxH    = captured.box("_bxH");
    const auto bxU    = captured.box("_bxU");
    const auto bxV    = captured.box("_bxV");
    const auto bxHstr = captured.box("_bxHstr");
    const auto bxQstr = captured.box("_bxQstr");

    hv_opts_CS_C CS{};
    CS.Laplacian              = captured.logical("_CS_Laplacian");
    CS.biharmonic             = captured.logical("_CS_biharmonic");
    CS.debug                  = captured.logical("_CS_debug");
    CS.no_slip                = captured.logical("_CS_no_slip");
    CS.bound_Kh               = captured.logical("_CS_bound_Kh");
    CS.EY24_EBT_BS            = captured.logical("_CS_EY24_EBT_BS");
    CS.bound_Ah               = captured.logical("_CS_bound_Ah");
    CS.Re_Ah                  = captured.real64("_CS_Re_Ah");
    CS.KS_coef                = captured.real64("_CS_KS_coef");
    CS.backscatter_underbound = captured.logical("_CS_backscatter_underbound");
    CS.Smagorinsky_Kh         = captured.logical("_CS_Smagorinsky_Kh");
    CS.Smagorinsky_Ah         = captured.logical("_CS_Smagorinsky_Ah");
    CS.Leith_Kh               = captured.logical("_CS_Leith_Kh");
    CS.Modified_Leith         = captured.logical("_CS_Modified_Leith");
    CS.use_beta_in_Leith      = captured.logical("_CS_use_beta_in_Leith");
    CS.Leith_Ah               = captured.logical("_CS_Leith_Ah");
    CS.use_Leithy             = captured.logical("_CS_use_Leithy");
    CS.c_K                    = captured.real64("_CS_c_K");
    CS.smooth_Ah              = captured.logical("_CS_smooth_Ah");
    CS.use_QG_Leith_visc      = captured.logical("_CS_use_QG_Leith_visc");
    CS.bound_Coriolis         = captured.logical("_CS_bound_Coriolis");
    CS.Kh_bg_min              = captured.real64("_CS_Kh_bg_min");
    CS.FrictWork_bug          = captured.logical("_CS_FrictWork_bug");
    CS.OBC_strain_bug         = captured.logical("_CS_OBC_strain_bug");
    CS.use_land_mask          = captured.logical("_CS_use_land_mask");
    CS.anisotropic            = captured.logical("_CS_anisotropic");
    CS.add_LES_viscosity      = captured.logical("_CS_add_LES_viscosity");
    CS.Kh_aniso               = captured.real64("_CS_Kh_aniso");
    CS.res_scale_MEKE         = captured.logical("_CS_res_scale_MEKE");
    CS.use_GME                = captured.logical("_CS_use_GME");
    CS.answer_date            = captured.integer("_CS_answer_date");
    CS.GME_h0                 = captured.real64("_CS_GME_h0");
    CS.GME_efficiency         = captured.real64("_CS_GME_efficiency");
    CS.GME_limiter            = captured.real64("_CS_GME_limiter");
    CS.min_grid_Kh            = captured.real64("_CS_min_grid_Kh");
    CS.min_grid_Ah            = captured.real64("_CS_min_grid_Ah");
    CS.use_cont_thick         = captured.logical("_CS_use_cont_thick");
    CS.use_ZB2020             = captured.logical("_CS_use_ZB2020");
    CS.use_circulation        = captured.logical("_CS_use_circulation");
    CS.do_grid_Re_Ah          = captured.logical("_CS_do_grid_Re_Ah");
    CS.do_grid_Re_Kh          = captured.logical("_CS_do_grid_Re_Kh");
    CS.do_Ah_h                = captured.logical("_CS_do_Ah_h");
    CS.do_Ah_q                = captured.logical("_CS_do_Ah_q");
    CS.do_Kh_h                = captured.logical("_CS_do_Kh_h");
    CS.do_Kh_q                = captured.logical("_CS_do_Kh_q");
    CS.do_GME_coeff_h         = captured.logical("_CS_do_GME_coeff_h");
    CS.do_GME_coeff_q         = captured.logical("_CS_do_GME_coeff_q");
    CS.do_vort_xy_q           = captured.logical("_CS_do_vort_xy_q");
    CS.do_div_xx_h            = captured.logical("_CS_do_div_xx_h");
    CS.do_sh_xy_q             = captured.logical("_CS_do_sh_xy_q");
    CS.do_sh_xx_h             = captured.logical("_CS_do_sh_xx_h");
    CS.do_FrictWork           = captured.logical("_CS_do_FrictWork");
    CS.do_FrictWorkIntz       = captured.logical("_CS_do_FrictWorkIntz");
    CS.do_FrictWork_bh        = captured.logical("_CS_do_FrictWork_bh");
    CS.do_FrictWorkIntz_bh    = captured.logical("_CS_do_FrictWorkIntz_bh");
    CS.do_normstress          = captured.logical("_CS_do_normstress");
    CS.do_shearstress         = captured.logical("_CS_do_shearstress");
    CS.do_BS_coeff_h          = captured.logical("_CS_do_BS_coeff_h");
    CS.do_BS_coeff_q          = captured.logical("_CS_do_BS_coeff_q");

    // MEKE source from lateral friction in the momentum equations [R Z L2 T-3 ~> W m-2]; may be unallocated
    const bool MEKE_mom_src_allocated = captured.logical("_MEKE_mom_src_allocated");
    amrex::FArrayBox MEKE_mom_src_fab;
    amrex::Array4<amrex::Real> MEKE_mom_src{};
    if (MEKE_mom_src_allocated) {
        MEKE_mom_src_fab = captured.fab_device("_MEKE_mom_src_before");
        MEKE_mom_src = MEKE_mom_src_fab.array();
    }

    // MEKE source from the biharmonic part of the lateral friction in the momentum equations [R Z L2 T-3 ~> W
    // m-2]; may be unallocated
    const bool MEKE_mom_src_bh_allocated = captured.logical("_MEKE_mom_src_bh_allocated");
    amrex::FArrayBox MEKE_mom_src_bh_fab;
    amrex::Array4<amrex::Real> MEKE_mom_src_bh{};
    if (MEKE_mom_src_bh_allocated) {
        MEKE_mom_src_bh_fab = captured.fab_device("_MEKE_mom_src_bh_before");
        MEKE_mom_src_bh = MEKE_mom_src_bh_fab.array();
    }

    // MEKE sink from GME backscatter in the momentum equations [R Z L2 T-3 ~> W m-2]; may be unallocated
    const bool MEKE_GME_snk_allocated = captured.logical("_MEKE_GME_snk_allocated");
    amrex::FArrayBox MEKE_GME_snk_fab;
    amrex::Array4<amrex::Real> MEKE_GME_snk{};
    if (MEKE_GME_snk_allocated) {
        MEKE_GME_snk_fab = captured.fab_device("_MEKE_GME_snk_before");
        MEKE_GME_snk = MEKE_GME_snk_fab.array();
    }

    // Dissipation rate used to set amplitude of SKEBS [L2 T-3 ~> m2 s-3]; may be unallocated
    const bool skeb_diss_allocated = captured.logical("_STOCH_skeb_diss_allocated");
    amrex::FArrayBox skeb_diss_fab;
    amrex::Array4<amrex::Real> skeb_diss{};
    if (skeb_diss_allocated) {
        skeb_diss_fab = captured.fab_device("_STOCH_skeb_diss_before");
        skeb_diss = skeb_diss_fab.array();
    }

    const double backscatter_Ro_pow = captured.real64("_MEKE_backscatter_Ro_pow");
    const double backscatter_Ro_c   = captured.real64("_MEKE_backscatter_Ro_c");
    const double skeb_frict_coef    = captured.real64("_STOCH_skeb_frict_coef");

    const auto mask2dBu   = captured.fab_device("_mask2dBu");
    const auto IdxCu      = captured.fab_device("_IdxCu");
    const auto IdyCu      = captured.fab_device("_IdyCu");
    const auto IareaCu    = captured.fab_device("_IareaCu");
    const auto IdyCv      = captured.fab_device("_IdyCv");
    const auto IdxCv      = captured.fab_device("_IdxCv");
    const auto IareaCv    = captured.fab_device("_IareaCv");
    const auto IdxT       = captured.fab_device("_IdxT");
    const auto IdyT       = captured.fab_device("_IdyT");
    const auto IdyBu      = captured.fab_device("_IdyBu");
    const auto IdxBu      = captured.fab_device("_IdxBu");
    const auto IareaT     = captured.fab_device("_IareaT");
    const auto dxCu       = captured.fab_device("_dxCu");
    const auto dyCv       = captured.fab_device("_dyCv");
    const auto CoriolisBu = captured.fab_device("_CoriolisBu");

    // The amount by which stresses through h points are reduced due to partial barriers [nondim]; may be unallocated
    amrex::FArrayBox reduction_xx_fab;
    amrex::Array4<const amrex::Real> reduction_xx{};
    if (captured.logical("_metrics_reduction_xx_allocated")) {
        reduction_xx_fab = captured.fab_device("_metrics_reduction_xx_before");
        reduction_xx = reduction_xx_fab.const_array();
    }

    // The amount by which stresses through q points are reduced due to partial barriers [nondim]; may be unallocated
    amrex::FArrayBox reduction_xy_fab;
    amrex::Array4<const amrex::Real> reduction_xy{};
    if (captured.logical("_metrics_reduction_xy_allocated")) {
        reduction_xy_fab = captured.fab_device("_metrics_reduction_xy_before");
        reduction_xy = reduction_xy_fab.const_array();
    }

    // Pre-calculated dx^2 at h points [L2 ~> m2]; may be unallocated
    amrex::FArrayBox dx2h_fab;
    amrex::Array4<const amrex::Real> dx2h{};
    if (captured.logical("_metrics_dx2h_allocated")) {
        dx2h_fab = captured.fab_device("_metrics_dx2h_before");
        dx2h = dx2h_fab.const_array();
    }

    // Pre-calculated dy^2 at h points [L2 ~> m2]; may be unallocated
    amrex::FArrayBox dy2h_fab;
    amrex::Array4<const amrex::Real> dy2h{};
    if (captured.logical("_metrics_dy2h_allocated")) {
        dy2h_fab = captured.fab_device("_metrics_dy2h_before");
        dy2h = dy2h_fab.const_array();
    }

    // Pre-calculated dx^2 at q points [L2 ~> m2]; may be unallocated
    amrex::FArrayBox dx2q_fab;
    amrex::Array4<const amrex::Real> dx2q{};
    if (captured.logical("_metrics_dx2q_allocated")) {
        dx2q_fab = captured.fab_device("_metrics_dx2q_before");
        dx2q = dx2q_fab.const_array();
    }

    // Pre-calculated dy^2 at q points [L2 ~> m2]; may be unallocated
    amrex::FArrayBox dy2q_fab;
    amrex::Array4<const amrex::Real> dy2q{};
    if (captured.logical("_metrics_dy2q_allocated")) {
        dy2q_fab = captured.fab_device("_metrics_dy2q_before");
        dy2q = dy2q_fab.const_array();
    }

    const int    kstart            = captured.integer("_kstart");
    const double h_neglect         = captured.real64("_h_neglect");
    const bool   apply_OBC         = captured.logical("_apply_OBC");
    const bool   find_FrictWork    = captured.logical("_find_FrictWork");
    const bool   find_FrictWork_bh = captured.logical("_find_FrictWork_bh");
    const double backscat_subround = captured.real64("_backscat_subround");
    const bool   skeb_use_frict    = captured.logical("_skeb_use_frict");
    const double H_to_RZ           = captured.real64("_H_to_RZ");
    const double s_to_T            = captured.real64("_s_to_T");

    const auto u                       = captured.fab_device("_u");
    const auto v                       = captured.fab_device("_v");
    const auto h                       = captured.fab_device("_h");
    const auto uh                      = captured.fab_device("_uh");
    const auto vh                      = captured.fab_device("_vh");
    const auto str_xx_GME              = captured.fab_device("_str_xx_GME_before");
    const auto str_xy_GME              = captured.fab_device("_str_xy_GME_before");
    const auto visc_limit_h_flag       = captured.fab_device("_visc_limit_h_flag_before");
    const auto h_u                     = captured.fab_device("_h_u_before");
    const auto h_v                     = captured.fab_device("_h_v_before");
    const auto sh_xx                   = captured.fab_device("_sh_xx_before");
    const auto bhstr_xx                = captured.fab_device("_bhstr_xx_before");
    const auto sh_xy                   = captured.fab_device("_sh_xy_before");
    const auto bhstr_xy                = captured.fab_device("_bhstr_xy_before");
    const auto hq                      = captured.fab_device("_hq_before");
    auto       diffu                   = captured.fab_device("_diffu_before");
    auto       diffv                   = captured.fab_device("_diffv_before");
    auto       FrictWork               = captured.fab_device("_FrictWork_before");
    auto       FrictWork_bh            = captured.fab_device("_FrictWork_bh_before");
    auto       FrictWork_GME           = captured.fab_device("_FrictWork_GME_before");
    auto       str_xx                  = captured.fab_device("_str_xx_before");
    auto       str_xy                  = captured.fab_device("_str_xy_before");
    const auto diffu_after             = captured.fab_host("_diffu_after");
    const auto diffv_after             = captured.fab_host("_diffv_after");
    const auto FrictWork_after         = captured.fab_host("_FrictWork_after");
    const auto FrictWork_bh_after      = captured.fab_host("_FrictWork_bh_after");
    const auto FrictWork_GME_after     = captured.fab_host("_FrictWork_GME_after");
    const auto str_xx_after            = captured.fab_host("_str_xx_after");
    const auto str_xy_after            = captured.fab_host("_str_xy_after");

    MOM::hor_visc_stress_div(bxH, bxU, bxV, bxHstr, bxQstr, CS, MEKE_mom_src, MEKE_mom_src_bh, MEKE_GME_snk,
                             backscatter_Ro_pow, backscatter_Ro_c, skeb_diss, skeb_frict_coef, mask2dBu.const_array(),
                             IdxCu.const_array(), IdyCu.const_array(), IareaCu.const_array(), IdyCv.const_array(),
                             IdxCv.const_array(), IareaCv.const_array(), IdxT.const_array(), IdyT.const_array(),
                             IdyBu.const_array(), IdxBu.const_array(), IareaT.const_array(), dxCu.const_array(),
                             dyCv.const_array(), CoriolisBu.const_array(), reduction_xx, reduction_xy, dx2h, dy2h,
                             dx2q, dy2q, nullptr, kstart, h_neglect, apply_OBC, find_FrictWork, find_FrictWork_bh,
                             backscat_subround, skeb_use_frict, u.const_array(), v.const_array(), h.const_array(),
                             uh.const_array(), vh.const_array(), diffu.array(), diffv.array(),
                             str_xx_GME.const_array(), str_xy_GME.const_array(), FrictWork.array(),
                             FrictWork_bh.array(), FrictWork_GME.array(), visc_limit_h_flag.const_array(),
                             h_u.const_array(), h_v.const_array(), sh_xx.const_array(), str_xx.array(),
                             bhstr_xx.const_array(), sh_xy.const_array(), str_xy.array(), bhstr_xy.const_array(),
                             hq.const_array(), H_to_RZ, s_to_T);
    amrex::Gpu::synchronize();

    expect_arrays_equal(diffu_after,             to_host_fab(diffu),             "diffu");
    expect_arrays_equal(diffv_after,             to_host_fab(diffv),             "diffv");
    expect_arrays_equal(FrictWork_after,         to_host_fab(FrictWork),         "FrictWork");
    expect_arrays_equal(FrictWork_bh_after,      to_host_fab(FrictWork_bh),      "FrictWork_bh");
    expect_arrays_equal(FrictWork_GME_after,     to_host_fab(FrictWork_GME),     "FrictWork_GME");
    expect_arrays_equal(str_xx_after,            to_host_fab(str_xx),            "str_xx");
    expect_arrays_equal(str_xy_after,            to_host_fab(str_xy),            "str_xy");
    if (MEKE_mom_src_allocated) {
        expect_arrays_equal(captured.fab_host("_MEKE_mom_src_after"), to_host_fab(MEKE_mom_src_fab), "MEKE_mom_src");
    }
    if (MEKE_mom_src_bh_allocated) {
        expect_arrays_equal(captured.fab_host("_MEKE_mom_src_bh_after"), to_host_fab(MEKE_mom_src_bh_fab),
                            "MEKE_mom_src_bh");
    }
    if (MEKE_GME_snk_allocated) {
        expect_arrays_equal(captured.fab_host("_MEKE_GME_snk_after"), to_host_fab(MEKE_GME_snk_fab), "MEKE_GME_snk");
    }
    if (skeb_diss_allocated) {
        expect_arrays_equal(captured.fab_host("_STOCH_skeb_diss_after"), to_host_fab(skeb_diss_fab), "skeb_diss");
    }
}

// -------------------------------------------------------------------------
// Field mapping: MOM6 src/parameterizations/lateral/MOM_hor_visc.F90:5632-5751 (hor_visc_q_coeffs capture block).
TEST(HorViscQCoeffs, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "hor_visc_q_coeffs");

    const auto bxQstr   = captured.box("_bxQstr");
    const auto bxHgrad  = captured.box("_bxHgrad");
    const auto bxQshear = captured.box("_bxQshear");

    hv_opts_CS_C CS{};
    CS.Laplacian              = captured.logical("_CS_Laplacian");
    CS.biharmonic             = captured.logical("_CS_biharmonic");
    CS.debug                  = captured.logical("_CS_debug");
    CS.no_slip                = captured.logical("_CS_no_slip");
    CS.bound_Kh               = captured.logical("_CS_bound_Kh");
    CS.EY24_EBT_BS            = captured.logical("_CS_EY24_EBT_BS");
    CS.bound_Ah               = captured.logical("_CS_bound_Ah");
    CS.Re_Ah                  = captured.real64("_CS_Re_Ah");
    CS.KS_coef                = captured.real64("_CS_KS_coef");
    CS.backscatter_underbound = captured.logical("_CS_backscatter_underbound");
    CS.Smagorinsky_Kh         = captured.logical("_CS_Smagorinsky_Kh");
    CS.Smagorinsky_Ah         = captured.logical("_CS_Smagorinsky_Ah");
    CS.Leith_Kh               = captured.logical("_CS_Leith_Kh");
    CS.Modified_Leith         = captured.logical("_CS_Modified_Leith");
    CS.use_beta_in_Leith      = captured.logical("_CS_use_beta_in_Leith");
    CS.Leith_Ah               = captured.logical("_CS_Leith_Ah");
    CS.use_Leithy             = captured.logical("_CS_use_Leithy");
    CS.c_K                    = captured.real64("_CS_c_K");
    CS.smooth_Ah              = captured.logical("_CS_smooth_Ah");
    CS.use_QG_Leith_visc      = captured.logical("_CS_use_QG_Leith_visc");
    CS.bound_Coriolis         = captured.logical("_CS_bound_Coriolis");
    CS.Kh_bg_min              = captured.real64("_CS_Kh_bg_min");
    CS.FrictWork_bug          = captured.logical("_CS_FrictWork_bug");
    CS.OBC_strain_bug         = captured.logical("_CS_OBC_strain_bug");
    CS.use_land_mask          = captured.logical("_CS_use_land_mask");
    CS.anisotropic            = captured.logical("_CS_anisotropic");
    CS.add_LES_viscosity      = captured.logical("_CS_add_LES_viscosity");
    CS.Kh_aniso               = captured.real64("_CS_Kh_aniso");
    CS.res_scale_MEKE         = captured.logical("_CS_res_scale_MEKE");
    CS.use_GME                = captured.logical("_CS_use_GME");
    CS.answer_date            = captured.integer("_CS_answer_date");
    CS.GME_h0                 = captured.real64("_CS_GME_h0");
    CS.GME_efficiency         = captured.real64("_CS_GME_efficiency");
    CS.GME_limiter            = captured.real64("_CS_GME_limiter");
    CS.min_grid_Kh            = captured.real64("_CS_min_grid_Kh");
    CS.min_grid_Ah            = captured.real64("_CS_min_grid_Ah");
    CS.use_cont_thick         = captured.logical("_CS_use_cont_thick");
    CS.use_ZB2020             = captured.logical("_CS_use_ZB2020");
    CS.use_circulation        = captured.logical("_CS_use_circulation");
    CS.do_grid_Re_Ah          = captured.logical("_CS_do_grid_Re_Ah");
    CS.do_grid_Re_Kh          = captured.logical("_CS_do_grid_Re_Kh");
    CS.do_Ah_h                = captured.logical("_CS_do_Ah_h");
    CS.do_Ah_q                = captured.logical("_CS_do_Ah_q");
    CS.do_Kh_h                = captured.logical("_CS_do_Kh_h");
    CS.do_Kh_q                = captured.logical("_CS_do_Kh_q");
    CS.do_GME_coeff_h         = captured.logical("_CS_do_GME_coeff_h");
    CS.do_GME_coeff_q         = captured.logical("_CS_do_GME_coeff_q");
    CS.do_vort_xy_q           = captured.logical("_CS_do_vort_xy_q");
    CS.do_div_xx_h            = captured.logical("_CS_do_div_xx_h");
    CS.do_sh_xy_q             = captured.logical("_CS_do_sh_xy_q");
    CS.do_sh_xx_h             = captured.logical("_CS_do_sh_xx_h");
    CS.do_FrictWork           = captured.logical("_CS_do_FrictWork");
    CS.do_FrictWorkIntz       = captured.logical("_CS_do_FrictWorkIntz");
    CS.do_FrictWork_bh        = captured.logical("_CS_do_FrictWork_bh");
    CS.do_FrictWorkIntz_bh    = captured.logical("_CS_do_FrictWorkIntz_bh");
    CS.do_normstress          = captured.logical("_CS_do_normstress");
    CS.do_shearstress         = captured.logical("_CS_do_shearstress");
    CS.do_BS_coeff_h          = captured.logical("_CS_do_BS_coeff_h");
    CS.do_BS_coeff_q          = captured.logical("_CS_do_BS_coeff_q");

    // The MEKE-derived lateral viscosity coefficient [L2 T-1 ~> m2 s-1]; may be unallocated
    amrex::FArrayBox MEKE_Ku_fab;
    amrex::Array4<const amrex::Real> MEKE_Ku{};
    if (captured.logical("_MEKE_Ku_allocated")) {
        MEKE_Ku_fab = captured.fab_device("_MEKE_Ku_before");
        MEKE_Ku = MEKE_Ku_fab.const_array();
    }

    // The MEKE-derived lateral biharmonic viscosity coefficient [L4 T-1 ~> m4 s-1]; may be unallocated
    amrex::FArrayBox MEKE_Au_fab;
    amrex::Array4<const amrex::Real> MEKE_Au{};
    if (captured.logical("_MEKE_Au_allocated")) {
        MEKE_Au_fab = captured.fab_device("_MEKE_Au_before");
        MEKE_Au = MEKE_Au_fab.const_array();
    }

    // Non-dimensional function of the ratio the first baroclinic deformation radius to the grid spacing at q
    // points [nondim]; may be unallocated
    amrex::FArrayBox Res_fn_q_fab;
    amrex::Array4<const amrex::Real> Res_fn_q{};
    if (captured.logical("_VarMix_Res_fn_q_allocated")) {
        Res_fn_q_fab = captured.fab_device("_VarMix_Res_fn_q_before");
        Res_fn_q = Res_fn_q_fab.const_array();
    }

    // Vertical structure function used in backscatter [nondim]; may be unallocated
    amrex::FArrayBox BS_struct_fab;
    amrex::Array4<const amrex::Real> BS_struct{};
    if (captured.logical("_VarMix_BS_struct_allocated")) {
        BS_struct_fab = captured.fab_device("_VarMix_BS_struct_before");
        BS_struct = BS_struct_fab.const_array();
    }

    // The background Laplacian viscosity at q points [L2 T-1 ~> m2 s-1]; may be unallocated
    amrex::FArrayBox Kh_bg_xy_fab;
    amrex::Array4<const amrex::Real> Kh_bg_xy{};
    if (captured.logical("_metrics_Kh_bg_xy_allocated")) {
        Kh_bg_xy_fab = captured.fab_device("_metrics_Kh_bg_xy_before");
        Kh_bg_xy = Kh_bg_xy_fab.const_array();
    }

    // The background biharmonic viscosity at q points [L4 T-1 ~> m4 s-1]; may be unallocated
    amrex::FArrayBox Ah_bg_xy_fab;
    amrex::Array4<const amrex::Real> Ah_bg_xy{};
    if (captured.logical("_metrics_Ah_bg_xy_allocated")) {
        Ah_bg_xy_fab = captured.fab_device("_metrics_Ah_bg_xy_before");
        Ah_bg_xy = Ah_bg_xy_fab.const_array();
    }

    // The amount by which stresses through q points are reduced due to partial barriers [nondim]; may be unallocated
    amrex::FArrayBox reduction_xy_fab;
    amrex::Array4<const amrex::Real> reduction_xy{};
    if (captured.logical("_metrics_reduction_xy_allocated")) {
        reduction_xy_fab = captured.fab_device("_metrics_reduction_xy_before");
        reduction_xy = reduction_xy_fab.const_array();
    }

    // The maximum permitted Laplacian viscosity [L2 T-1 ~> m2 s-1]; may be unallocated
    amrex::FArrayBox Kh_Max_xy_fab;
    amrex::Array4<const amrex::Real> Kh_Max_xy{};
    if (captured.logical("_metrics_Kh_Max_xy_allocated")) {
        Kh_Max_xy_fab = captured.fab_device("_metrics_Kh_Max_xy_before");
        Kh_Max_xy = Kh_Max_xy_fab.const_array();
    }

    // The maximum permitted biharmonic viscosity [L4 T-1 ~> m4 s-1]; may be unallocated
    amrex::FArrayBox Ah_Max_xy_fab;
    amrex::Array4<const amrex::Real> Ah_Max_xy{};
    if (captured.logical("_metrics_Ah_Max_xy_allocated")) {
        Ah_Max_xy_fab = captured.fab_device("_metrics_Ah_Max_xy_before");
        Ah_Max_xy = Ah_Max_xy_fab.const_array();
    }

    // The maximum permitted biharmonic viscosity for the kill switch [L4 T-1 ~> m4 s-1]; may be unallocated
    amrex::FArrayBox Ah_Max_xy_KS_fab;
    amrex::Array4<const amrex::Real> Ah_Max_xy_KS{};
    if (captured.logical("_metrics_Ah_Max_xy_KS_allocated")) {
        Ah_Max_xy_KS_fab = captured.fab_device("_metrics_Ah_Max_xy_KS_before");
        Ah_Max_xy_KS = Ah_Max_xy_KS_fab.const_array();
    }

    // Factor n1*n2 in the anisotropic direction tensor at q-points [nondim]; may be unallocated
    amrex::FArrayBox n1n2_q_fab;
    amrex::Array4<const amrex::Real> n1n2_q{};
    if (captured.logical("_metrics_n1n2_q_allocated")) {
        n1n2_q_fab = captured.fab_device("_metrics_n1n2_q_before");
        n1n2_q = n1n2_q_fab.const_array();
    }

    // Factor n1**2-n2**2 in the anisotropic direction tensor at q-points [nondim]; may be unallocated
    amrex::FArrayBox n1n1_m_n2n2_q_fab;
    amrex::Array4<const amrex::Real> n1n1_m_n2n2_q{};
    if (captured.logical("_metrics_n1n1_m_n2n2_q_allocated")) {
        n1n1_m_n2n2_q_fab = captured.fab_device("_metrics_n1n1_m_n2n2_q_before");
        n1n1_m_n2n2_q = n1n1_m_n2n2_q_fab.const_array();
    }

    // Laplacian metric-dependent constants [L2 ~> m2]; may be unallocated
    amrex::FArrayBox Laplac2_const_xy_fab;
    amrex::Array4<const amrex::Real> Laplac2_const_xy{};
    if (captured.logical("_metrics_Laplac2_const_xy_allocated")) {
        Laplac2_const_xy_fab = captured.fab_device("_metrics_Laplac2_const_xy_before");
        Laplac2_const_xy = Laplac2_const_xy_fab.const_array();
    }

    // Biharmonic metric-dependent constants [L6 ~> m6]; may be unallocated
    amrex::FArrayBox Biharm6_const_xy_fab;
    amrex::Array4<const amrex::Real> Biharm6_const_xy{};
    if (captured.logical("_metrics_Biharm6_const_xy_allocated")) {
        Biharm6_const_xy_fab = captured.fab_device("_metrics_Biharm6_const_xy_before");
        Biharm6_const_xy = Biharm6_const_xy_fab.const_array();
    }

    // Laplacian metric-dependent constants [L3 ~> m3]; may be unallocated
    amrex::FArrayBox Laplac3_const_xy_fab;
    amrex::Array4<const amrex::Real> Laplac3_const_xy{};
    if (captured.logical("_metrics_Laplac3_const_xy_allocated")) {
        Laplac3_const_xy_fab = captured.fab_device("_metrics_Laplac3_const_xy_before");
        Laplac3_const_xy = Laplac3_const_xy_fab.const_array();
    }

    // Biharmonic metric-dependent constants [L4 ~> m4]; may be unallocated
    amrex::FArrayBox Biharm_const_xy_fab;
    amrex::Array4<const amrex::Real> Biharm_const_xy{};
    if (captured.logical("_metrics_Biharm_const_xy_allocated")) {
        Biharm_const_xy_fab = captured.fab_device("_metrics_Biharm_const_xy_before");
        Biharm_const_xy = Biharm_const_xy_fab.const_array();
    }

    // Biharmonic metric-dependent constants [T L4 ~> s m4]; may be unallocated
    amrex::FArrayBox Biharm_const2_xy_fab;
    amrex::Array4<const amrex::Real> Biharm_const2_xy{};
    if (captured.logical("_metrics_Biharm_const2_xy_allocated")) {
        Biharm_const2_xy_fab = captured.fab_device("_metrics_Biharm_const2_xy_before");
        Biharm_const2_xy = Biharm_const2_xy_fab.const_array();
    }

    // Biharmonic metric-dependent constants [L3 ~> m3]; may be unallocated
    amrex::FArrayBox Re_Ah_const_xy_fab;
    amrex::Array4<const amrex::Real> Re_Ah_const_xy{};
    if (captured.logical("_metrics_Re_Ah_const_xy_allocated")) {
        Re_Ah_const_xy_fab = captured.fab_device("_metrics_Re_Ah_const_xy_before");
        Re_Ah_const_xy = Re_Ah_const_xy_fab.const_array();
    }

    const int    kstart         = captured.integer("_kstart");
    const double inv_PI3        = captured.real64("_inv_PI3");
    const double inv_PI6        = captured.real64("_inv_PI6");
    const bool   find_FrictWork = captured.logical("_find_FrictWork");
    const bool   use_kh_struct  = captured.logical("_use_kh_struct");
    const bool   use_MEKE_Ku    = captured.logical("_use_MEKE_Ku");
    const bool   use_MEKE_Au    = captured.logical("_use_MEKE_Au");
    const bool   rescale_Kh     = captured.logical("_rescale_Kh");

    const auto mask2dBu                 = captured.fab_device("_mask2dBu");
    const auto u                        = captured.fab_device("_u");
    const auto v                        = captured.fab_device("_v");
    const auto sh_xx_bt                 = captured.fab_device("_sh_xx_bt");
    const auto sh_xy_bt                 = captured.fab_device("_sh_xy_bt");
    const auto GME_effic_h              = captured.fab_device("_GME_effic_h");
    const auto GME_effic_q              = captured.fab_device("_GME_effic_q");
    const auto KH_u_GME                 = captured.fab_device("_KH_u_GME");
    const auto KH_v_GME                 = captured.fab_device("_KH_v_GME");
    const auto Ah_h                     = captured.fab_device("_Ah_h_before");
    const auto Kh_h                     = captured.fab_device("_Kh_h_before");
    const auto sh_xx                    = captured.fab_device("_sh_xx_before");
    const auto dDel2vdx                 = captured.fab_device("_dDel2vdx_before");
    const auto dDel2udy                 = captured.fab_device("_dDel2udy_before");
    const auto sh_xy                    = captured.fab_device("_sh_xy_before");
    const auto sh_xy_smooth             = captured.fab_device("_sh_xy_smooth_before");
    const auto vort_xy                  = captured.fab_device("_vort_xy_before");
    const auto grad_vort_mag_q          = captured.fab_device("_grad_vort_mag_q_before");
    const auto grad_vort_mag_q_2d       = captured.fab_device("_grad_vort_mag_q_2d_before");
    const auto Del2vort_q               = captured.fab_device("_Del2vort_q_before");
    const auto grad_div_mag_q           = captured.fab_device("_grad_div_mag_q_before");
    const auto hq                       = captured.fab_device("_hq_before");
    const auto Shear_mag                = captured.fab_device("_Shear_mag_before");
    const auto hrat_min                 = captured.fab_device("_hrat_min_before");
    auto       GME_coeff_h              = captured.fab_device("_GME_coeff_h_before");
    auto       GME_coeff_q              = captured.fab_device("_GME_coeff_q_before");
    auto       str_xx_GME               = captured.fab_device("_str_xx_GME_before");
    auto       str_xy_GME               = captured.fab_device("_str_xy_GME_before");
    auto       Ah_q                     = captured.fab_device("_Ah_q_before");
    auto       Kh_q                     = captured.fab_device("_Kh_q_before");
    auto       vort_xy_q                = captured.fab_device("_vort_xy_q_before");
    auto       sh_xy_q                  = captured.fab_device("_sh_xy_q_before");
    auto       visc_limit_q             = captured.fab_device("_visc_limit_q_before");
    auto       visc_limit_q_flag        = captured.fab_device("_visc_limit_q_flag_before");
    auto       visc_limit_q_frac        = captured.fab_device("_visc_limit_q_frac_before");
    auto       BS_coeff_q               = captured.fab_device("_BS_coeff_q_before");
    auto       str_xy                   = captured.fab_device("_str_xy_before");
    auto       bhstr_xy                 = captured.fab_device("_bhstr_xy_before");
    auto       Ah                       = captured.fab_device("_Ah_before");
    auto       Kh                       = captured.fab_device("_Kh_before");
    auto       vert_vort_mag            = captured.fab_device("_vert_vort_mag_before");
    auto       visc_bound_rem           = captured.fab_device("_visc_bound_rem_before");
    const auto GME_coeff_h_after        = captured.fab_host("_GME_coeff_h_after");
    const auto GME_coeff_q_after        = captured.fab_host("_GME_coeff_q_after");
    const auto str_xx_GME_after         = captured.fab_host("_str_xx_GME_after");
    const auto str_xy_GME_after         = captured.fab_host("_str_xy_GME_after");
    const auto Ah_q_after               = captured.fab_host("_Ah_q_after");
    const auto Kh_q_after               = captured.fab_host("_Kh_q_after");
    const auto vort_xy_q_after          = captured.fab_host("_vort_xy_q_after");
    const auto sh_xy_q_after            = captured.fab_host("_sh_xy_q_after");
    const auto visc_limit_q_after       = captured.fab_host("_visc_limit_q_after");
    const auto visc_limit_q_flag_after  = captured.fab_host("_visc_limit_q_flag_after");
    const auto visc_limit_q_frac_after  = captured.fab_host("_visc_limit_q_frac_after");
    const auto BS_coeff_q_after         = captured.fab_host("_BS_coeff_q_after");
    const auto str_xy_after             = captured.fab_host("_str_xy_after");
    const auto bhstr_xy_after           = captured.fab_host("_bhstr_xy_after");
    const auto Ah_after                 = captured.fab_host("_Ah_after");
    const auto Kh_after                 = captured.fab_host("_Kh_after");
    const auto vert_vort_mag_after      = captured.fab_host("_vert_vort_mag_after");
    const auto visc_bound_rem_after     = captured.fab_host("_visc_bound_rem_after");

    MOM::hor_visc_q_coeffs(bxQstr, bxHgrad, bxQshear, CS, MEKE_Ku, MEKE_Au, Res_fn_q, BS_struct, mask2dBu.const_array(),
                           Kh_bg_xy, Ah_bg_xy, reduction_xy, Kh_Max_xy, Ah_Max_xy, Ah_Max_xy_KS, n1n2_q, n1n1_m_n2n2_q,
                           Laplac2_const_xy, Biharm6_const_xy, Laplac3_const_xy, Biharm_const_xy, Biharm_const2_xy,
                           Re_Ah_const_xy, kstart, inv_PI3, inv_PI6, find_FrictWork, use_kh_struct, use_MEKE_Ku,
                           use_MEKE_Au, rescale_Kh, u.const_array(), v.const_array(), sh_xx_bt.const_array(),
                           sh_xy_bt.const_array(), GME_effic_h.const_array(), GME_effic_q.const_array(),
                           KH_u_GME.const_array(), KH_v_GME.const_array(), GME_coeff_h.array(), GME_coeff_q.array(),
                           str_xx_GME.array(), str_xy_GME.array(), Ah_h.const_array(), Ah_q.array(),
                           Kh_h.const_array(), Kh_q.array(), vort_xy_q.array(), sh_xy_q.array(),
                           visc_limit_q.array(), visc_limit_q_flag.array(), visc_limit_q_frac.array(),
                           BS_coeff_q.array(), sh_xx.const_array(), dDel2vdx.const_array(), dDel2udy.const_array(),
                           sh_xy.const_array(), sh_xy_smooth.const_array(), str_xy.array(), bhstr_xy.array(),
                           vort_xy.const_array(), grad_vort_mag_q.const_array(), grad_vort_mag_q_2d.const_array(),
                           Del2vort_q.const_array(), grad_div_mag_q.const_array(), hq.const_array(), Ah.array(),
                           Kh.array(), Shear_mag.const_array(), vert_vort_mag.array(), hrat_min.const_array(),
                           visc_bound_rem.array());
    amrex::Gpu::synchronize();

    expect_arrays_equal(GME_coeff_h_after,        to_host_fab(GME_coeff_h),        "GME_coeff_h");
    expect_arrays_equal(GME_coeff_q_after,        to_host_fab(GME_coeff_q),        "GME_coeff_q");
    expect_arrays_equal(str_xx_GME_after,         to_host_fab(str_xx_GME),         "str_xx_GME");
    expect_arrays_equal(str_xy_GME_after,         to_host_fab(str_xy_GME),         "str_xy_GME");
    expect_arrays_equal(Ah_q_after,               to_host_fab(Ah_q),               "Ah_q");
    expect_arrays_equal(Kh_q_after,               to_host_fab(Kh_q),               "Kh_q");
    expect_arrays_equal(vort_xy_q_after,          to_host_fab(vort_xy_q),          "vort_xy_q");
    expect_arrays_equal(sh_xy_q_after,            to_host_fab(sh_xy_q),            "sh_xy_q");
    expect_arrays_equal(visc_limit_q_after,       to_host_fab(visc_limit_q),       "visc_limit_q");
    expect_arrays_equal(visc_limit_q_flag_after,  to_host_fab(visc_limit_q_flag),  "visc_limit_q_flag");
    expect_arrays_equal(visc_limit_q_frac_after,  to_host_fab(visc_limit_q_frac),  "visc_limit_q_frac");
    expect_arrays_equal(BS_coeff_q_after,         to_host_fab(BS_coeff_q),         "BS_coeff_q");
    expect_arrays_equal(str_xy_after,             to_host_fab(str_xy),             "str_xy");
    expect_arrays_equal(bhstr_xy_after,           to_host_fab(bhstr_xy),           "bhstr_xy");
    expect_arrays_equal(Ah_after,                 to_host_fab(Ah),                 "Ah");
    expect_arrays_equal(Kh_after,                 to_host_fab(Kh),                 "Kh");
    expect_arrays_equal(vert_vort_mag_after,      to_host_fab(vert_vort_mag),      "vert_vort_mag");
    expect_arrays_equal(visc_bound_rem_after,     to_host_fab(visc_bound_rem),     "visc_bound_rem");
}

// -------------------------------------------------------------------------
// Field mapping: MOM6 src/parameterizations/lateral/MOM_hor_visc.F90:2665-2813 (horizontal_viscosity_TR capture block).
TEST(HorizontalViscosityTr, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "horizontal_viscosity_tr");

    const auto bxH0 = captured.box("_bxH0");
    const auto bxQ0 = captured.box("_bxQ0");

    hv_opts_CS_C CS{};
    CS.Laplacian              = captured.logical("_CS_Laplacian");
    CS.biharmonic             = captured.logical("_CS_biharmonic");
    CS.debug                  = captured.logical("_CS_debug");
    CS.no_slip                = captured.logical("_CS_no_slip");
    CS.bound_Kh               = captured.logical("_CS_bound_Kh");
    CS.EY24_EBT_BS            = captured.logical("_CS_EY24_EBT_BS");
    CS.bound_Ah               = captured.logical("_CS_bound_Ah");
    CS.Re_Ah                  = captured.real64("_CS_Re_Ah");
    CS.KS_coef                = captured.real64("_CS_KS_coef");
    CS.backscatter_underbound = captured.logical("_CS_backscatter_underbound");
    CS.Smagorinsky_Kh         = captured.logical("_CS_Smagorinsky_Kh");
    CS.Smagorinsky_Ah         = captured.logical("_CS_Smagorinsky_Ah");
    CS.Leith_Kh               = captured.logical("_CS_Leith_Kh");
    CS.Modified_Leith         = captured.logical("_CS_Modified_Leith");
    CS.use_beta_in_Leith      = captured.logical("_CS_use_beta_in_Leith");
    CS.Leith_Ah               = captured.logical("_CS_Leith_Ah");
    CS.use_Leithy             = captured.logical("_CS_use_Leithy");
    CS.c_K                    = captured.real64("_CS_c_K");
    CS.smooth_Ah              = captured.logical("_CS_smooth_Ah");
    CS.use_QG_Leith_visc      = captured.logical("_CS_use_QG_Leith_visc");
    CS.bound_Coriolis         = captured.logical("_CS_bound_Coriolis");
    CS.Kh_bg_min              = captured.real64("_CS_Kh_bg_min");
    CS.FrictWork_bug          = captured.logical("_CS_FrictWork_bug");
    CS.OBC_strain_bug         = captured.logical("_CS_OBC_strain_bug");
    CS.use_land_mask          = captured.logical("_CS_use_land_mask");
    CS.anisotropic            = captured.logical("_CS_anisotropic");
    CS.add_LES_viscosity      = captured.logical("_CS_add_LES_viscosity");
    CS.Kh_aniso               = captured.real64("_CS_Kh_aniso");
    CS.res_scale_MEKE         = captured.logical("_CS_res_scale_MEKE");
    CS.use_GME                = captured.logical("_CS_use_GME");
    CS.answer_date            = captured.integer("_CS_answer_date");
    CS.GME_h0                 = captured.real64("_CS_GME_h0");
    CS.GME_efficiency         = captured.real64("_CS_GME_efficiency");
    CS.GME_limiter            = captured.real64("_CS_GME_limiter");
    CS.min_grid_Kh            = captured.real64("_CS_min_grid_Kh");
    CS.min_grid_Ah            = captured.real64("_CS_min_grid_Ah");
    CS.use_cont_thick         = captured.logical("_CS_use_cont_thick");
    CS.use_ZB2020             = captured.logical("_CS_use_ZB2020");
    CS.use_circulation        = captured.logical("_CS_use_circulation");
    CS.do_grid_Re_Ah          = captured.logical("_CS_do_grid_Re_Ah");
    CS.do_grid_Re_Kh          = captured.logical("_CS_do_grid_Re_Kh");
    CS.do_Ah_h                = captured.logical("_CS_do_Ah_h");
    CS.do_Ah_q                = captured.logical("_CS_do_Ah_q");
    CS.do_Kh_h                = captured.logical("_CS_do_Kh_h");
    CS.do_Kh_q                = captured.logical("_CS_do_Kh_q");
    CS.do_GME_coeff_h         = captured.logical("_CS_do_GME_coeff_h");
    CS.do_GME_coeff_q         = captured.logical("_CS_do_GME_coeff_q");
    CS.do_vort_xy_q           = captured.logical("_CS_do_vort_xy_q");
    CS.do_div_xx_h            = captured.logical("_CS_do_div_xx_h");
    CS.do_sh_xy_q             = captured.logical("_CS_do_sh_xy_q");
    CS.do_sh_xx_h             = captured.logical("_CS_do_sh_xx_h");
    CS.do_FrictWork           = captured.logical("_CS_do_FrictWork");
    CS.do_FrictWorkIntz       = captured.logical("_CS_do_FrictWorkIntz");
    CS.do_FrictWork_bh        = captured.logical("_CS_do_FrictWork_bh");
    CS.do_FrictWorkIntz_bh    = captured.logical("_CS_do_FrictWorkIntz_bh");
    CS.do_normstress          = captured.logical("_CS_do_normstress");
    CS.do_shearstress         = captured.logical("_CS_do_shearstress");
    CS.do_BS_coeff_h          = captured.logical("_CS_do_BS_coeff_h");
    CS.do_BS_coeff_q          = captured.logical("_CS_do_BS_coeff_q");

    // Configurations horizontal_viscosity_tr does not support yet
    if (captured.logical("_OBC_associated")) {
        FAIL() << "capture has open boundary conditions, which horizontal_viscosity_tr does not support";
    }
    if (CS.use_GME) {
        FAIL() << "capture uses GME (USE_GME), which horizontal_viscosity_tr does not support";
    }
    if (CS.use_QG_Leith_visc) {
        FAIL() << "capture uses USE_QG_LEITH_VISC, which horizontal_viscosity_tr does not support";
    }
    if (CS.use_ZB2020) {
        FAIL() << "capture uses ZB2020, which horizontal_viscosity_tr does not support";
    }
    if (CS.biharmonic && CS.use_Leithy && CS.smooth_Ah) {
        FAIL() << "capture uses SMOOTH_AH with Leith+E, which horizontal_viscosity_tr does not support";
    }

    // Core ocean grid fields
    const auto mask2dT_fab    = captured.fab_device("_mask2dT");
    const auto IareaT_fab     = captured.fab_device("_IareaT");
    const auto IdxT_fab       = captured.fab_device("_IdxT");
    const auto IdyT_fab       = captured.fab_device("_IdyT");
    const auto mask2dCu_fab   = captured.fab_device("_mask2dCu");
    const auto IdxCu_fab      = captured.fab_device("_IdxCu");
    const auto dxCu_fab       = captured.fab_device("_dxCu");
    const auto mask2dCv_fab   = captured.fab_device("_mask2dCv");
    const auto IdyCv_fab      = captured.fab_device("_IdyCv");
    const auto dyCv_fab       = captured.fab_device("_dyCv");
    const auto CoriolisBu_fab = captured.fab_device("_CoriolisBu");
    const auto IdyCu_fab      = captured.fab_device("_IdyCu");
    const auto IdxCv_fab      = captured.fab_device("_IdxCv");
    const auto mask2dBu_fab   = captured.fab_device("_mask2dBu");
    const auto dxBu_fab       = captured.fab_device("_dxBu");
    const auto dyBu_fab       = captured.fab_device("_dyBu");
    const auto IareaBu_fab    = captured.fab_device("_IareaBu");
    const auto df_dx_fab      = captured.fab_device("_dF_dx");
    const auto df_dy_fab      = captured.fab_device("_dF_dy");
    const auto IareaCu_fab    = captured.fab_device("_IareaCu");
    const auto IareaCv_fab    = captured.fab_device("_IareaCv");
    const auto IdxBu_fab      = captured.fab_device("_IdxBu");
    const auto IdyBu_fab      = captured.fab_device("_IdyBu");
    MOM::hv_grid_dev G;
    G.mask2dT = mask2dT_fab.const_array();
    G.IareaT = IareaT_fab.const_array();
    G.IdxT = IdxT_fab.const_array();
    G.IdyT = IdyT_fab.const_array();
    G.mask2dCu = mask2dCu_fab.const_array();
    G.IdxCu = IdxCu_fab.const_array();
    G.dxCu = dxCu_fab.const_array();
    G.mask2dCv = mask2dCv_fab.const_array();
    G.IdyCv = IdyCv_fab.const_array();
    G.dyCv = dyCv_fab.const_array();
    G.CoriolisBu = CoriolisBu_fab.const_array();
    G.IdyCu = IdyCu_fab.const_array();
    G.IdxCv = IdxCv_fab.const_array();
    G.mask2dBu = mask2dBu_fab.const_array();
    G.dxBu = dxBu_fab.const_array();
    G.dyBu = dyBu_fab.const_array();
    G.IareaBu = IareaBu_fab.const_array();
    G.df_dx = df_dx_fab.const_array();
    G.df_dy = df_dy_fab.const_array();
    G.IareaCu = IareaCu_fab.const_array();
    G.IareaCv = IareaCv_fab.const_array();
    G.IdxBu = IdxBu_fab.const_array();
    G.IdyBu = IdyBu_fab.const_array();

    // Horizontal viscosity metric and coefficient arrays; each may be unallocated
    MOM::hv_metrics_dev M;
    std::vector<amrex::FArrayBox> metrics_fabs;
    metrics_fabs.reserve(44);
    const auto load_metric = [&](const char* name, amrex::Array4<const amrex::Real>& view) {
        const std::string field = std::string("_metrics_") + name;
        if (captured.logical(field + "_allocated")) {
            metrics_fabs.push_back(captured.fab_device(field + "_before"));
            view = metrics_fabs.back().const_array();
        }
    };
    load_metric("Kh_bg_xx", M.Kh_bg_xx);
    load_metric("Ah_bg_xx", M.Ah_bg_xx);
    load_metric("reduction_xx", M.reduction_xx);
    load_metric("Kh_Max_xx", M.Kh_Max_xx);
    load_metric("Ah_Max_xx", M.Ah_Max_xx);
    load_metric("Ah_Max_xx_KS", M.Ah_Max_xx_KS);
    load_metric("n1n2_h", M.n1n2_h);
    load_metric("n1n1_m_n2n2_h", M.n1n1_m_n2n2_h);
    load_metric("grid_sp_h2", M.grid_sp_h2);
    load_metric("grid_sp_h3", M.grid_sp_h3);
    load_metric("Kh_bg_xy", M.Kh_bg_xy);
    load_metric("Ah_bg_xy", M.Ah_bg_xy);
    load_metric("reduction_xy", M.reduction_xy);
    load_metric("Kh_Max_xy", M.Kh_Max_xy);
    load_metric("Ah_Max_xy", M.Ah_Max_xy);
    load_metric("Ah_Max_xy_KS", M.Ah_Max_xy_KS);
    load_metric("n1n2_q", M.n1n2_q);
    load_metric("n1n1_m_n2n2_q", M.n1n1_m_n2n2_q);
    load_metric("dx2h", M.dx2h);
    load_metric("dy2h", M.dy2h);
    load_metric("dx_dyT", M.dx_dyT);
    load_metric("dy_dxT", M.dy_dxT);
    load_metric("m_const_leithy", M.m_const_leithy);
    load_metric("m_leithy_max", M.m_leithy_max);
    load_metric("dx2q", M.dx2q);
    load_metric("dy2q", M.dy2q);
    load_metric("dx_dyBu", M.dx_dyBu);
    load_metric("dy_dxBu", M.dy_dxBu);
    load_metric("Idx2dyCu", M.Idx2dyCu);
    load_metric("Idxdy2u", M.Idxdy2u);
    load_metric("Idx2dyCv", M.Idx2dyCv);
    load_metric("Idxdy2v", M.Idxdy2v);
    load_metric("Laplac2_const_xx", M.Laplac2_const_xx);
    load_metric("Biharm6_const_xx", M.Biharm6_const_xx);
    load_metric("Laplac3_const_xx", M.Laplac3_const_xx);
    load_metric("Biharm_const_xx", M.Biharm_const_xx);
    load_metric("Biharm_const2_xx", M.Biharm_const2_xx);
    load_metric("Re_Ah_const_xx", M.Re_Ah_const_xx);
    load_metric("Laplac2_const_xy", M.Laplac2_const_xy);
    load_metric("Biharm6_const_xy", M.Biharm6_const_xy);
    load_metric("Laplac3_const_xy", M.Laplac3_const_xy);
    load_metric("Biharm_const_xy", M.Biharm_const_xy);
    load_metric("Biharm_const2_xy", M.Biharm_const2_xy);
    load_metric("Re_Ah_const_xy", M.Re_Ah_const_xy);

    // MEKE, VarMix and STOCH fields; each array may be unallocated
    MOM::hv_mixing_dev X;
    amrex::FArrayBox MEKE_Ku_fab;
    if (captured.logical("_MEKE_Ku_allocated")) {
        MEKE_Ku_fab = captured.fab_device("_MEKE_Ku_before");
        X.MEKE_Ku = MEKE_Ku_fab.const_array();
    }
    amrex::FArrayBox MEKE_Au_fab;
    if (captured.logical("_MEKE_Au_allocated")) {
        MEKE_Au_fab = captured.fab_device("_MEKE_Au_before");
        X.MEKE_Au = MEKE_Au_fab.const_array();
    }
    amrex::FArrayBox Res_fn_h_fab;
    if (captured.logical("_VarMix_Res_fn_h_allocated")) {
        Res_fn_h_fab = captured.fab_device("_VarMix_Res_fn_h_before");
        X.Res_fn_h = Res_fn_h_fab.const_array();
    }
    amrex::FArrayBox Res_fn_q_fab;
    if (captured.logical("_VarMix_Res_fn_q_allocated")) {
        Res_fn_q_fab = captured.fab_device("_VarMix_Res_fn_q_before");
        X.Res_fn_q = Res_fn_q_fab.const_array();
    }
    amrex::FArrayBox BS_struct_fab;
    if (captured.logical("_VarMix_BS_struct_allocated")) {
        BS_struct_fab = captured.fab_device("_VarMix_BS_struct_before");
        X.BS_struct = BS_struct_fab.const_array();
    }
    const bool MEKE_mom_src_allocated = captured.logical("_MEKE_mom_src_allocated");
    amrex::FArrayBox MEKE_mom_src_fab;
    if (MEKE_mom_src_allocated) {
        MEKE_mom_src_fab = captured.fab_device("_MEKE_mom_src_before");
        X.MEKE_mom_src = MEKE_mom_src_fab.array();
    }
    const bool MEKE_mom_src_bh_allocated = captured.logical("_MEKE_mom_src_bh_allocated");
    amrex::FArrayBox MEKE_mom_src_bh_fab;
    if (MEKE_mom_src_bh_allocated) {
        MEKE_mom_src_bh_fab = captured.fab_device("_MEKE_mom_src_bh_before");
        X.MEKE_mom_src_bh = MEKE_mom_src_bh_fab.array();
    }
    const bool MEKE_GME_snk_allocated = captured.logical("_MEKE_GME_snk_allocated");
    amrex::FArrayBox MEKE_GME_snk_fab;
    if (MEKE_GME_snk_allocated) {
        MEKE_GME_snk_fab = captured.fab_device("_MEKE_GME_snk_before");
        X.MEKE_GME_snk = MEKE_GME_snk_fab.array();
    }
    const bool skeb_diss_allocated = captured.logical("_STOCH_skeb_diss_allocated");
    amrex::FArrayBox skeb_diss_fab;
    if (skeb_diss_allocated) {
        skeb_diss_fab = captured.fab_device("_STOCH_skeb_diss_before");
        X.skeb_diss = skeb_diss_fab.array();
    }
    X.backscatter_Ro_pow  = captured.real64("_MEKE_backscatter_Ro_pow");
    X.backscatter_Ro_c    = captured.real64("_MEKE_backscatter_Ro_c");
    X.use_variable_mixing = captured.logical("_VarMix_use_variable_mixing");
    X.Resoln_scaled_Kh    = captured.logical("_VarMix_Resoln_scaled_Kh");
    X.skeb_use_frict      = captured.logical("_STOCH_skeb_use_frict");
    X.skeb_frict_coef     = captured.real64("_STOCH_skeb_frict_coef");

    const int    nkblock       = captured.integer("_nkblock");
    const double H_subroundoff = captured.real64("_H_subroundoff");
    const double H_to_RZ       = captured.real64("_H_to_RZ");
    const double s_to_T        = captured.real64("_s_to_T");

    // Layer thickness at u-points [H ~> m or kg m-2]; may be unallocated
    amrex::FArrayBox hu_cont_fab;
    amrex::Array4<const amrex::Real> hu_cont{};
    if (captured.logical("_hu_cont_allocated")) {
        hu_cont_fab = captured.fab_device("_hu_cont_before");
        hu_cont = hu_cont_fab.const_array();
    }
    // Layer thickness at v-points [H ~> m or kg m-2]; may be unallocated
    amrex::FArrayBox hv_cont_fab;
    amrex::Array4<const amrex::Real> hv_cont{};
    if (captured.logical("_hv_cont_allocated")) {
        hv_cont_fab = captured.fab_device("_hv_cont_before");
        hv_cont = hv_cont_fab.const_array();
    }

    const auto u                       = captured.fab_device("_u");
    const auto v                       = captured.fab_device("_v");
    const auto h                       = captured.fab_device("_h");
    const auto uh                      = captured.fab_device("_uh");
    const auto vh                      = captured.fab_device("_vh");
    const auto u_smooth                = captured.fab_device("_u_smooth");
    const auto v_smooth                = captured.fab_device("_v_smooth");
    const auto sh_xx_bt                = captured.fab_device("_sh_xx_bt");
    const auto sh_xy_bt                = captured.fab_device("_sh_xy_bt");
    const auto GME_effic_h             = captured.fab_device("_GME_effic_h");
    const auto GME_effic_q             = captured.fab_device("_GME_effic_q");
    const auto KH_u_GME                = captured.fab_device("_KH_u_GME");
    const auto KH_v_GME                = captured.fab_device("_KH_v_GME");
    auto       diffu                   = captured.fab_device("_diffu_before");
    auto       diffv                   = captured.fab_device("_diffv_before");
    auto       GME_coeff_h             = captured.fab_device("_GME_coeff_h_before");
    auto       GME_coeff_q             = captured.fab_device("_GME_coeff_q_before");
    auto       str_xx_GME              = captured.fab_device("_str_xx_GME_before");
    auto       str_xy_GME              = captured.fab_device("_str_xy_GME_before");
    auto       NoSt                    = captured.fab_device("_NoSt_before");
    auto       ShSt                    = captured.fab_device("_ShSt_before");
    auto       FrictWork               = captured.fab_device("_FrictWork_before");
    auto       FrictWork_bh            = captured.fab_device("_FrictWork_bh_before");
    auto       FrictWork_GME           = captured.fab_device("_FrictWork_GME_before");
    auto       Ah_h                    = captured.fab_device("_Ah_h_before");
    auto       Ah_q                    = captured.fab_device("_Ah_q_before");
    auto       Kh_h                    = captured.fab_device("_Kh_h_before");
    auto       Kh_q                    = captured.fab_device("_Kh_q_before");
    auto       grid_Re_Kh              = captured.fab_device("_grid_Re_Kh_before");
    auto       grid_Re_Ah              = captured.fab_device("_grid_Re_Ah_before");
    auto       div_xx_h                = captured.fab_device("_div_xx_h_before");
    auto       sh_xx_h                 = captured.fab_device("_sh_xx_h_before");
    auto       vort_xy_q               = captured.fab_device("_vort_xy_q_before");
    auto       sh_xy_q                 = captured.fab_device("_sh_xy_q_before");
    auto       visc_limit_h            = captured.fab_device("_visc_limit_h_before");
    auto       visc_limit_q            = captured.fab_device("_visc_limit_q_before");
    auto       visc_limit_h_flag       = captured.fab_device("_visc_limit_h_flag_before");
    auto       visc_limit_q_flag       = captured.fab_device("_visc_limit_q_flag_before");
    auto       visc_limit_h_frac       = captured.fab_device("_visc_limit_h_frac_before");
    auto       visc_limit_q_frac       = captured.fab_device("_visc_limit_q_frac_before");
    auto       BS_coeff_h              = captured.fab_device("_BS_coeff_h_before");
    auto       BS_coeff_q              = captured.fab_device("_BS_coeff_q_before");
    const auto diffu_after             = captured.fab_host("_diffu_after");
    const auto diffv_after             = captured.fab_host("_diffv_after");
    const auto GME_coeff_h_after       = captured.fab_host("_GME_coeff_h_after");
    const auto GME_coeff_q_after       = captured.fab_host("_GME_coeff_q_after");
    const auto str_xx_GME_after        = captured.fab_host("_str_xx_GME_after");
    const auto str_xy_GME_after        = captured.fab_host("_str_xy_GME_after");
    const auto NoSt_after              = captured.fab_host("_NoSt_after");
    const auto ShSt_after              = captured.fab_host("_ShSt_after");
    const auto FrictWork_after         = captured.fab_host("_FrictWork_after");
    const auto FrictWork_bh_after      = captured.fab_host("_FrictWork_bh_after");
    const auto FrictWork_GME_after     = captured.fab_host("_FrictWork_GME_after");
    const auto Ah_h_after              = captured.fab_host("_Ah_h_after");
    const auto Ah_q_after              = captured.fab_host("_Ah_q_after");
    const auto Kh_h_after              = captured.fab_host("_Kh_h_after");
    const auto Kh_q_after              = captured.fab_host("_Kh_q_after");
    const auto grid_Re_Kh_after        = captured.fab_host("_grid_Re_Kh_after");
    const auto grid_Re_Ah_after        = captured.fab_host("_grid_Re_Ah_after");
    const auto div_xx_h_after          = captured.fab_host("_div_xx_h_after");
    const auto sh_xx_h_after           = captured.fab_host("_sh_xx_h_after");
    const auto vort_xy_q_after         = captured.fab_host("_vort_xy_q_after");
    const auto sh_xy_q_after           = captured.fab_host("_sh_xy_q_after");
    const auto visc_limit_h_after      = captured.fab_host("_visc_limit_h_after");
    const auto visc_limit_q_after      = captured.fab_host("_visc_limit_q_after");
    const auto visc_limit_h_flag_after = captured.fab_host("_visc_limit_h_flag_after");
    const auto visc_limit_q_flag_after = captured.fab_host("_visc_limit_q_flag_after");
    const auto visc_limit_h_frac_after = captured.fab_host("_visc_limit_h_frac_after");
    const auto visc_limit_q_frac_after = captured.fab_host("_visc_limit_q_frac_after");
    const auto BS_coeff_h_after        = captured.fab_host("_BS_coeff_h_after");
    const auto BS_coeff_q_after        = captured.fab_host("_BS_coeff_q_after");

    MOM::horizontal_viscosity_tr(bxH0, bxQ0, CS, G, M, X, nullptr, nkblock, H_subroundoff, H_to_RZ, s_to_T,
                                 u.const_array(), v.const_array(), h.const_array(), uh.const_array(),
                                 vh.const_array(), diffu.array(), diffv.array(), hu_cont, hv_cont,
                                 u_smooth.const_array(), v_smooth.const_array(), sh_xx_bt.const_array(),
                                 sh_xy_bt.const_array(), GME_effic_h.const_array(), GME_effic_q.const_array(),
                                 KH_u_GME.const_array(), KH_v_GME.const_array(), GME_coeff_h.array(),
                                 GME_coeff_q.array(), str_xx_GME.array(), str_xy_GME.array(), NoSt.array(),
                                 ShSt.array(), FrictWork.array(), FrictWork_bh.array(), FrictWork_GME.array(),
                                 Ah_h.array(), Ah_q.array(), Kh_h.array(), Kh_q.array(), grid_Re_Kh.array(),
                                 grid_Re_Ah.array(), div_xx_h.array(), sh_xx_h.array(), vort_xy_q.array(),
                                 sh_xy_q.array(), visc_limit_h.array(), visc_limit_q.array(),
                                 visc_limit_h_flag.array(), visc_limit_q_flag.array(), visc_limit_h_frac.array(),
                                 visc_limit_q_frac.array(), BS_coeff_h.array(), BS_coeff_q.array());
    amrex::Gpu::synchronize();

    expect_arrays_equal(diffu_after,             to_host_fab(diffu),             "diffu");
    expect_arrays_equal(diffv_after,             to_host_fab(diffv),             "diffv");
    expect_arrays_equal(GME_coeff_h_after,       to_host_fab(GME_coeff_h),       "GME_coeff_h");
    expect_arrays_equal(GME_coeff_q_after,       to_host_fab(GME_coeff_q),       "GME_coeff_q");
    expect_arrays_equal(str_xx_GME_after,        to_host_fab(str_xx_GME),        "str_xx_GME");
    expect_arrays_equal(str_xy_GME_after,        to_host_fab(str_xy_GME),        "str_xy_GME");
    expect_arrays_equal(NoSt_after,              to_host_fab(NoSt),              "NoSt");
    expect_arrays_equal(ShSt_after,              to_host_fab(ShSt),              "ShSt");
    expect_arrays_equal(FrictWork_after,         to_host_fab(FrictWork),         "FrictWork");
    expect_arrays_equal(FrictWork_bh_after,      to_host_fab(FrictWork_bh),      "FrictWork_bh");
    expect_arrays_equal(FrictWork_GME_after,     to_host_fab(FrictWork_GME),     "FrictWork_GME");
    expect_arrays_equal(Ah_h_after,              to_host_fab(Ah_h),              "Ah_h");
    expect_arrays_equal(Ah_q_after,              to_host_fab(Ah_q),              "Ah_q");
    expect_arrays_equal(Kh_h_after,              to_host_fab(Kh_h),              "Kh_h");
    expect_arrays_equal(Kh_q_after,              to_host_fab(Kh_q),              "Kh_q");
    expect_arrays_equal(grid_Re_Kh_after,        to_host_fab(grid_Re_Kh),        "grid_Re_Kh");
    expect_arrays_equal(grid_Re_Ah_after,        to_host_fab(grid_Re_Ah),        "grid_Re_Ah");
    expect_arrays_equal(div_xx_h_after,          to_host_fab(div_xx_h),          "div_xx_h");
    expect_arrays_equal(sh_xx_h_after,           to_host_fab(sh_xx_h),           "sh_xx_h");
    expect_arrays_equal(vort_xy_q_after,         to_host_fab(vort_xy_q),         "vort_xy_q");
    expect_arrays_equal(sh_xy_q_after,           to_host_fab(sh_xy_q),           "sh_xy_q");
    expect_arrays_equal(visc_limit_h_after,      to_host_fab(visc_limit_h),      "visc_limit_h");
    expect_arrays_equal(visc_limit_q_after,      to_host_fab(visc_limit_q),      "visc_limit_q");
    expect_arrays_equal(visc_limit_h_flag_after, to_host_fab(visc_limit_h_flag), "visc_limit_h_flag");
    expect_arrays_equal(visc_limit_q_flag_after, to_host_fab(visc_limit_q_flag), "visc_limit_q_flag");
    expect_arrays_equal(visc_limit_h_frac_after, to_host_fab(visc_limit_h_frac), "visc_limit_h_frac");
    expect_arrays_equal(visc_limit_q_frac_after, to_host_fab(visc_limit_q_frac), "visc_limit_q_frac");
    expect_arrays_equal(BS_coeff_h_after,        to_host_fab(BS_coeff_h),        "BS_coeff_h");
    expect_arrays_equal(BS_coeff_q_after,        to_host_fab(BS_coeff_q),        "BS_coeff_q");
    if (MEKE_mom_src_allocated) {
        expect_arrays_equal(captured.fab_host("_MEKE_mom_src_after"), to_host_fab(MEKE_mom_src_fab), "MEKE_mom_src");
    }
    if (MEKE_mom_src_bh_allocated) {
        expect_arrays_equal(captured.fab_host("_MEKE_mom_src_bh_after"), to_host_fab(MEKE_mom_src_bh_fab),
                            "MEKE_mom_src_bh");
    }
    if (MEKE_GME_snk_allocated) {
        expect_arrays_equal(captured.fab_host("_MEKE_GME_snk_after"), to_host_fab(MEKE_GME_snk_fab), "MEKE_GME_snk");
    }
    if (skeb_diss_allocated) {
        expect_arrays_equal(captured.fab_host("_STOCH_skeb_diss_after"), to_host_fab(skeb_diss_fab), "skeb_diss");
    }
}
