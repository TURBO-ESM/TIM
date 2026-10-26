// Unit tests for MOM::smooth_x9_h / smooth_x9_uv / smooth_gme_h_pass / smooth_gme_q_pass /
// hor_visc_backscatter_h / hor_visc_backscatter_q / hor_visc_gme_bt / hor_visc_leith_grad_a /
// hor_visc_leith_grad_b / hor_visc_leithy_m / hor_visc_leithy_ah.
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

    const auto bxQ = captured.box("_bxQ");

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

    MOM::hor_visc_backscatter_h(bxQ, CS,
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

    const auto bxH = captured.box("_bxH");

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

    MOM::hor_visc_backscatter_q(bxH, CS,
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

    const auto bxH    = captured.box("_bxH");
    const auto bxKh   = captured.box("_bxKh");
    const auto bxKh_i = captured.box("_bxKh_i");
    const auto bxKh_j = captured.box("_bxKh_j");

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

    MOM::hor_visc_leith_grad_a(bxH, bxKh, bxKh_i, bxKh_j, CS,
                               dyBu.const_array(), IdxBu.const_array(), IdyCu.const_array(), dxBu.const_array(), IdyBu.const_array(),
                               IdxCv.const_array(), IdyCv.const_array(), IdxCu.const_array(), dF_dx.const_array(), dF_dy.const_array(),
                               dudx.const_array(), dvdy.const_array(), vort_xy.const_array(), vort_xy_smooth.const_array(),
                               grad_vort_mag_h_2d.array(), grad_div_mag_h.array(), grad_vort_mag_q_2d.array(), grad_div_mag_q.array(),
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

    const auto bxH  = captured.box("_bxH");
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

    MOM::hor_visc_leith_grad_b(bxH, bxKh, CS,
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
