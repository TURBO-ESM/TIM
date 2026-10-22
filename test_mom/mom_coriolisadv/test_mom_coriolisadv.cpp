// Unit tests for MOM::gradke / coradv_setup / coradv_common_terms / coradv_arakawa /
// coradv_robust_enstro / coradv_sadourny / coradv_weno / coradcalc_tr.
//
// Each test loads a captured Fortran (input, expected-output) pair from
// <data-dir>/<name>.{bin,meta}, runs the C++ kernel over equivalent AMReX
// containers, and compares the result against the captured "after" arrays.

// SKILLS: 0.3.2

#include <gtest/gtest.h>

#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>

#include "amrex_assertions.hpp"
#include "captured_io.hpp"
#include "data_dir.hpp"
#include "mom_coriolisadv.hpp"

using test_mom::expect_arrays_equal;
using test_mom::to_host_fab;

namespace {

// An acceleration diagnostic that is only present when it was allocated at
// capture time ("<prefix>_allocated"). When present, its pre- and post-call
// states are "<prefix>_before" / "<prefix>_after"; when absent, `arr` stays
// default-constructed (.p == nullptr).
struct OptionalDiagnostic {
    amrex::FArrayBox before_fab;
    amrex::FArrayBox after_fab;
    amrex::Array4<amrex::Real> arr{};
    bool present = false;
};

OptionalDiagnostic bind_optional_diagnostic(const test_mom::CapturedFile& captured,
                                            const std::string& prefix) {
    OptionalDiagnostic d;
    d.present = captured.logical(prefix + "_allocated");
    if (d.present) {
        d.before_fab = captured.fab_device(prefix + "_before");
        d.arr = d.before_fab.array();
        d.after_fab = captured.fab_host(prefix + "_after");
    }
    return d;
}

} // namespace

// -------------------------------------------------------------------------
// gradKE
// -------------------------------------------------------------------------
// Field mapping: MOM6 src/core/MOM_CoriolisAdv.F90:3843–3868 (gradKE capture block).
TEST(Gradke, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "gradke");

    const auto   bxQ            = captured.box("_bxQ");
    const auto   bxU            = captured.box("_bxU");
    const auto   bxV            = captured.box("_bxV");
    const auto   u              = captured.fab_device("_u");
    const auto   v              = captured.fab_device("_v");
    auto         KE             = captured.fab_device("_KE_before");
    auto         KEx            = captured.fab_device("_KEx_before");
    auto         KEy            = captured.fab_device("_KEy_before");
    const auto   areaCu         = captured.fab_device("_areaCu");
    const auto   areaCv         = captured.fab_device("_areaCv");
    const auto   IareaT         = captured.fab_device("_IareaT");
    const auto   mask2dCu       = captured.fab_device("_mask2dCu");
    const auto   mask2dCv       = captured.fab_device("_mask2dCv");
    const auto   IdxCu_OBCmask  = captured.fab_device("_IdxCu_OBCmask");
    const auto   IdyCv_OBCmask  = captured.fab_device("_IdyCv_OBCmask");
    const auto   KE_after       = captured.fab_host("_KE_after");
    const auto   KEx_after      = captured.fab_host("_KEx_after");
    const auto   KEy_after      = captured.fab_host("_KEy_after");
    const int    KE_Scheme      = captured.integer("_KE_Scheme");
    const bool   KE_use_limiter = captured.logical("_KE_use_limiter");

    MOM::gradke(bxQ, bxU, bxV,
                u.const_array(),
                v.const_array(),
                KE.array(),
                KEx.array(),
                KEy.array(),
                areaCu.const_array(),
                areaCv.const_array(),
                IareaT.const_array(),
                mask2dCu.const_array(),
                mask2dCv.const_array(),
                IdxCu_OBCmask.const_array(),
                IdyCv_OBCmask.const_array(),
                KE_Scheme,
                KE_use_limiter);
    amrex::Gpu::synchronize();

    expect_arrays_equal(KE_after,  to_host_fab(KE),  "KE");
    expect_arrays_equal(KEx_after, to_host_fab(KEx), "KEx");
    expect_arrays_equal(KEy_after, to_host_fab(KEy), "KEy");
}

// -------------------------------------------------------------------------
// CorAdv_setup
// -------------------------------------------------------------------------
// Field mapping: MOM6 src/core/MOM_CoriolisAdv.F90:1869-1924 (CorAdv_setup capture block).
TEST(CoradvSetup, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "coradv_setup");

    const auto   bxH               = captured.box("_bxH");
    const auto   bxQ               = captured.box("_bxQ");
    const auto   bxU               = captured.box("_bxU");
    const auto   bxV               = captured.box("_bxV");
    const auto   bxQs              = captured.box("_bxQs");
    const auto   u                 = captured.fab_device("_u");
    const auto   v                 = captured.fab_device("_v");
    const auto   h                 = captured.fab_device("_h");
    const bool   rv_x_v_associated = captured.logical("_AD_rv_x_v_allocated");
    const bool   rv_x_u_associated = captured.logical("_AD_rv_x_u_allocated");
    const auto   por_face_areaU    = captured.fab_device("_por_face_areaU");
    const auto   por_face_areaV    = captured.fab_device("_por_face_areaV");
    const auto   dyCv              = captured.fab_device("_dyCv");
    const auto   dxCu              = captured.fab_device("_dxCu");
    const auto   dy_Cu             = captured.fab_device("_dy_Cu");
    const auto   dx_Cv             = captured.fab_device("_dx_Cv");
    const auto   dyBu              = captured.fab_device("_dyBu");
    const auto   dxBu              = captured.fab_device("_dxBu");
    const auto   mask2dBu          = captured.fab_device("_mask2dBu");
    const auto   IareaBu           = captured.fab_device("_IareaBu");
    const auto   CoriolisBu        = captured.fab_device("_CoriolisBu");
    const auto   Area_h            = captured.fab_device("_Area_h");
    const auto   Area_q            = captured.fab_device("_Area_q");
    const bool   Stokes_VF         = captured.logical("_Stokes_VF");
    const bool   use_weno          = captured.logical("_use_weno");
    const double vol_neglect       = captured.real64("_vol_neglect");
    const double area_neglect      = captured.real64("_area_neglect");
    auto         q                 = captured.fab_device("_q_before");
    auto         Ih_q              = captured.fab_device("_Ih_q_before");
    auto         abs_vort          = captured.fab_device("_abs_vort_before");
    auto         h_q               = captured.fab_device("_h_q_before");
    auto         qS                = captured.fab_device("_qS_before");
    auto         q2                = captured.fab_device("_q2_before");
    auto         uh_center         = captured.fab_device("_uh_center_before");
    auto         vh_center         = captured.fab_device("_vh_center_before");
    auto         RV                = captured.fab_device("_RV_before");
    auto         PV                = captured.fab_device("_PV_before");
    const auto   q_after           = captured.fab_host("_q_after");
    const auto   Ih_q_after        = captured.fab_host("_Ih_q_after");
    const auto   abs_vort_after    = captured.fab_host("_abs_vort_after");
    const auto   h_q_after         = captured.fab_host("_h_q_after");
    const auto   qS_after          = captured.fab_host("_qS_after");
    const auto   q2_after          = captured.fab_host("_q2_after");
    const auto   uh_center_after   = captured.fab_host("_uh_center_after");
    const auto   vh_center_after   = captured.fab_host("_vh_center_after");
    const auto   RV_after          = captured.fab_host("_RV_after");
    const auto   PV_after          = captured.fab_host("_PV_after");

    // Options used in the Coriolis and momentum advection calculations
    CorAdv_opts_CS_C CS{};
    CS.Coriolis_scheme_CS.Coriolis_Scheme = captured.integer("_CS_Coriolis_scheme_CS_Coriolis_Scheme");
    CS.Coriolis_scheme_CS.Coriolis_En_Dis = captured.logical("_CS_Coriolis_scheme_CS_Coriolis_En_Dis");
    CS.KE_Scheme            = captured.integer("_CS_KE_Scheme");
    CS.KE_use_limiter       = captured.logical("_CS_KE_use_limiter");
    CS.PV_Adv_Scheme        = captured.integer("_CS_PV_Adv_Scheme");
    CS.nkblock              = captured.integer("_CS_nkblock");
    CS.F_eff_max_blend      = captured.real64("_CS_F_eff_max_blend");
    CS.wt_lin_blend         = captured.real64("_CS_wt_lin_blend");
    CS.no_slip              = captured.logical("_CS_no_slip");
    CS.bound_Coriolis       = captured.logical("_CS_bound_Coriolis");
    CS.weno_velocity_smooth = captured.logical("_CS_weno_velocity_smooth");
    CS.do_CAS_diag          = captured.logical("_CS_do_CAS_diag");
    CS.do_RV_diag           = captured.logical("_CS_do_RV_diag");
    CS.do_PV_diag           = captured.logical("_CS_do_PV_diag");

    MOM::coradv_setup(bxH, bxQ, bxU, bxV, bxQs,
                      u.const_array(),
                      v.const_array(),
                      h.const_array(),
                      nullptr,
                      rv_x_v_associated,
                      rv_x_u_associated,
                      por_face_areaU.const_array(),
                      por_face_areaV.const_array(),
                      nullptr,
                      dyCv.const_array(),
                      dxCu.const_array(),
                      dy_Cu.const_array(),
                      dx_Cv.const_array(),
                      dyBu.const_array(),
                      dxBu.const_array(),
                      mask2dBu.const_array(),
                      IareaBu.const_array(),
                      CoriolisBu.const_array(),
                      Area_h.const_array(),
                      Area_q.const_array(),
                      Stokes_VF,
                      use_weno,
                      vol_neglect,
                      area_neglect,
                      q.array(),
                      Ih_q.array(),
                      abs_vort.array(),
                      h_q.array(),
                      qS.array(),
                      q2.array(),
                      uh_center.array(),
                      vh_center.array(),
                      RV.array(),
                      PV.array(),
                      CS);
    amrex::Gpu::synchronize();

    expect_arrays_equal(q_after,           to_host_fab(q),         "q");
    expect_arrays_equal(Ih_q_after,        to_host_fab(Ih_q),      "Ih_q");
    expect_arrays_equal(abs_vort_after,    to_host_fab(abs_vort),  "abs_vort");
    expect_arrays_equal(h_q_after,         to_host_fab(h_q),       "h_q");
    expect_arrays_equal(qS_after,          to_host_fab(qS),        "qS");
    expect_arrays_equal(q2_after,          to_host_fab(q2),        "q2");
    expect_arrays_equal(uh_center_after,   to_host_fab(uh_center), "uh_center");
    expect_arrays_equal(vh_center_after,   to_host_fab(vh_center), "vh_center");
    expect_arrays_equal(RV_after,          to_host_fab(RV),        "RV");
    expect_arrays_equal(PV_after,          to_host_fab(PV),        "PV");
}

// -------------------------------------------------------------------------
// CorAdv_common_terms
// -------------------------------------------------------------------------
// Field mapping: MOM6 src/core/MOM_CoriolisAdv.F90:3472-3504 (CorAdv_common_terms capture block).
TEST(CoradvCommonTerms, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "coradv_common_terms");

    const auto   bxU            = captured.box("_bxU");
    const auto   bxV            = captured.box("_bxV");
    const auto   u              = captured.fab_device("_u");
    const auto   v              = captured.fab_device("_v");
    const auto   uh             = captured.fab_device("_uh");
    const auto   vh             = captured.fab_device("_vh");
    const auto   abs_vort       = captured.fab_device("_abs_vort");
    const auto   qS             = captured.fab_device("_qS");
    const auto   q2             = captured.fab_device("_q2");
    const auto   KEx            = captured.fab_device("_KEx");
    const auto   KEy            = captured.fab_device("_KEy");
    const bool   Stokes_VF      = captured.logical("_Stokes_VF");
    const auto   IdxCu          = captured.fab_device("_IdxCu");
    const auto   IdyCv          = captured.fab_device("_IdyCv");
    auto         CAu            = captured.fab_device("_CAu_before");
    auto         CAv            = captured.fab_device("_CAv_before");
    auto         CAuS           = captured.fab_device("_CAuS_before");
    auto         CAvS           = captured.fab_device("_CAvS_before");
    const auto   CAu_after      = captured.fab_host("_CAu_after");
    const auto   CAv_after      = captured.fab_host("_CAv_after");
    const auto   CAuS_after     = captured.fab_host("_CAuS_after");
    const auto   CAvS_after     = captured.fab_host("_CAvS_after");
    auto         gradKEu        = bind_optional_diagnostic(captured, "_AD_gradKEu");
    auto         gradKEv        = bind_optional_diagnostic(captured, "_AD_gradKEv");
    auto         rv_x_v         = bind_optional_diagnostic(captured, "_AD_rv_x_v");
    auto         rv_x_u         = bind_optional_diagnostic(captured, "_AD_rv_x_u");

    // Options used in the Coriolis and momentum advection calculations
    CorAdv_opts_CS_C CS{};
    CS.Coriolis_scheme_CS.Coriolis_Scheme = captured.integer("_CS_Coriolis_scheme_CS_Coriolis_Scheme");
    CS.Coriolis_scheme_CS.Coriolis_En_Dis = captured.logical("_CS_Coriolis_scheme_CS_Coriolis_En_Dis");
    CS.KE_Scheme            = captured.integer("_CS_KE_Scheme");
    CS.KE_use_limiter       = captured.logical("_CS_KE_use_limiter");
    CS.PV_Adv_Scheme        = captured.integer("_CS_PV_Adv_Scheme");
    CS.nkblock              = captured.integer("_CS_nkblock");
    CS.F_eff_max_blend      = captured.real64("_CS_F_eff_max_blend");
    CS.wt_lin_blend         = captured.real64("_CS_wt_lin_blend");
    CS.no_slip              = captured.logical("_CS_no_slip");
    CS.bound_Coriolis       = captured.logical("_CS_bound_Coriolis");
    CS.weno_velocity_smooth = captured.logical("_CS_weno_velocity_smooth");
    CS.do_CAS_diag          = captured.logical("_CS_do_CAS_diag");
    CS.do_RV_diag           = captured.logical("_CS_do_RV_diag");
    CS.do_PV_diag           = captured.logical("_CS_do_PV_diag");

    MOM::coradv_common_terms(bxU, bxV,
                             u.const_array(),
                             v.const_array(),
                             uh.const_array(),
                             vh.const_array(),
                             abs_vort.const_array(),
                             qS.const_array(),
                             q2.const_array(),
                             KEx.const_array(),
                             KEy.const_array(),
                             Stokes_VF,
                             IdxCu.const_array(),
                             IdyCv.const_array(),
                             CAu.array(),
                             CAv.array(),
                             CAuS.array(),
                             CAvS.array(),
                             gradKEu.arr,
                             gradKEv.arr,
                             rv_x_v.arr,
                             rv_x_u.arr,
                             CS);
    amrex::Gpu::synchronize();

    expect_arrays_equal(CAu_after,  to_host_fab(CAu),  "CAu");
    expect_arrays_equal(CAv_after,  to_host_fab(CAv),  "CAv");
    expect_arrays_equal(CAuS_after, to_host_fab(CAuS), "CAuS");
    expect_arrays_equal(CAvS_after, to_host_fab(CAvS), "CAvS");
    if (gradKEu.present) {
        expect_arrays_equal(gradKEu.after_fab, to_host_fab(gradKEu.before_fab), "gradKEu");
    }
    if (gradKEv.present) {
        expect_arrays_equal(gradKEv.after_fab, to_host_fab(gradKEv.before_fab), "gradKEv");
    }
    if (rv_x_v.present) {
        expect_arrays_equal(rv_x_v.after_fab, to_host_fab(rv_x_v.before_fab), "rv_x_v");
    }
    if (rv_x_u.present) {
        expect_arrays_equal(rv_x_u.after_fab, to_host_fab(rv_x_u.before_fab), "rv_x_u");
    }
}

// -------------------------------------------------------------------------
// CorAdv_arakawa
// -------------------------------------------------------------------------
// Field mapping: MOM6 src/core/MOM_CoriolisAdv.F90:2545-2568 (CorAdv_arakawa capture block).
TEST(CoradvArakawa, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "coradv_arakawa");

    const auto   bxH          = captured.box("_bxH");
    const auto   bxQ          = captured.box("_bxQ");
    const auto   bxU          = captured.box("_bxU");
    const auto   bxV          = captured.box("_bxV");
    const auto   uh           = captured.fab_device("_uh");
    const auto   vh           = captured.fab_device("_vh");
    const auto   q            = captured.fab_device("_q");
    const auto   Ih_q         = captured.fab_device("_Ih_q");
    const double Fe_m2        = captured.real64("_Fe_m2");
    const double rat_lin      = captured.real64("_rat_lin");
    const auto   IdxCu        = captured.fab_device("_IdxCu");
    const auto   IdyCv        = captured.fab_device("_IdyCv");
    auto         CAu          = captured.fab_device("_CAu_before");
    auto         CAv          = captured.fab_device("_CAv_before");
    const auto   CAu_after    = captured.fab_host("_CAu_after");
    const auto   CAv_after    = captured.fab_host("_CAv_after");
    const double wt_lin_blend = captured.real64("_wt_lin_blend");

    // Options selecting the Coriolis discretization
    Coriolis_scheme_CS_C CS{};
    CS.Coriolis_Scheme = captured.integer("_CS_Coriolis_Scheme");
    CS.Coriolis_En_Dis = captured.logical("_CS_Coriolis_En_Dis");

    MOM::coradv_arakawa(bxH, bxQ, bxU, bxV,
                        uh.const_array(),
                        vh.const_array(),
                        q.const_array(),
                        Ih_q.const_array(),
                        Fe_m2,
                        rat_lin,
                        IdxCu.const_array(),
                        IdyCv.const_array(),
                        CAu.array(),
                        CAv.array(),
                        wt_lin_blend,
                        CS);
    amrex::Gpu::synchronize();

    expect_arrays_equal(CAu_after, to_host_fab(CAu), "CAu");
    expect_arrays_equal(CAv_after, to_host_fab(CAv), "CAv");
}

// -------------------------------------------------------------------------
// CorAdv_robust_enstro
// -------------------------------------------------------------------------
// Field mapping: MOM6 src/core/MOM_CoriolisAdv.F90:2760-2784 (CorAdv_robust_enstro capture block).
TEST(CoradvRobustEnstro, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "coradv_robust_enstro");

    const auto   bxU           = captured.box("_bxU");
    const auto   bxV           = captured.box("_bxV");
    const auto   u             = captured.fab_device("_u");
    const auto   v             = captured.fab_device("_v");
    const auto   h             = captured.fab_device("_h");
    const auto   uh            = captured.fab_device("_uh");
    const auto   vh            = captured.fab_device("_vh");
    const auto   abs_vort      = captured.fab_device("_abs_vort");
    const double eps_vel       = captured.real64("_eps_vel");
    const double h_tiny        = captured.real64("_h_tiny");
    const auto   IdxCv         = captured.fab_device("_IdxCv");
    const auto   IdxCu         = captured.fab_device("_IdxCu");
    const auto   IdyCu         = captured.fab_device("_IdyCu");
    const auto   IdyCv         = captured.fab_device("_IdyCv");
    auto         CAu           = captured.fab_device("_CAu_before");
    auto         CAv           = captured.fab_device("_CAv_before");
    const auto   CAu_after     = captured.fab_host("_CAu_after");
    const auto   CAv_after     = captured.fab_host("_CAv_after");
    const int    PV_Adv_Scheme = captured.integer("_PV_Adv_Scheme");

    MOM::coradv_robust_enstro(bxU, bxV,
                              u.const_array(),
                              v.const_array(),
                              h.const_array(),
                              uh.const_array(),
                              vh.const_array(),
                              abs_vort.const_array(),
                              eps_vel,
                              h_tiny,
                              IdxCv.const_array(),
                              IdxCu.const_array(),
                              IdyCu.const_array(),
                              IdyCv.const_array(),
                              CAu.array(),
                              CAv.array(),
                              PV_Adv_Scheme);
    amrex::Gpu::synchronize();

    expect_arrays_equal(CAu_after, to_host_fab(CAu), "CAu");
    expect_arrays_equal(CAv_after, to_host_fab(CAv), "CAv");
}

// -------------------------------------------------------------------------
// CorAdv_sadourny
// -------------------------------------------------------------------------
// Field mapping: MOM6 src/core/MOM_CoriolisAdv.F90:2246-2270 (CorAdv_sadourny capture block).
TEST(CoradvSadourny, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "coradv_sadourny");

    const auto   bxH       = captured.box("_bxH");
    const auto   bxU       = captured.box("_bxU");
    const auto   bxV       = captured.box("_bxV");
    const auto   u         = captured.fab_device("_u");
    const auto   v         = captured.fab_device("_v");
    const auto   uh        = captured.fab_device("_uh");
    const auto   vh        = captured.fab_device("_vh");
    const auto   q         = captured.fab_device("_q");
    const auto   uh_center = captured.fab_device("_uh_center");
    const auto   vh_center = captured.fab_device("_vh_center");
    const auto   dy_Cu     = captured.fab_device("_dy_Cu");
    const auto   dx_Cv     = captured.fab_device("_dx_Cv");
    const auto   IdxCu     = captured.fab_device("_IdxCu");
    const auto   IdyCv     = captured.fab_device("_IdyCv");
    auto         CAu       = captured.fab_device("_CAu_before");
    auto         CAv       = captured.fab_device("_CAv_before");
    const auto   CAu_after = captured.fab_host("_CAu_after");
    const auto   CAv_after = captured.fab_host("_CAv_after");

    // Options selecting the Coriolis discretization
    Coriolis_scheme_CS_C CS{};
    CS.Coriolis_Scheme = captured.integer("_CS_Coriolis_Scheme");
    CS.Coriolis_En_Dis = captured.logical("_CS_Coriolis_En_Dis");

    MOM::coradv_sadourny(bxH, bxU, bxV,
                         u.const_array(),
                         v.const_array(),
                         uh.const_array(),
                         vh.const_array(),
                         q.const_array(),
                         uh_center.const_array(),
                         vh_center.const_array(),
                         dy_Cu.const_array(),
                         dx_Cv.const_array(),
                         IdxCu.const_array(),
                         IdyCv.const_array(),
                         CAu.array(),
                         CAv.array(),
                         CS);
    amrex::Gpu::synchronize();

    expect_arrays_equal(CAu_after, to_host_fab(CAu), "CAu");
    expect_arrays_equal(CAv_after, to_host_fab(CAv), "CAv");
}

// -------------------------------------------------------------------------
// CorAdv_weno
// -------------------------------------------------------------------------
// Field mapping: MOM6 src/core/MOM_CoriolisAdv.F90:3163-3188 (CorAdv_weno capture block).
TEST(CoradvWeno, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "coradv_weno");

    const auto   bxU                  = captured.box("_bxU");
    const auto   bxV                  = captured.box("_bxV");
    const auto   u                    = captured.fab_device("_u");
    const auto   v                    = captured.fab_device("_v");
    const auto   uh                   = captured.fab_device("_uh");
    const auto   vh                   = captured.fab_device("_vh");
    const auto   q                    = captured.fab_device("_q");
    const auto   abs_vort             = captured.fab_device("_abs_vort");
    const auto   h_q                  = captured.fab_device("_h_q");
    const auto   IdxCu                = captured.fab_device("_IdxCu");
    const auto   mask2dCu             = captured.fab_device("_mask2dCu");
    const auto   IdyCv                = captured.fab_device("_IdyCv");
    const auto   mask2dCv             = captured.fab_device("_mask2dCv");
    const double H_subroundoff        = captured.real64("_H_subroundoff");
    auto         CAu                  = captured.fab_device("_CAu_before");
    auto         CAv                  = captured.fab_device("_CAv_before");
    const auto   CAu_after            = captured.fab_host("_CAu_after");
    const auto   CAv_after            = captured.fab_host("_CAv_after");
    const bool   weno_velocity_smooth = captured.logical("_weno_velocity_smooth");

    // Options selecting the Coriolis discretization
    Coriolis_scheme_CS_C CS{};
    CS.Coriolis_Scheme = captured.integer("_CS_Coriolis_Scheme");
    CS.Coriolis_En_Dis = captured.logical("_CS_Coriolis_En_Dis");

    MOM::coradv_weno(bxU, bxV,
                     u.const_array(),
                     v.const_array(),
                     uh.const_array(),
                     vh.const_array(),
                     q.const_array(),
                     abs_vort.const_array(),
                     h_q.const_array(),
                     IdxCu.const_array(),
                     mask2dCu.const_array(),
                     IdyCv.const_array(),
                     mask2dCv.const_array(),
                     H_subroundoff,
                     CAu.array(),
                     CAv.array(),
                     weno_velocity_smooth,
                     CS);
    amrex::Gpu::synchronize();

    expect_arrays_equal(CAu_after, to_host_fab(CAu), "CAu");
    expect_arrays_equal(CAv_after, to_host_fab(CAv), "CAv");
}

// -------------------------------------------------------------------------
// CorAdCalc_TR
// -------------------------------------------------------------------------
// Field mapping: MOM6 src/core/MOM_CoriolisAdv.F90:1250-1310 (CorAdCalc_TR capture block).
TEST(CoradcalcTr, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "coradcalc_tr");

    const auto   bxH0              = captured.box("_bxH0");
    const auto   bxQ0              = captured.box("_bxQ0");
    const auto   u                 = captured.fab_device("_u");
    const auto   v                 = captured.fab_device("_v");
    const auto   h                 = captured.fab_device("_h");
    const auto   uh                = captured.fab_device("_uh");
    const auto   vh                = captured.fab_device("_vh");
    auto         CAu               = captured.fab_device("_CAu_before");
    auto         CAv               = captured.fab_device("_CAv_before");
    auto         RV                = captured.fab_device("_RV_before");
    auto         PV                = captured.fab_device("_PV_before");
    auto         CAuS              = captured.fab_device("_CAuS_before");
    auto         CAvS              = captured.fab_device("_CAvS_before");
    const auto   CAu_after         = captured.fab_host("_CAu_after");
    const auto   CAv_after         = captured.fab_host("_CAv_after");
    const auto   RV_after          = captured.fab_host("_RV_after");
    const auto   PV_after          = captured.fab_host("_PV_after");
    const auto   CAuS_after        = captured.fab_host("_CAuS_after");
    const auto   CAvS_after        = captured.fab_host("_CAvS_after");
    auto         gradKEu           = bind_optional_diagnostic(captured, "_AD_gradKEu");
    auto         gradKEv           = bind_optional_diagnostic(captured, "_AD_gradKEv");
    auto         rv_x_v            = bind_optional_diagnostic(captured, "_AD_rv_x_v");
    auto         rv_x_u            = bind_optional_diagnostic(captured, "_AD_rv_x_u");
    const auto   mask2dT           = captured.fab_device("_mask2dT");
    const auto   areaT             = captured.fab_device("_areaT");
    const auto   areaCu            = captured.fab_device("_areaCu");
    const auto   areaCv            = captured.fab_device("_areaCv");
    const auto   IareaT            = captured.fab_device("_IareaT");
    const auto   mask2dCu          = captured.fab_device("_mask2dCu");
    const auto   mask2dCv          = captured.fab_device("_mask2dCv");
    const auto   dyCv              = captured.fab_device("_dyCv");
    const auto   dxCu              = captured.fab_device("_dxCu");
    const auto   dy_Cu             = captured.fab_device("_dy_Cu");
    const auto   dx_Cv             = captured.fab_device("_dx_Cv");
    const auto   dyBu              = captured.fab_device("_dyBu");
    const auto   dxBu              = captured.fab_device("_dxBu");
    const auto   mask2dBu          = captured.fab_device("_mask2dBu");
    const auto   IareaBu           = captured.fab_device("_IareaBu");
    const auto   CoriolisBu        = captured.fab_device("_CoriolisBu");
    const auto   IdxCu             = captured.fab_device("_IdxCu");
    const auto   IdyCv             = captured.fab_device("_IdyCv");
    const auto   IdxCv             = captured.fab_device("_IdxCv");
    const auto   IdyCu             = captured.fab_device("_IdyCu");
    const auto   IdxCu_OBCmask     = captured.fab_device("_IdxCu_OBCmask");
    const auto   IdyCv_OBCmask     = captured.fab_device("_IdyCv_OBCmask");
    const int    nz                = captured.integer("_nz");
    const double H_subroundoff     = captured.real64("_H_subroundoff");
    const double Angstrom_H        = captured.real64("_Angstrom_H");
    const double m_to_L            = captured.real64("_m_to_L");
    const double m_s_to_L_T        = captured.real64("_m_s_to_L_T");
    const auto   por_face_areaU    = captured.fab_device("_por_face_areaU");
    const auto   por_face_areaV    = captured.fab_device("_por_face_areaV");
    const bool   Stokes_VF         = captured.logical("_Stokes_VF");

    // Options used in the Coriolis and momentum advection calculations
    CorAdv_opts_CS_C CS{};
    CS.Coriolis_scheme_CS.Coriolis_Scheme = captured.integer("_CS_Coriolis_scheme_CS_Coriolis_Scheme");
    CS.Coriolis_scheme_CS.Coriolis_En_Dis = captured.logical("_CS_Coriolis_scheme_CS_Coriolis_En_Dis");
    CS.KE_Scheme            = captured.integer("_CS_KE_Scheme");
    CS.KE_use_limiter       = captured.logical("_CS_KE_use_limiter");
    CS.PV_Adv_Scheme        = captured.integer("_CS_PV_Adv_Scheme");
    CS.nkblock              = captured.integer("_CS_nkblock");
    CS.F_eff_max_blend      = captured.real64("_CS_F_eff_max_blend");
    CS.wt_lin_blend         = captured.real64("_CS_wt_lin_blend");
    CS.no_slip              = captured.logical("_CS_no_slip");
    CS.bound_Coriolis       = captured.logical("_CS_bound_Coriolis");
    CS.weno_velocity_smooth = captured.logical("_CS_weno_velocity_smooth");
    CS.do_CAS_diag          = captured.logical("_CS_do_CAS_diag");
    CS.do_RV_diag           = captured.logical("_CS_do_RV_diag");
    CS.do_PV_diag           = captured.logical("_CS_do_PV_diag");

    MOM::coradcalc_tr(bxH0, bxQ0,
                      u.const_array(),
                      v.const_array(),
                      h.const_array(),
                      uh.const_array(),
                      vh.const_array(),
                      CAu.array(),
                      CAv.array(),
                      RV.array(),
                      PV.array(),
                      CAuS.array(),
                      CAvS.array(),
                      nullptr,
                      gradKEu.arr,
                      gradKEv.arr,
                      rv_x_v.arr,
                      rv_x_u.arr,
                      mask2dT.const_array(),
                      areaT.const_array(),
                      areaCu.const_array(),
                      areaCv.const_array(),
                      IareaT.const_array(),
                      mask2dCu.const_array(),
                      mask2dCv.const_array(),
                      dyCv.const_array(),
                      dxCu.const_array(),
                      dy_Cu.const_array(),
                      dx_Cv.const_array(),
                      dyBu.const_array(),
                      dxBu.const_array(),
                      mask2dBu.const_array(),
                      IareaBu.const_array(),
                      CoriolisBu.const_array(),
                      IdxCu.const_array(),
                      IdyCv.const_array(),
                      IdxCv.const_array(),
                      IdyCu.const_array(),
                      IdxCu_OBCmask.const_array(),
                      IdyCv_OBCmask.const_array(),
                      nz,
                      H_subroundoff,
                      Angstrom_H,
                      m_to_L,
                      m_s_to_L_T,
                      CS,
                      por_face_areaU.const_array(),
                      por_face_areaV.const_array(),
                      nullptr,
                      Stokes_VF);
    amrex::Gpu::synchronize();

    expect_arrays_equal(CAu_after,  to_host_fab(CAu),  "CAu");
    expect_arrays_equal(CAv_after,  to_host_fab(CAv),  "CAv");
    expect_arrays_equal(RV_after,   to_host_fab(RV),   "RV");
    expect_arrays_equal(PV_after,   to_host_fab(PV),   "PV");
    expect_arrays_equal(CAuS_after, to_host_fab(CAuS), "CAuS");
    expect_arrays_equal(CAvS_after, to_host_fab(CAvS), "CAvS");
    if (gradKEu.present) {
        expect_arrays_equal(gradKEu.after_fab, to_host_fab(gradKEu.before_fab), "gradKEu");
    }
    if (gradKEv.present) {
        expect_arrays_equal(gradKEv.after_fab, to_host_fab(gradKEv.before_fab), "gradKEv");
    }
    if (rv_x_v.present) {
        expect_arrays_equal(rv_x_v.after_fab, to_host_fab(rv_x_v.before_fab), "rv_x_v");
    }
    if (rv_x_u.present) {
        expect_arrays_equal(rv_x_u.after_fab, to_host_fab(rv_x_u.before_fab), "rv_x_u");
    }
}
