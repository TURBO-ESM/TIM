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
