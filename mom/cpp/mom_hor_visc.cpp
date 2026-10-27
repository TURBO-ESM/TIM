// SKILLS: 0.3.3
///
/// @file mom_hor_visc.cpp
/// @brief Box-level AMReX kernels for the MOM6 horizontal viscosity terms.
///
#include <AMReX.H>
#include <AMReX_Algorithm.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <AMReX_Math.H>

#include "mom_hor_visc.hpp"

#include <cmath>
#include <string>

/// @brief Abort with @p msg annotated by the source file and line number.
/// @param msg Message to report before aborting.
#define AMREX_ABORT_LOC(msg) \
	amrex::Abort(std::string(msg) + " [" + __FILE__ + ":" + std::to_string(__LINE__) + "]")

using namespace amrex;
using namespace amrex::literals;

namespace MOM {

///
/// @brief Grow the i and j extents of a box by separate amounts on its low and high sides.
///
/// @param bx   The box to grow
/// @param lo_i Number of cells to grow the low end in i
/// @param lo_j Number of cells to grow the low end in j
/// @param hi_i Number of cells to grow the high end in i
/// @param hi_j Number of cells to grow the high end in j
///
static Box grow_by(Box bx, const int lo_i, const int lo_j, const int hi_i, const int hi_j)
{
    bx.growLo(0, lo_i).growLo(1, lo_j).growHi(0, hi_i).growHi(1, hi_j);
    return bx;
}

///
/// @brief One pass of the 9-point smoothing filter applied to a field.
///
/// @tparam zero_land If false, return the average of the surrounding ocean points when
///                   smoothing, otherwise use a value of 0 for land points and include them
///                   in the averages.
/// @param bx         The box of this pass, grown by s cells for pass s
/// @param mask       0 for land points and 1 for ocean points at the field's points [nondim]
/// @param f_prev     The value of the field at the previous iteration [arbitrary]
/// @param field      The field to be smoothed [arbitrary]
///
template <bool zero_land>
void smooth_x9_pass(const Box& bx,
                    Array4<const Real> const& mask,
                    Array4<const Real> const& f_prev,
                    Array4<Real>       const& field)
{
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        if (mask(i,j,0) > 0.0_rt) {
            Real Iwts = 0.0625_rt;  ///< The inverse of the sum of the weights [nondim]
            if constexpr (!zero_land) {
                Iwts = 1.0_rt / ( (4.0_rt*mask(i,j,0) +
                                ( 2.0_rt*((mask(i-1,j,0) + mask(i+1,j,0)) +
                                       (mask(i,j-1,0) + mask(i,j+1,0))) +
                                 ((mask(i-1,j-1,0) + mask(i+1,j+1,0)) +
                                  (mask(i-1,j+1,0) + mask(i+1,j-1,0))) ) ) + 1.0e-16_rt );
            }
            field(i,j,0) = Iwts * ( 4.0_rt*mask(i,j,0) * f_prev(i,j,0)
                                  + (2.0_rt*((mask(i-1,j,0) * f_prev(i-1,j,0) + mask(i+1,j,0) * f_prev(i+1,j,0)) +
                                          (mask(i,j-1,0) * f_prev(i,j-1,0) + mask(i,j+1,0) * f_prev(i,j+1,0)))
                                    + ((mask(i-1,j-1,0) * f_prev(i-1,j-1,0) + mask(i+1,j+1,0) * f_prev(i+1,j+1,0)) +
                                       (mask(i-1,j+1,0) * f_prev(i-1,j+1,0) + mask(i+1,j-1,0) * f_prev(i-1,j-1,0))) ));
        }
    });
}

void smooth_x9_h(const Box& bxH0,
                 Array4<const Real> const& mask2dT,
                 Array4<Real>       const& field_h,
                 const bool zero_land)
{
    BL_PROFILE("smooth_x9_h");

    Box bxH0_2d = bxH0;
    bxH0_2d.setSmall(2, 0);
    bxH0_2d.setBig(2, 0);

    /// The value of the h-point field at the previous iteration [arbitrary]
    const Box fh_bx = Box(field_h);
    FArrayBox fh_prev_fab(fh_bx, 1, The_Async_Arena());
    Array4<Real> const fh_prev = fh_prev_fab.array();

    for (int s = 1; s >= 0; --s) {
        /// h points smoothed in pass s, [is-s:ie+s, js-s:je+s]: pass 1 also smooths the halo row that
        /// pass 0 reads
        const Box bxH0pass = grow_by(bxH0_2d, s, s, s, s);

        fh_prev_fab.copy<RunOn::Device>(FArrayBox(field_h));
        /// apply smoothing on field_h using rotationally symmetric expressions.
        if (zero_land) {
            smooth_x9_pass<true>(bxH0pass, mask2dT, fh_prev, field_h);
        } else {
            smooth_x9_pass<false>(bxH0pass, mask2dT, fh_prev, field_h);
        }
    }
}

void smooth_x9_uv(const Box& bxU0,
                  const Box& bxV0,
                  Array4<const Real> const& mask2dCu,
                  Array4<const Real> const& mask2dCv,
                  Array4<Real>       const& field_u,
                  Array4<Real>       const& field_v,
                  const bool zero_land)
{
    BL_PROFILE("smooth_x9_uv");

    Box bxU0_2d = bxU0;
    bxU0_2d.setSmall(2, 0);
    bxU0_2d.setBig(2, 0);
    Box bxV0_2d = bxV0;
    bxV0_2d.setSmall(2, 0);
    bxV0_2d.setBig(2, 0);

    /// The value of the u-point field at the previous iteration [arbitrary]
    const Box fu_bx = Box(field_u);
    FArrayBox fu_prev_fab(fu_bx, 1, The_Async_Arena());
    Array4<Real> const fu_prev = fu_prev_fab.array();
    /// The value of the v-point field at the previous iteration [arbitrary]
    const Box fv_bx = Box(field_v);
    FArrayBox fv_prev_fab(fv_bx, 1, The_Async_Arena());
    Array4<Real> const fv_prev = fv_prev_fab.array();

    for (int s = 1; s >= 0; --s) {
        /// u points smoothed in pass s, [Isq-s:Ieq+s, js-s:je+s]: pass 1 also smooths the halo row that
        /// pass 0 reads
        const Box bxU0pass = grow_by(bxU0_2d, s, s, s, s);
        /// v points smoothed in pass s, [is-s:ie+s, Jsq-s:Jeq+s]: pass 1 also smooths the halo row that
        /// pass 0 reads
        const Box bxV0pass = grow_by(bxV0_2d, s, s, s, s);

        fu_prev_fab.copy<RunOn::Device>(FArrayBox(field_u));
        /// apply smoothing on field_u using the original non-rotationally symmetric expressions.
        if (zero_land) {
            smooth_x9_pass<true>(bxU0pass, mask2dCu, fu_prev, field_u);
        } else {
            smooth_x9_pass<false>(bxU0pass, mask2dCu, fu_prev, field_u);
        }

        fv_prev_fab.copy<RunOn::Device>(FArrayBox(field_v));
        /// apply smoothing on field_v using the original non-rotationally symmetric expressions.
        if (zero_land) {
            smooth_x9_pass<true>(bxV0pass, mask2dCv, fv_prev, field_v);
        } else {
            smooth_x9_pass<false>(bxV0pass, mask2dCv, fv_prev, field_v);
        }
    }
}

///
/// @brief One pass of the 1-1-4-1-1 Laplacian filter applied to a GME diffusive flux.
///
/// @param bx     The range of this smoothing pass
/// @param mask   0 for land points and 1 for ocean points at the flux's points [nondim]
/// @param f_orig The previous value of the GME diffusive flux [L2 T-2 ~> m2 s-2]
/// @param flux   GME diffusive flux [L2 T-2 ~> m2 s-2]
///
static void smooth_gme_pass(const Box& bx,
                            Array4<const Real> const& mask,
                            Array4<const Real> const& f_orig,
                            Array4<Real>       const& flux)
{
    /// apply smoothing on GME
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        /// skip land points
        if (mask(i,j,0) == 0.0_rt) return;
        /// compute weights
        Real const ww = 0.125_rt * mask(i-1,j,0);
        Real const we = 0.125_rt * mask(i+1,j,0);
        Real const ws = 0.125_rt * mask(i,j-1,0);
        Real const wn = 0.125_rt * mask(i,j+1,0);
        Real const wc = 1.0_rt - ((ww+we)+(wn+ws));
        flux(i,j,0) =  wc * f_orig(i,j,0)
                     + ((ww * f_orig(i-1,j,0) + we * f_orig(i+1,j,0))
                      + (ws * f_orig(i,j-1,0) + wn * f_orig(i,j+1,0)));
    });
}

void smooth_gme_h_pass(const Box& bxPass,
                       Array4<const Real> const& mask2dT,
                       Array4<Real>       const& GME_flux_h)
{
    BL_PROFILE("smooth_gme_h_pass");

    Box bxPass_2d = bxPass;
    bxPass_2d.setSmall(2, 0);
    bxPass_2d.setBig(2, 0);

    /// The previous value of GME_flux_h [L2 T-2 ~> m2 s-2]
    const Box orig_bx = Box(GME_flux_h);
    FArrayBox GME_flux_h_original_fab(orig_bx, 1, The_Async_Arena());
    Array4<Real> const GME_flux_h_original = GME_flux_h_original_fab.array();

    GME_flux_h_original_fab.copy<RunOn::Device>(FArrayBox(GME_flux_h));
    smooth_gme_pass(bxPass_2d, mask2dT, GME_flux_h_original, GME_flux_h);
}

void smooth_gme_q_pass(const Box& bxPass,
                       Array4<const Real> const& mask2dBu,
                       Array4<Real>       const& GME_flux_q)
{
    BL_PROFILE("smooth_gme_q_pass");

    Box bxPass_2d = bxPass;
    bxPass_2d.setSmall(2, 0);
    bxPass_2d.setBig(2, 0);

    /// The previous value of GME_flux_q [L2 T-2 ~> m2 s-2]
    const Box orig_bx = Box(GME_flux_q);
    FArrayBox GME_flux_q_original_fab(orig_bx, 1, The_Async_Arena());
    Array4<Real> const GME_flux_q_original = GME_flux_q_original_fab.array();

    GME_flux_q_original_fab.copy<RunOn::Device>(FArrayBox(GME_flux_q));
    smooth_gme_pass(bxPass_2d, mask2dBu, GME_flux_q_original, GME_flux_q);
}

///
/// @brief Adds the MEKE-based (EY24_EBT_BS) backscatter contribution to the diagonal term of
///        the stress tensor at h-points.
///
/// @tparam use_kh_struct    If true, shape the backscatter coefficient with VarMix%BS_struct
/// @tparam do_BS_coeff_h    If true, the BS_coeff_h diagnostic is requested.
/// @param bxHstr            The h points of this tile where the normal stresses are computed,
///                         [IscB:IecB+1, JscB:JecB+1, 1:kmax]
/// @param MEKE_Ku           The MEKE-derived lateral viscosity coefficient [L2 T-1 ~> m2 s-1]
/// @param BS_struct         Vertical structure function used in backscatter [nondim]
/// @param kstart            The first absolute k-layer of the current k-block
/// @param visc_limit_h_flag determines whether backscatter is shut off [nondim]
/// @param sh_xx             horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
/// @param str_xx            The diagonal term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2]
/// @param BS_coeff_h        A diagnostic array of the backscatter coefficient [L2 T-1 ~> m2 s-1]
///
template <bool use_kh_struct, bool do_BS_coeff_h>
void hor_visc_backscatter_h_loop(const Box& bxHstr,
                                 Array4<const Real> const& MEKE_Ku,
                                 Array4<const Real> const& BS_struct,
                                 const int kstart,
                                 Array4<const Real> const& visc_limit_h_flag,
                                 Array4<const Real> const& sh_xx,
                                 Array4<Real>       const& str_xx,
                                 Array4<Real>       const& BS_coeff_h)
{
    ParallelFor(bxHstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
    {
        const int k = kstart + kk - 1;

        Real BS_struct_k = 0.0_rt;
        if constexpr (use_kh_struct) {
            BS_struct_k = BS_struct(i,j,k);
        }

        Real Kh_BS;  ///< Laplacian antiviscosity [L2 T-1 ~> m2 s-1]
        hor_visc_backscatter_h_point<use_kh_struct>(str_xx(i,j,kk), Kh_BS,
                                                    visc_limit_h_flag(i,j,k),
                                                    MEKE_Ku(i,j,0),
                                                    BS_struct_k,
                                                    sh_xx(i,j,kk));

        if constexpr (do_BS_coeff_h) {
            BS_coeff_h(i,j,k) = Kh_BS;
        }
    });
}

void hor_visc_backscatter_h(const Box& bxHstr,
                            const hv_opts_CS_C& CS,
                            Array4<const Real> const& MEKE_Ku,
                            Array4<const Real> const& BS_struct,
                            const bool use_kh_struct,
                            const int kstart,
                            Array4<const Real> const& visc_limit_h_flag,
                            Array4<const Real> const& sh_xx,
                            Array4<Real>       const& str_xx,
                            Array4<Real>       const& BS_coeff_h)
{
    BL_PROFILE("hor_visc_backscatter_h");

    if (use_kh_struct) {
        if (CS.do_BS_coeff_h) {
            hor_visc_backscatter_h_loop<true, true>(bxHstr, MEKE_Ku, BS_struct, kstart,
                                                    visc_limit_h_flag, sh_xx, str_xx, BS_coeff_h);
        } else {
            hor_visc_backscatter_h_loop<true, false>(bxHstr, MEKE_Ku, BS_struct, kstart,
                                                     visc_limit_h_flag, sh_xx, str_xx, BS_coeff_h);
        }
    } else {
        if (CS.do_BS_coeff_h) {
            hor_visc_backscatter_h_loop<false, true>(bxHstr, MEKE_Ku, BS_struct, kstart,
                                                     visc_limit_h_flag, sh_xx, str_xx, BS_coeff_h);
        } else {
            hor_visc_backscatter_h_loop<false, false>(bxHstr, MEKE_Ku, BS_struct, kstart,
                                                      visc_limit_h_flag, sh_xx, str_xx, BS_coeff_h);
        }
    }
}

///
/// @brief Adds the MEKE-based (EY24_EBT_BS) backscatter contribution to the cross term of
///        the stress tensor at q-points.
///
/// @tparam use_kh_struct    If true, shape the backscatter coefficient with VarMix%BS_struct
/// @tparam do_BS_coeff_q    If true, the BS_coeff_q diagnostic is requested.
/// @param bxQstr            The q points of this tile where the shear stresses are computed,
///                          [isc-1:iec, jsc-1:jec, 1:kmax]
/// @param MEKE_Ku           The MEKE-derived lateral viscosity coefficient [L2 T-1 ~> m2 s-1]
/// @param BS_struct         Vertical structure function used in backscatter [nondim]
/// @param kstart            The first absolute k-layer of the current k-block
/// @param visc_limit_q_flag determines whether backscatter is shut off [nondim]
/// @param sh_xy             horizontal shearing strain (du/dy + dv/dx) including metric terms [T-1 ~> s-1]
/// @param str_xy            The cross term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2]
/// @param BS_coeff_q        A diagnostic array of the backscatter coefficient [L2 T-1 ~> m2 s-1]
///
template <bool use_kh_struct, bool do_BS_coeff_q>
void hor_visc_backscatter_q_loop(const Box& bxQstr,
                                 Array4<const Real> const& MEKE_Ku,
                                 Array4<const Real> const& BS_struct,
                                 const int kstart,
                                 Array4<const Real> const& visc_limit_q_flag,
                                 Array4<const Real> const& sh_xy,
                                 Array4<Real>       const& str_xy,
                                 Array4<Real>       const& BS_coeff_q)
{
    ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
    {
        const int k = kstart + kk - 1;

        Real Kh_BS;  ///< Laplacian antiviscosity [L2 T-1 ~> m2 s-1]
        if (visc_limit_q_flag(i,j,k) > 0.0_rt) {
            Kh_BS = 0._rt;
        } else {
            if constexpr (use_kh_struct) {
                Kh_BS = 0.25_rt*( ((MEKE_Ku(i,j,0)*BS_struct(i,j,k)) +
                                (MEKE_Ku(i+1,j+1,0)*BS_struct(i+1,j+1,k))) +
                               ((MEKE_Ku(i+1,j,0)*BS_struct(i+1,j,k)) +
                                (MEKE_Ku(i,j+1,0)*BS_struct(i,j+1,k))) );
            } else {
                Kh_BS = 0.25_rt*( (MEKE_Ku(i,j,0) + MEKE_Ku(i+1,j+1,0)) +
                               (MEKE_Ku(i+1,j,0) + MEKE_Ku(i,j+1,0)) );
            }
        }

        /// The cross term in the stress tensor due to backscatter [H L2 T-2 ~> m3 s-2 or kg s-2]
        Real const str_xy_BS = -Kh_BS * (sh_xy(i,j,kk));

        if constexpr (do_BS_coeff_q) {
            BS_coeff_q(i,j,k) = Kh_BS;
        }

        str_xy(i,j,kk) = str_xy(i,j,kk) + str_xy_BS;
    });
}

void hor_visc_backscatter_q(const Box& bxQstr,
                            const hv_opts_CS_C& CS,
                            Array4<const Real> const& MEKE_Ku,
                            Array4<const Real> const& BS_struct,
                            const bool use_kh_struct,
                            const int kstart,
                            Array4<const Real> const& visc_limit_q_flag,
                            Array4<const Real> const& sh_xy,
                            Array4<Real>       const& str_xy,
                            Array4<Real>       const& BS_coeff_q)
{
    BL_PROFILE("hor_visc_backscatter_q");

    if (use_kh_struct) {
        if (CS.do_BS_coeff_q) {
            hor_visc_backscatter_q_loop<true, true>(bxQstr, MEKE_Ku, BS_struct, kstart,
                                                    visc_limit_q_flag, sh_xy, str_xy, BS_coeff_q);
        } else {
            hor_visc_backscatter_q_loop<true, false>(bxQstr, MEKE_Ku, BS_struct, kstart,
                                                     visc_limit_q_flag, sh_xy, str_xy, BS_coeff_q);
        }
    } else {
        if (CS.do_BS_coeff_q) {
            hor_visc_backscatter_q_loop<false, true>(bxQstr, MEKE_Ku, BS_struct, kstart,
                                                     visc_limit_q_flag, sh_xy, str_xy, BS_coeff_q);
        } else {
            hor_visc_backscatter_q_loop<false, false>(bxQstr, MEKE_Ku, BS_struct, kstart,
                                                      visc_limit_q_flag, sh_xy, str_xy, BS_coeff_q);
        }
    }
}

void hor_visc_gme_bt(const Box& bxH0,
                     const Box& bxQ0,
                     const hv_opts_CS_C& CS,
                     Array4<const Real> const& IdyCu,
                     Array4<const Real> const& IdxCv,
                     Array4<const Real> const& IdyCv,
                     Array4<const Real> const& IdxCu,
                     Array4<const Real> const& mask2dBu,
                     Array4<const Real> const& mask2dCu,
                     Array4<const Real> const& mask2dCv,
                     Array4<const Real> const& mask2dT,
                     Array4<const Real> const& dx_dyT,
                     Array4<const Real> const& dy_dxT,
                     Array4<const Real> const& dx_dyBu,
                     Array4<const Real> const& dy_dxBu,
                     Array4<const Real> const& h,
                     Array4<const Real> const& ubtav,
                     Array4<const Real> const& vbtav,
                     Array4<Real>       const& dudx_bt,
                     Array4<Real>       const& dvdy_bt,
                     Array4<Real>       const& dvdx_bt,
                     Array4<Real>       const& dudy_bt,
                     Array4<Real>       const& sh_xx_bt,
                     Array4<Real>       const& sh_xy_bt,
                     Array4<Real>       const& GME_effic_h,
                     Array4<Real>       const& GME_effic_q)
{
    BL_PROFILE("hor_visc_gme_bt");

    const int nz = bxH0.bigEnd(2) + 1;

    Box bxH0_2d = bxH0;
    bxH0_2d.setSmall(2, 0);
    bxH0_2d.setBig(2, 0);
    Box bxQ0_2d = bxQ0;
    bxQ0_2d.setSmall(2, 0);
    bxQ0_2d.setBig(2, 0);

    /// q points [Isq-2:Ieq+2, Jsq-2:Jeq+2] where dvdx_bt and dudy_bt are computed: GME_effic_h over
    /// bxH0grad averages them from the four surrounding q points
    const Box bxQ0dvdx = grow_by(bxQ0_2d, 2, 2, 2, 2);
    /// h points [Isq-1:Ieq+2, Jsq-1:Jeq+2] where sh_xx_bt and GME_effic_h are computed: the range
    /// over which hor_visc_q_coeffs builds str_xx_GME (bxHgrad); built from bxQ0 because the range
    /// is set in q indices
    const Box bxH0grad = grow_by(bxQ0_2d, 1, 1, 2, 2);
    /// h points [is-2:ie+2, js-2:je+2] where dudx_bt, dvdy_bt and htot are computed: GME_effic_q
    /// over bxQ0shear averages them from the four surrounding h points
    const Box bxH0dudx = grow_by(bxH0_2d, 2, 2, 2, 2);
    /// q points [is-2:ie+1, js-2:je+1] where sh_xy_bt and GME_effic_q are computed: the range over
    /// which hor_visc_q_coeffs builds str_xy_GME (bxQshear); built from bxH0 because the range
    /// starts at is-2 in both memory modes
    const Box bxQ0shear = grow_by(bxH0_2d, 2, 2, 1, 1);

    /// The total thickness of all layers [H ~> m or kg m-2]
    Box htot_bx = Box(h);
    htot_bx.setSmall(2, 0);
    htot_bx.setBig(2, 0);
    FArrayBox htot_fab(htot_bx, 1, The_Async_Arena());
    Array4<Real> const htot = htot_fab.array();

    /// Calculate the barotropic horizontal tension
    ParallelFor(bxH0dudx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        dudx_bt(i,j,0) = dy_dxT(i,j,0)*((IdyCu(i,j,0) * ubtav(i,j,0)) -
                                        (IdyCu(i-1,j,0) * ubtav(i-1,j,0)));
        dvdy_bt(i,j,0) = dx_dyT(i,j,0)*((IdxCv(i,j,0) * vbtav(i,j,0)) -
                                        (IdxCv(i,j-1,0) * vbtav(i,j-1,0)));
    });
    ParallelFor(bxH0grad, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        sh_xx_bt(i,j,0) = dudx_bt(i,j,0) - dvdy_bt(i,j,0);
    });

    /// Components for the barotropic shearing strain
    ParallelFor(bxQ0dvdx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        dvdx_bt(i,j,0) = dy_dxBu(i,j,0)*((vbtav(i+1,j,0)*IdyCv(i+1,j,0))
                                         - (vbtav(i,j,0)*IdyCv(i,j,0)));
        dudy_bt(i,j,0) = dx_dyBu(i,j,0)*((ubtav(i,j+1,0)*IdxCu(i,j+1,0))
                                         - (ubtav(i,j,0)*IdxCu(i,j,0)));
    });

    if (CS.no_slip) {
        ParallelFor(bxQ0shear, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            sh_xy_bt(i,j,0) = (2.0_rt-mask2dBu(i,j,0)) * ( dvdx_bt(i,j,0) + dudy_bt(i,j,0) );
        });
    } else {
        ParallelFor(bxQ0shear, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            sh_xy_bt(i,j,0) = mask2dBu(i,j,0) * ( dvdx_bt(i,j,0) + dudy_bt(i,j,0) );
        });
    }

    ParallelFor(bxH0dudx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        htot(i,j,0) = 0.0_rt;
        for (int k = 0; k < nz; ++k) {
            htot(i,j,0) = htot(i,j,0) + h(i,j,k);
        }
    });

    /// The inverse of GME tapering scale [H-1 ~> m-1 or m2 kg-1]
    const Real I_GME_h0 = 1.0_rt / CS.GME_h0;
    /// The nondimensional prefactor multiplying the GME coefficient [nondim]
    const Real GME_efficiency = CS.GME_efficiency;

    ParallelFor(bxH0grad, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        /// A mask that zeroes out cells with at least one land edge [nondim]
        Real const boundary_mask_h = (mask2dCu(i,j,0) * mask2dCu(i-1,j,0)) * (mask2dCv(i,j,0) * mask2dCv(i,j-1,0));
        Real const dvdx_bt_h = 0.25_rt*((dvdx_bt(i,j,0)+dvdx_bt(i-1,j-1,0)) + (dvdx_bt(i,j-1,0)+dvdx_bt(i-1,j,0)));
        Real const dudy_bt_h = 0.25_rt*((dudy_bt(i,j,0)+dudy_bt(i-1,j-1,0)) + (dudy_bt(i,j-1,0)+dudy_bt(i-1,j,0)));
        /// Magnitude of the barotropic velocity gradient tensor squared at h-points [T-2 ~> s-2]
        Real const grad_vel_mag_bt_h = mask2dT(i,j,0) * boundary_mask_h * (dudx_bt(i,j,0)*dudx_bt(i,j,0) +
              dvdy_bt(i,j,0)*dvdy_bt(i,j,0) + dvdx_bt_h*dvdx_bt_h + dudy_bt_h*dudy_bt_h);
        /// Probably the following test could be simplified to
        /// if (boundary_mask_h * G%mask2dT(I,J) > 0.0) then
        if (grad_vel_mag_bt_h > 0.0_rt) {
            Real const htot_frac = amrex::min(htot(i,j,0) * I_GME_h0, 1.0_rt);
            GME_effic_h(i,j,0) = GME_efficiency * mask2dT(i,j,0) * (htot_frac*htot_frac);
        } else {
            GME_effic_h(i,j,0) = 0.0_rt;
        }
    });

    ParallelFor(bxQ0shear, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        /// A mask that zeroes out cells with at least one land edge [nondim]
        Real const boundary_mask_q = (mask2dCv(i,j,0) * mask2dCv(i+1,j,0)) * (mask2dCu(i,j,0) * mask2dCu(i,j+1,0));
        Real const dudx_bt_q = 0.25_rt*((dudx_bt(i,j,0)+dudx_bt(i+1,j+1,0)) + (dudx_bt(i,j+1,0)+dudx_bt(i+1,j,0)));
        Real const dvdy_bt_q = 0.25_rt*((dvdy_bt(i,j,0)+dvdy_bt(i+1,j+1,0)) + (dvdy_bt(i,j+1,0)+dvdy_bt(i+1,j,0)));
        /// Magnitude of the barotropic velocity gradient tensor squared at q-points [T-2 ~> s-2]
        Real const grad_vel_mag_bt_q = mask2dBu(i,j,0) * boundary_mask_q * (dvdx_bt(i,j,0)*dvdx_bt(i,j,0) +
              dudy_bt(i,j,0)*dudy_bt(i,j,0) + dudx_bt_q*dudx_bt_q + dvdy_bt_q*dvdy_bt_q);
        /// Probably the following test could be simplified to
        /// if (boundary_mask_q * G%mask2dBu(I,J) > 0.0) then
        if (grad_vel_mag_bt_q > 0.0_rt) {
            /// The arithmetic mean total thickness at q points [H ~> m or kg m-2]
            Real const h_arith_q = 0.25_rt * ((htot(i,j,0) + htot(i+1,j+1,0)) + (htot(i+1,j,0) + htot(i,j+1,0)));
            Real const h_arith_frac = amrex::min(h_arith_q * I_GME_h0, 1.0_rt);
            GME_effic_q(i,j,0) = GME_efficiency * mask2dBu(i,j,0) * (h_arith_frac*h_arith_frac);
        } else {
            GME_effic_q(i,j,0) = 0.0_rt;
        }
    });
}

void hor_visc_leith_grad_a(const Box& bxH,
                           const Box& bxKh,
                           const Box& bxKh_i,
                           const Box& bxKh_j,
                           const Box& bxQstr,
                           const Box& bxQvsmooth,
                           const hv_opts_CS_C& CS,
                           Array4<const Real> const& dyBu,
                           Array4<const Real> const& IdxBu,
                           Array4<const Real> const& IdyCu,
                           Array4<const Real> const& dxBu,
                           Array4<const Real> const& IdyBu,
                           Array4<const Real> const& IdxCv,
                           Array4<const Real> const& IdyCv,
                           Array4<const Real> const& IdxCu,
                           Array4<const Real> const& dF_dx,
                           Array4<const Real> const& dF_dy,
                           Array4<const Real> const& dudx,
                           Array4<const Real> const& dvdy,
                           Array4<const Real> const& vort_xy,
                           Array4<const Real> const& vort_xy_smooth,
                           Array4<Real>       const& grad_vort_mag_h_2d,
                           Array4<Real>       const& grad_div_mag_h,
                           Array4<Real>       const& grad_vort_mag_q_2d,
                           Array4<Real>       const& grad_div_mag_q,
                           Array4<Real>       const& Del2vort_q,
                           Array4<Real>       const& vort_xy_dx,
                           Array4<Real>       const& vort_xy_dx_smooth,
                           Array4<Real>       const& div_xx_dy,
                           Array4<Real>       const& vort_xy_dy,
                           Array4<Real>       const& vort_xy_dy_smooth,
                           Array4<Real>       const& div_xx_dx)
{
    BL_PROFILE("hor_visc_leith_grad_a");

    /// h points [is_Kh-1:ie_Kh+1, js_Kh-1:je_Kh+1] where div_xx is computed: the divergence gradients
    /// at u and v points difference it
    const Box bxHdiv = grow_by(bxKh, 1, 1, 1, 1);
    /// v points [is_Kh-1:ie_Kh+1, js-2:je+1] where vort_xy_dx is computed: the gradient magnitudes
    /// over bxKh and over the q points read it on both sides
    const Box bxVvortdx = grow_by(bxKh_i, 1, 2, 1, 1);
    /// u points [is_Kh-1:ie_Kh, js-1:je+1] where div_xx_dx is computed: the gradient magnitudes
    /// over bxKh and bxQstr read it on both sides
    const Box bxUdivdx = grow_by(bxKh_i, 1, 1, 0, 1);
    /// u points [is-2:ie+1, js_Kh-1:je_Kh+1] where vort_xy_dy is computed: the gradient magnitudes
    /// over bxKh and over the q points read it on both sides
    const Box bxUvortdy = grow_by(bxKh_j, 2, 1, 1, 1);
    /// v points [is-1:ie+1, js_Kh-1:je_Kh] where div_xx_dy is computed: the gradient magnitudes
    /// over bxKh and bxQstr read it on both sides
    const Box bxVdivdy = grow_by(bxKh_j, 1, 1, 1, 0);
    /// u points [is_Kh-1:ie_Kh, js_Kh:je_Kh] where vort_xy_dy_smooth is computed:
    /// hor_visc_leith_grad_b reads it on both sides of each bxKh point in i
    const Box bxUvortdysm = grow_by(bxKh, 1, 0, 0, 0);
    /// v points [is_Kh:ie_Kh, js_Kh-1:je_Kh] where vort_xy_dx_smooth is computed:
    /// hor_visc_leith_grad_b reads it on both sides of each bxKh point in j
    const Box bxVvortdxsm = grow_by(bxKh, 0, 1, 0, 0);
    /// u points [is-2:ie+1, js-1:je+1] where beta is added to vort_xy_dy; built from bxH because
    /// this range does not follow the option-dependent halo of bxKh
    const Box bxUbeta = grow_by(bxH, 2, 1, 1, 1);
    /// v points [is-1:ie+1, js-2:je+1] where beta is added to vort_xy_dx; built from bxH because
    /// this range does not follow the option-dependent halo of bxKh
    const Box bxVbeta = grow_by(bxH, 1, 2, 1, 1);

    /// Vorticity gradient
    ParallelFor(bxVvortdx, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
    {
        /// Ratio of meridional over zonal grid spacing at vertices [nondim]
        Real const DY_dxBu = dyBu(i,j,0) * IdxBu(i,j,0);
        vort_xy_dx(i,j,kk) = DY_dxBu * ((vort_xy(i,j,kk) * IdyCu(i,j,0)) - (vort_xy(i-1,j,kk) * IdyCu(i-1,j,0)));
    });

    ParallelFor(bxUvortdy, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
    {
        /// Ratio of zonal over meridional grid spacing at vertices [nondim]
        Real const DX_dyBu = dxBu(i,j,0) * IdyBu(i,j,0);
        vort_xy_dy(i,j,kk) = DX_dyBu * ((vort_xy(i,j,kk) * IdxCv(i,j,0)) - (vort_xy(i,j-1,kk) * IdxCv(i,j-1,0)));
    });

    if (CS.use_Leithy) {
        /// Gradient of smoothed vorticity
        ParallelFor(bxVvortdxsm, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            Real const DY_dxBu = dyBu(i,j,0) * IdxBu(i,j,0);
            vort_xy_dx_smooth(i,j,kk) = DY_dxBu *
                        ((vort_xy_smooth(i,j,kk) * IdyCu(i,j,0)) - (vort_xy_smooth(i-1,j,kk) * IdyCu(i-1,j,0)));
        });

        ParallelFor(bxUvortdysm, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            Real const DX_dyBu = dxBu(i,j,0) * IdyBu(i,j,0);
            vort_xy_dy_smooth(i,j,kk) = DX_dyBu *
                        ((vort_xy_smooth(i,j,kk) * IdxCv(i,j,0)) - (vort_xy_smooth(i,j-1,kk) * IdxCv(i,j-1,0)));
        });
    } /// If Leithy

    /// Laplacian of vorticity
    /// if (CS%Leith_Ah .or. CS%use_Leithy) then
    ParallelFor(bxQvsmooth, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
    {
        Real const DY_dxBu = dyBu(i,j,0) * IdxBu(i,j,0);
        Real const DX_dyBu = dxBu(i,j,0) * IdyBu(i,j,0);

        Del2vort_q(i,j,kk) = DY_dxBu * ((vort_xy_dx(i+1,j,kk) * IdyCv(i+1,j,0)) - (vort_xy_dx(i,j,kk) * IdyCv(i,j,0))) +
                             DX_dyBu * ((vort_xy_dy(i,j+1,kk) * IdyCu(i,j+1,0)) - (vort_xy_dy(i,j,kk) * IdyCu(i,j,0)));
    });
    /// endif

    if (CS.Modified_Leith) {

        /// Estimate of horizontal divergence at h-points [T-1 ~> s-1]
        FArrayBox div_xx_fab(Box(dudx), 1, The_Async_Arena());
        Array4<Real> const div_xx = div_xx_fab.array();

        /// Divergence
        ParallelFor(bxHdiv, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            div_xx(i,j,kk) = dudx(i,j,kk) + dvdy(i,j,kk);
        });

        /// Divergence gradient
        ParallelFor(bxUdivdx, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            div_xx_dx(i,j,kk) = IdxCu(i,j,0)*(div_xx(i+1,j,kk) - div_xx(i,j,kk));
        });
        ParallelFor(bxVdivdy, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            div_xx_dy(i,j,kk) = IdyCv(i,j,0)*(div_xx(i,j+1,kk) - div_xx(i,j,kk));
        });

        /// Magnitude of divergence gradient
        ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            Real const div_xx_dx_h = 0.5_rt*(div_xx_dx(i,j,kk) + div_xx_dx(i-1,j,kk));
            Real const div_xx_dy_h = 0.5_rt*(div_xx_dy(i,j,kk) + div_xx_dy(i,j-1,kk));
            grad_div_mag_h(i,j,kk) = std::sqrt((div_xx_dx_h*div_xx_dx_h) + (div_xx_dy_h*div_xx_dy_h));
        });
        ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            Real const div_xx_dx_q = 0.5_rt*(div_xx_dx(i,j,kk) + div_xx_dx(i,j+1,kk));
            Real const div_xx_dy_q = 0.5_rt*(div_xx_dy(i,j,kk) + div_xx_dy(i+1,j,kk));
            grad_div_mag_q(i,j,kk) = std::sqrt((div_xx_dx_q*div_xx_dx_q) + (div_xx_dy_q*div_xx_dy_q));
        });

    } else {

        FArrayBox(div_xx_dx).setVal<RunOn::Device>(0.0_rt, bxUdivdx);
        FArrayBox(div_xx_dy).setVal<RunOn::Device>(0.0_rt, bxVdivdy);
        FArrayBox(grad_div_mag_h).setVal<RunOn::Device>(0.0_rt, bxKh);
        FArrayBox(grad_div_mag_q).setVal<RunOn::Device>(0.0_rt, bxQstr);

    }

    /// Add in beta for the Leith viscosity
    if (CS.use_beta_in_Leith) {
        ParallelFor(bxVbeta, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            vort_xy_dx(i,j,kk) = vort_xy_dx(i,j,kk) + 0.5_rt * ( dF_dx(i,j,0) + dF_dx(i,j+1,0));
        });
        ParallelFor(bxUbeta, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            vort_xy_dy(i,j,kk) = vort_xy_dy(i,j,kk) + 0.5_rt * ( dF_dy(i,j,0) + dF_dy(i+1,j,0));
        });
    }

    if (CS.use_QG_Leith_visc) {

        ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            Real const vort_xy_dx_h = 0.5_rt*(vort_xy_dx(i,j,kk) + vort_xy_dx(i,j-1,kk));
            Real const vort_xy_dy_h = 0.5_rt*(vort_xy_dy(i,j,kk) + vort_xy_dy(i-1,j,kk));
            grad_vort_mag_h_2d(i,j,kk) = std::sqrt((vort_xy_dx_h*vort_xy_dx_h) + (vort_xy_dy_h*vort_xy_dy_h));
        });
        ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            Real const vort_xy_dx_q = 0.5_rt*(vort_xy_dx(i,j,kk) + vort_xy_dx(i+1,j,kk));
            Real const vort_xy_dy_q = 0.5_rt*(vort_xy_dy(i,j,kk) + vort_xy_dy(i,j+1,kk));
            grad_vort_mag_q_2d(i,j,kk) = std::sqrt((vort_xy_dx_q*vort_xy_dx_q) + (vort_xy_dy_q*vort_xy_dy_q));
        });
    }
}

void hor_visc_leith_grad_b(const Box& bxKh,
                           const Box& bxQstr,
                           const hv_opts_CS_C& CS,
                           Array4<Real>       const& grad_vort_mag_h,
                           Array4<Real>       const& grad_vort_mag_q,
                           Array4<Real>       const& vert_vort_mag_smooth,
                           Array4<const Real> const& vort_xy_dx,
                           Array4<const Real> const& vort_xy_dx_smooth,
                           Array4<const Real> const& vort_xy_dy,
                           Array4<const Real> const& vort_xy_dy_smooth)
{
    BL_PROFILE("hor_visc_leith_grad_b");

    ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
    {
        Real const vort_xy_dx_h = 0.5_rt*(vort_xy_dx(i,j,kk) + vort_xy_dx(i,j-1,kk));
        Real const vort_xy_dy_h = 0.5_rt*(vort_xy_dy(i,j,kk) + vort_xy_dy(i-1,j,kk));
        grad_vort_mag_h(i,j,kk) = std::sqrt((vort_xy_dx_h*vort_xy_dx_h) + (vort_xy_dy_h*vort_xy_dy_h));
    });
    ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
    {
        Real const vort_xy_dx_q = 0.5_rt*(vort_xy_dx(i,j,kk) + vort_xy_dx(i+1,j,kk));
        Real const vort_xy_dy_q = 0.5_rt*(vort_xy_dy(i,j,kk) + vort_xy_dy(i,j+1,kk));
        grad_vort_mag_q(i,j,kk) = std::sqrt((vort_xy_dx_q*vort_xy_dx_q) + (vort_xy_dy_q*vort_xy_dy_q));
    });

    if (CS.use_Leithy) {
        ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            Real const vort_xy_dx_smooth_h = 0.5_rt*(vort_xy_dx_smooth(i,j,kk) + vort_xy_dx_smooth(i,j-1,kk));
            Real const vort_xy_dy_smooth_h = 0.5_rt*(vort_xy_dy_smooth(i,j,kk) + vort_xy_dy_smooth(i-1,j,kk));
            vert_vort_mag_smooth(i,j,kk) = std::sqrt((vort_xy_dx_smooth_h*vort_xy_dx_smooth_h) +
                                                     (vort_xy_dy_smooth_h*vort_xy_dy_smooth_h));
        });
    }
}

void hor_visc_leithy_m(const Box& bxKh,
                       const hv_opts_CS_C& CS,
                       Array4<const Real> const& mask2dBu,
                       Array4<const Real> const& Ah_bg_xx,
                       Array4<const Real> const& m_const_leithy,
                       Array4<const Real> const& m_leithy_max,
                       Array4<const Real> const& Biharm6_const_xx,
                       const Real inv_PI6,
                       Array4<const Real> const& Del2vort_q,
                       Array4<const Real> const& vert_vort_mag,
                       Array4<const Real> const& vort_xy_smooth,
                       Array4<Real>       const& m_leithy)
{
    BL_PROFILE("hor_visc_leithy_m");

    const Real c_K = CS.c_K;

    /// Get m_leithy
    if (CS.smooth_Ah) {
        Box m_leithy_bx = Box(m_leithy);
        m_leithy_bx.setSmall(2, bxKh.smallEnd(2));
        m_leithy_bx.setBig(2, bxKh.bigEnd(2));
        /// This is here to initialize domain edge halo values.
        FArrayBox(m_leithy).setVal<RunOn::Device>(0.0_rt, m_leithy_bx);
    }
    ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
    {
        /// Laplacian of vorticity at h-points [L-2 T-1 ~> m-2 s-1]
        Real const Del2vort_h = 0.25_rt * ((Del2vort_q(i,j,kk) + Del2vort_q(i-1,j-1,kk)) +
                                        (Del2vort_q(i-1,j,kk) + Del2vort_q(i,j-1,kk)));
        /// 2D Leith biharmonic viscosity [L4 T-1 ~> m4 s-1]
        Real const AhLth = Biharm6_const_xx(i,j,0) * inv_PI6 * amrex::Math::abs(Del2vort_h);
        if (AhLth <= Ah_bg_xx(i,j,0)) {
            m_leithy(i,j,kk) = 0.0_rt;
        } else {
            if ((m_const_leithy(i,j,0)*vert_vort_mag(i,j,kk)) < amrex::Math::abs(vort_xy_smooth(i,j,kk))) {
                Real const vort_ratio = vert_vort_mag(i,j,kk) / vort_xy_smooth(i,j,kk);
                m_leithy(i,j,kk) = c_K * (vort_ratio*vort_ratio);
            } else {
                m_leithy(i,j,kk) = m_leithy_max(i,j,0);
            }
            m_leithy(i,j,kk) = mask2dBu(i,j,0) * m_leithy(i,j,kk);
        }
    });
}

void hor_visc_leithy_ah(const Box& bxKh,
                        const hv_opts_CS_C& CS,
                        Array4<const Real> const& Ah_bg_xx,
                        Array4<const Real> const& Biharm6_const_xx,
                        const Real inv_PI6,
                        Array4<const Real> const& Del2vort_q,
                        Array4<const Real> const& vert_vort_mag_smooth,
                        Array4<Real>       const& Ah,
                        Array4<const Real> const& m_leithy,
                        Array4<Real>       const& Ah_sq)
{
    BL_PROFILE("hor_visc_leithy_ah");

    /// Get Ah
    ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
    {
        /// Laplacian of vorticity at h-points [L-2 T-1 ~> m-2 s-1]
        Real const Del2vort_h = 0.25_rt * ((Del2vort_q(i,j,kk) + Del2vort_q(i-1,j-1,kk)) +
                                        (Del2vort_q(i-1,j,kk) + Del2vort_q(i,j-1,kk)));
        /// 2D Leith+E biharmonic viscosity [L4 T-1 ~> m4 s-1]
        Real const AhLthy = Biharm6_const_xx(i,j,0) * inv_PI6 *
                std::sqrt(amrex::max(0._rt, Del2vort_h*Del2vort_h -
                                         m_leithy(i,j,kk)*(vert_vort_mag_smooth(i,j,kk)*vert_vort_mag_smooth(i,j,kk))));
        Ah(i,j,kk) = amrex::max(Ah_bg_xx(i,j,0), AhLthy);
    });
    if (CS.smooth_Ah) {
        /// Smooth Ah before applying upper bound.  Square Ah, then smooth, then take its square root.
        Box Ah_sq_bx = Box(Ah_sq);
        Ah_sq_bx.setSmall(2, bxKh.smallEnd(2));
        Ah_sq_bx.setBig(2, bxKh.bigEnd(2));
        /// This is here to initialize domain edge halo values.
        FArrayBox(Ah_sq).setVal<RunOn::Device>(0.0_rt, Ah_sq_bx);
        ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            Ah_sq(i,j,kk) = Ah(i,j,kk)*Ah(i,j,kk);
        });
    }
}

void hor_visc_leithy_ah_h(const Box& bxKh,
                          const hv_opts_CS_C& CS,
                          Array4<const Real> const& Ah_bg_xx,
                          const int kstart,
                          Array4<Real>       const& Ah,
                          Array4<Real>       const& Ah_h,
                          Array4<const Real> const& Ah_sq)
{
    BL_PROFILE("hor_visc_leithy_ah_h");

    if (CS.smooth_Ah) {
        ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            const int k = kstart + kk - 1;
            hor_visc_leithy_ah_h_smooth_point(Ah(i,j,kk), Ah_h(i,j,k), Ah_bg_xx(i,j,0), Ah_sq(i,j,kk));
        });
    } else {
        ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            const int k = kstart + kk - 1;
            hor_visc_leithy_ah_h_nosmooth_point(Ah_h(i,j,k), Ah(i,j,kk));
        });
    }
}

void hor_visc_strain(const Box& bxH,
                     const Box& bxQ,
                     const Box& bxU,
                     const Box& bxV,
                     const Box& bxKh,
                     const Box& bxHstr,
                     const Box& bxQstr,
                     const Box& bxHgrad,
                     const Box& bxQshear,
                     const Box& bxQvsmooth,
                     const hv_opts_CS_C& CS,
                     Array4<const Real> const& IdyCu,
                     Array4<const Real> const& IdxCv,
                     Array4<const Real> const& IdyCv,
                     Array4<const Real> const& IdxCu,
                     Array4<const Real> const& mask2dT,
                     Array4<const Real> const& dxBu,
                     Array4<const Real> const& mask2dBu,
                     Array4<const Real> const& IareaBu,
                     Array4<const Real> const& dyCv,
                     Array4<const Real> const& dxCu,
                     Array4<const Real> const& dx2h,
                     Array4<const Real> const& dy2h,
                     Array4<const Real> const& dx_dyT,
                     Array4<const Real> const& dy_dxT,
                     Array4<const Real> const& dx2q,
                     Array4<const Real> const& dy2q,
                     Array4<const Real> const& dx_dyBu,
                     Array4<const Real> const& dy_dxBu,
                     Array4<const Real> const& Idx2dyCu,
                     Array4<const Real> const& Idxdy2u,
                     Array4<const Real> const& Idx2dyCv,
                     Array4<const Real> const& Idxdy2v,
                     OceanOBC* OBC,
                     const int kstart,
                     const bool use_vort_xy,
                     const bool apply_OBC,
                     const bool apply_OBC_strain,
                     const bool use_cont_huv,
                     Array4<const Real> const& u,
                     Array4<const Real> const& v,
                     Array4<const Real> const& h,
                     Array4<const Real> const& hu_cont,
                     Array4<const Real> const& hv_cont,
                     Array4<const Real> const& u_smooth,
                     Array4<const Real> const& v_smooth,
                     Array4<Real>       const& NoSt,
                     Array4<Real>       const& ShSt,
                     Array4<Real>       const& Del2u,
                     Array4<Real>       const& h_u,
                     Array4<Real>       const& Del2v,
                     Array4<Real>       const& h_v,
                     Array4<Real>       const& sh_xx,
                     Array4<Real>       const& sh_xx_smooth,
                     Array4<Real>       const& dudx,
                     Array4<Real>       const& dvdy,
                     Array4<Real>       const& dvdx,
                     Array4<Real>       const& dudy,
                     Array4<Real>       const& dvdx_smooth,
                     Array4<Real>       const& dudy_smooth,
                     Array4<Real>       const& sh_xy,
                     Array4<Real>       const& sh_xy_smooth,
                     Array4<Real>       const& vort_xy,
                     Array4<Real>       const& vort_xy_smooth)
{
    BL_PROFILE("hor_visc_strain");

    if (OBC != nullptr) {
        AMREX_ABORT_LOC("OBC pointer provided but not yet implemented");
    }

    /// u points [Isq-1:Ieq+1, js-1:je+1] where Del2u is computed: one point beyond bxU,
    /// because the biharmonic pass in hor_visc_h_coeffs_b differences Del2u between neighbours
    const Box bxUdel2 = grow_by(bxU, 1, 1, 1, 1);
    /// v points [is-1:ie+1, Jsq-1:Jeq+1] where Del2v is computed: one point beyond bxV,
    /// because the biharmonic pass in hor_visc_h_coeffs_b differences Del2v between neighbours
    const Box bxVdel2 = grow_by(bxV, 1, 1, 1, 1);
    /// v points [is-2:ie+2, js-2:Jeq+1] where h_v is set: hrat_min in
    /// hor_visc_h_coeffs_a reads it on both sides of bxKh, and two points out in i
    /// for the OBC cross-corner projections; built from bxH because the range starts
    /// at js-2 in both memory modes
    const Box bxVthick = grow_by(bxH, 2, 2, 2, 1);
    /// u points [is-2:Ieq+1, js-2:je+2] where h_u is set: hrat_min in
    /// hor_visc_h_coeffs_a reads it on both sides of bxKh, and two points out in j
    /// for the OBC cross-corner projections; built from bxH because the range starts
    /// at is-2 in both memory modes
    const Box bxUthick = grow_by(bxH, 2, 2, 1, 2);

    /// q points of this tile where dudy, dvdx and vort_xy are computed,
    /// [is_vort:ie_vort, js_vort:je_vort, 1:kmax]: two points below and one above the h points; with
    /// a Leith scheme, around bxKh instead, because Leith takes gradients of the vorticity
    Box bxVort;
    /// h points where dudx, dvdy and sh_xx are computed: bxHgrad, or with a Leith
    /// scheme [is-2:ie+2, js-2:je+2], because the divergence in hor_visc_leith_grad_a
    /// reads dudx and dvdy one point around bxKh; built from bxH because that range
    /// starts at is-2 in both memory modes
    Box bxHtens;
    /// Set the halo sizes used for the vorticity calculations.
    if ((CS.Leith_Kh) || (CS.Leith_Ah) || (CS.use_Leithy)) {
        bxVort = grow_by(bxKh, 2, 2, 1, 1);
        bxHtens = grow_by(bxH, 2, 2, 2, 2);
    } else {
        bxVort = grow_by(bxH, 2, 2, 1, 1);
        bxHtens = bxHgrad;
    }

    /// The following are the forms of the horizontal tension and horizontal
    /// shearing strain advocated by Smagorinsky (1993) and discussed in
    /// Griffies and Hallberg (2000).

    /// Calculate horizontal tension
    ParallelFor(bxHtens, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
    {
        const int k = kstart + kk - 1;
        dudx(i,j,kk) = dy_dxT(i,j,0)*((IdyCu(i,j,0) * u(i,j,k)) -
                                      (IdyCu(i-1,j,0) * u(i-1,j,k)));
        dvdy(i,j,kk) = dx_dyT(i,j,0)*((IdxCv(i,j,0) * v(i,j,k)) -
                                      (IdxCv(i,j-1,0) * v(i,j-1,k)));
        sh_xx(i,j,kk) = dudx(i,j,kk) - dvdy(i,j,kk);
    });

    /// Components for the shearing strain
    ParallelFor(bxVort, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
    {
        const int k = kstart + kk - 1;
        dvdx(i,j,kk) = dy_dxBu(i,j,0)*((v(i+1,j,k)*IdyCv(i+1,j,0)) - (v(i,j,k)*IdyCv(i,j,0)));
        dudy(i,j,kk) = dx_dyBu(i,j,0)*((u(i,j+1,k)*IdxCu(i,j+1,0)) - (u(i,j,k)*IdxCu(i,j,0)));
    });

    if (CS.use_Leithy) {
        /// Calculate horizontal tension from smoothed velocity
        ParallelFor(bxHstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            const int k = kstart + kk - 1;
            /// components in the horizontal tension from smoothed velocity [T-1 ~> s-1]
            Real const dudx_smooth = dy_dxT(i,j,0)*((IdyCu(i,j,0) * u_smooth(i,j,k)) -
                                                    (IdyCu(i-1,j,0) * u_smooth(i-1,j,k)));
            /// components in the horizontal tension from smoothed velocity [T-1 ~> s-1]
            Real const dvdy_smooth = dx_dyT(i,j,0)*((IdxCv(i,j,0) * v_smooth(i,j,k)) -
                                                    (IdxCv(i,j-1,0) * v_smooth(i,j-1,k)));
            sh_xx_smooth(i,j,kk) = dudx_smooth - dvdy_smooth;
        });

        /// Components for the shearing strain from smoothed velocity
        ParallelFor(bxQvsmooth, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            const int k = kstart + kk - 1;
            dvdx_smooth(i,j,kk) = dy_dxBu(i,j,0) *
                            ((v_smooth(i+1,j,k)*IdyCv(i+1,j,0)) - (v_smooth(i,j,k)*IdyCv(i,j,0)));
            dudy_smooth(i,j,kk) = dx_dyBu(i,j,0) *
                            ((u_smooth(i,j+1,k)*IdxCu(i,j+1,0)) - (u_smooth(i,j,k)*IdxCu(i,j,0)));
        });
    } /// use Leith+E

    if (CS.do_normstress) {
        ParallelFor(bxH, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            const int k = kstart + kk - 1;
            NoSt(i,j,k) = sh_xx(i,j,kk);
        });
    }

    /// Interpolate the thicknesses to velocity points.
    /// The extra wide halos are to accommodate the cross-corner-point projections
    /// in OBCs, which are not ordinarily be necessary, and might not be necessary
    /// even with OBCs if the accelerations are zeroed at OBC points, in which
    /// case the j-loop for h_u could collapse to j=js=1,je+1. -RWH
    if (use_cont_huv) {
        ParallelFor(bxUthick, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            const int k = kstart + kk - 1;
            h_u(i,j,kk) = hu_cont(i,j,k);
        });
        ParallelFor(bxVthick, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            const int k = kstart + kk - 1;
            h_v(i,j,kk) = hv_cont(i,j,k);
        });
    } else if (CS.use_land_mask) {
        ParallelFor(bxUthick, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            const int k = kstart + kk - 1;
            h_u(i,j,kk) = 0.5_rt * (mask2dT(i,j,0)*h(i,j,k) + mask2dT(i+1,j,0)*h(i+1,j,k));
        });
        ParallelFor(bxVthick, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            const int k = kstart + kk - 1;
            h_v(i,j,kk) = 0.5_rt * (mask2dT(i,j,0)*h(i,j,k) + mask2dT(i,j+1,0)*h(i,j+1,k));
        });
    } else {
        ParallelFor(bxUthick, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            const int k = kstart + kk - 1;
            h_u(i,j,kk) = 0.5_rt * (h(i,j,k) + h(i+1,j,k));
        });
        ParallelFor(bxVthick, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            const int k = kstart + kk - 1;
            h_v(i,j,kk) = 0.5_rt * (h(i,j,k) + h(i,j+1,k));
        });
    }

    /// Adjust contributions to shearing strain and interpolated values of
    /// thicknesses on open boundaries.
    /*
    if (apply_OBC) {
        const int is  = bxH.smallEnd(0);
        const int ie  = bxH.bigEnd(0);
        const int js  = bxH.smallEnd(1);
        const int je  = bxH.bigEnd(1);
        const int Ieq = bxQ.bigEnd(0);
        const int Jeq = bxQ.bigEnd(1);
        const int klo = bxH.smallEnd(2);
        const int khi = bxH.bigEnd(2);
        /// Loop ranges for vorticity terms
        const int is_vort = bxVort.smallEnd(0);
        const int ie_vort = bxVort.bigEnd(0);
        const int js_vort = bxVort.smallEnd(1);
        const int je_vort = bxVort.bigEnd(1);
        const bool use_Leithy = CS.use_Leithy;
        const int strain_config = OBC->strain_config;
        for (int n = 0; n < OBC->number_of_segments; ++n) {
            auto const& seg = OBC->segment[n];
            const int J = seg.HI.JsdB;
            const int I = seg.HI.IsdB;
            const int direction = seg.direction;
            auto const& tangential_vel  = seg.tangential_vel;
            auto const& tangential_grad = seg.tangential_grad;
            if (apply_OBC_strain) {
                if (seg.is_N_or_S && (J >= js_vort) && (J <= je_vort)) {
                    const Box seg_bx(IntVect(amrex::max(seg.HI.IsdB, is_vort), J, klo),
                                     IntVect(amrex::min(seg.HI.IedB, ie_vort), J, khi));
                    ParallelFor(seg_bx, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                    {
                        const int k = kstart + kk - 1;
                        switch (strain_config) {
                            case OBC_STRAIN_ZERO:
                                dvdx(i,j,kk) = 0._rt; dudy(i,j,kk) = 0._rt;
                                break;
                            case OBC_STRAIN_FREESLIP:
                                dudy(i,j,kk) = 0._rt;
                                break;
                            case OBC_STRAIN_COMPUTED:
                                if (direction == OBC_DIRECTION_N) {
                                    dudy(i,j,kk) = 2.0_rt*dx_dyBu(i,j,0)*
                                                (tangential_vel(i,j,k) - u(i,j,k))*IdxCu(i,j,0);
                                } else {
                                    dudy(i,j,kk) = 2.0_rt*dx_dyBu(i,j,0)*
                                                (u(i,j+1,k) - tangential_vel(i,j,k))*IdxCu(i,j+1,0);
                                }
                                break;
                            case OBC_STRAIN_SPECIFIED:
                                if (direction == OBC_DIRECTION_N) {
                                    dudy(i,j,kk) = dx_dyBu(i,j,0)*tangential_grad(i,j,k)*IdxCu(i,j,0)*dxBu(i,j,0);
                                } else {
                                    dudy(i,j,kk) = dx_dyBu(i,j,0)*tangential_grad(i,j,k)*IdxCu(i,j+1,0)*dxBu(i,j,0);
                                }
                                break;
                        }
                        if (use_Leithy) {
                            dvdx_smooth(i,j,kk) = dvdx(i,j,kk);
                            dudy_smooth(i,j,kk) = dudy(i,j,kk);
                        }
                    });
                } else if (seg.is_E_or_W && (I >= is_vort) && (I <= ie_vort)) {
                    const Box seg_bx(IntVect(I, amrex::max(seg.HI.JsdB, js_vort), klo),
                                     IntVect(I, amrex::min(seg.HI.JedB, je_vort), khi));
                    ParallelFor(seg_bx, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                    {
                        const int k = kstart + kk - 1;
                        switch (strain_config) {
                            case OBC_STRAIN_ZERO:
                                dvdx(i,j,kk) = 0._rt; dudy(i,j,kk) = 0._rt;
                                break;
                            case OBC_STRAIN_FREESLIP:
                                dvdx(i,j,kk) = 0._rt;
                                break;
                            case OBC_STRAIN_COMPUTED:
                                if (direction == OBC_DIRECTION_E) {
                                    dvdx(i,j,kk) = 2.0_rt*dy_dxBu(i,j,0)*
                                                (tangential_vel(i,j,k) - v(i,j,k))*IdyCv(i,j,0);
                                } else {
                                    dvdx(i,j,kk) = 2.0_rt*dy_dxBu(i,j,0)*
                                                (v(i+1,j,k) - tangential_vel(i,j,k))*IdyCv(i+1,j,0);
                                }
                                break;
                            case OBC_STRAIN_SPECIFIED:
                                if (direction == OBC_DIRECTION_E) {
                                    dvdx(i,j,kk) = dy_dxBu(i,j,0)*tangential_grad(i,j,k)*IdyCv(i,j,0)*dxBu(i,j,0);
                                } else {
                                    dvdx(i,j,kk) = dy_dxBu(i,j,0)*tangential_grad(i,j,k)*IdyCv(i+1,j,0)*dxBu(i,j,0);
                                }
                                break;
                        }
                        if (use_Leithy) {
                            dvdx_smooth(i,j,kk) = dvdx(i,j,kk);
                            dudy_smooth(i,j,kk) = dudy(i,j,kk);
                        }
                    });
                }
            }

            if (direction == OBC_DIRECTION_N) {
                /// There are extra wide halos here to accommodate the cross-corner-point
                /// OBC projections, but they might not be necessary if the accelerations
                /// are always zeroed out at OBC points, in which case the i-loop below
                /// becomes do i=is-1,ie+1. -RWH
                if ((J >= js-2) && (J <= Jeq+1)) {
                    const Box seg_bx(IntVect(amrex::max(is-2, seg.HI.isd), J, klo),
                                     IntVect(amrex::min(ie+2, seg.HI.ied), J, khi));
                    ParallelFor(seg_bx, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                    {
                        const int k = kstart + kk - 1;
                        h_v(i,j,kk) = h(i,j,k);
                    });
                }
            } else if (direction == OBC_DIRECTION_S) {
                if ((J >= js-2) && (J <= Jeq+1)) {
                    const Box seg_bx(IntVect(amrex::max(is-2, seg.HI.isd), J, klo),
                                     IntVect(amrex::min(ie+2, seg.HI.ied), J, khi));
                    ParallelFor(seg_bx, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                    {
                        const int k = kstart + kk - 1;
                        h_v(i,j,kk) = h(i,j+1,k);
                    });
                }
            } else if (direction == OBC_DIRECTION_E) {
                if ((I >= is-2) && (I <= Ieq+1)) {
                    const Box seg_bx(IntVect(I, amrex::max(js-2, seg.HI.jsd), klo),
                                     IntVect(I, amrex::min(je+2, seg.HI.jed), khi));
                    ParallelFor(seg_bx, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                    {
                        const int k = kstart + kk - 1;
                        h_u(i,j,kk) = h(i,j,k);
                    });
                }
            } else if (direction == OBC_DIRECTION_W) {
                if ((I >= is-2) && (I <= Ieq+1)) {
                    const Box seg_bx(IntVect(I, amrex::max(js-2, seg.HI.jsd), klo),
                                     IntVect(I, amrex::min(je+2, seg.HI.jed), khi));
                    ParallelFor(seg_bx, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                    {
                        const int k = kstart + kk - 1;
                        h_u(i,j,kk) = h(i+1,j,k);
                    });
                }
            }
        }
        /// Now project thicknesses across corner points on OBCs.
        for (int n = 0; n < OBC->number_of_segments; ++n) {
            auto const& seg = OBC->segment[n];
            const int J = seg.HI.JsdB;
            const int I = seg.HI.IsdB;
            if (seg.direction == OBC_DIRECTION_N) {
                if ((J >= js-2) && (J <= je)) {
                    const Box seg_bx(IntVect(amrex::max(is-2, seg.HI.IsdB), J, klo),
                                     IntVect(amrex::min(Ieq+1, seg.HI.IedB), J, khi));
                    ParallelFor(seg_bx, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                    {
                        h_u(i,j+1,kk) = h_u(i,j,kk);
                    });
                }
            } else if (seg.direction == OBC_DIRECTION_S) {
                if ((J >= js-1) && (J <= je+1)) {
                    const Box seg_bx(IntVect(amrex::max(is-2, seg.HI.isd), J, klo),
                                     IntVect(amrex::min(Ieq+1, seg.HI.ied), J, khi));
                    ParallelFor(seg_bx, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                    {
                        h_u(i,j,kk) = h_u(i,j+1,kk);
                    });
                }
            } else if (seg.direction == OBC_DIRECTION_E) {
                if ((I >= is-2) && (I <= ie)) {
                    const Box seg_bx(IntVect(I, amrex::max(js-2, seg.HI.jsd), klo),
                                     IntVect(I, amrex::min(Jeq+1, seg.HI.jed), khi));
                    ParallelFor(seg_bx, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                    {
                        h_v(i+1,j,kk) = h_v(i,j,kk);
                    });
                }
            } else if (seg.direction == OBC_DIRECTION_W) {
                if ((I >= is-1) && (I <= ie+1)) {
                    const Box seg_bx(IntVect(I, amrex::max(js-2, seg.HI.jsd), klo),
                                     IntVect(I, amrex::min(Jeq+1, seg.HI.jed), khi));
                    ParallelFor(seg_bx, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                    {
                        h_v(i,j,kk) = h_v(i+1,j,kk);
                    });
                }
            }
        }
    }
    */

    /// Shearing strain (including no-slip boundary conditions at the 2-D land-sea mask).
    /// dudy and dvdx include modifications at OBCs from above.
    if (CS.no_slip) {
        ParallelFor(bxQshear, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            sh_xy(i,j,kk) = (2.0_rt-mask2dBu(i,j,0)) * ( dvdx(i,j,kk) + dudy(i,j,kk) );
        });
    } else {
        ParallelFor(bxQshear, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            sh_xy(i,j,kk) = mask2dBu(i,j,0) * ( dvdx(i,j,kk) + dudy(i,j,kk) );
        });
    }
    if (CS.do_shearstress) {
        ParallelFor(bxQshear, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            const int k = kstart + kk - 1;
            ShSt(i,j,k) = sh_xy(i,j,kk);
        });
    }

    if (CS.use_Leithy) {
        /// Shearing strain (including no-slip boundary conditions at the 2-D land-sea mask).
        /// dudy_smooth and dvdx_smooth do not (yet) include modifications at OBCs from above.
        if (CS.no_slip) {
            ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                sh_xy_smooth(i,j,kk) = (2.0_rt-mask2dBu(i,j,0)) * ( dvdx_smooth(i,j,kk) + dudy_smooth(i,j,kk) );
            });
        } else {
            ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                sh_xy_smooth(i,j,kk) = mask2dBu(i,j,0) * ( dvdx_smooth(i,j,kk) + dudy_smooth(i,j,kk) );
            });
        }
    } /// use Leith+E

    ///  Evaluate Del2u = x.Div(Grad u) and Del2v = y.Div( Grad u)
    if (CS.biharmonic) {
        ParallelFor(bxUdel2, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            Del2u(i,j,kk) = Idx2dyCu(i,j,0) * ((dx2q(i,j,0)*sh_xy(i,j,kk)) - (dx2q(i,j-1,0)*sh_xy(i,j-1,kk))) +
                            Idxdy2u(i,j,0) * ((dy2h(i+1,j,0)*sh_xx(i+1,j,kk)) - (dy2h(i,j,0)*sh_xx(i,j,kk)));
        });
        ParallelFor(bxVdel2, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            Del2v(i,j,kk) = Idxdy2v(i,j,0) * ((dy2q(i,j,0)*sh_xy(i,j,kk)) - (dy2q(i-1,j,0)*sh_xy(i-1,j,kk))) -
                            Idx2dyCv(i,j,0) * ((dx2h(i,j+1,0)*sh_xx(i,j+1,kk)) - (dx2h(i,j,0)*sh_xx(i,j,kk)));
        });
        /*
        if (apply_OBC) {
            if (OBC->zero_biharmonic) {
                const int Isq = bxQ.smallEnd(0);
                const int Ieq = bxQ.bigEnd(0);
                const int Jsq = bxQ.smallEnd(1);
                const int Jeq = bxQ.bigEnd(1);
                const int klo = bxH.smallEnd(2);
                const int khi = bxH.bigEnd(2);
                for (int n = 0; n < OBC->number_of_segments; ++n) {
                    auto const& seg = OBC->segment[n];
                    const int I = seg.HI.IsdB;
                    const int J = seg.HI.JsdB;
                    if (seg.is_N_or_S && (J >= Jsq-1) && (J <= Jeq+1)) {
                        const Box seg_bx(IntVect(seg.HI.isd, J, klo), IntVect(seg.HI.ied, J, khi));
                        FArrayBox(Del2v).setVal<RunOn::Device>(0._rt, seg_bx);
                    } else if (seg.is_E_or_W && (I >= Isq-1) && (I <= Ieq+1)) {
                        const Box seg_bx(IntVect(I, seg.HI.jsd, klo), IntVect(I, seg.HI.jed, khi));
                        FArrayBox(Del2u).setVal<RunOn::Device>(0._rt, seg_bx);
                    }
                }
            }
        }
        */
    }

    /// Vorticity
    if (use_vort_xy) {
        if (CS.no_slip) {
            ParallelFor(bxVort, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                vort_xy(i,j,kk) = (2.0_rt-mask2dBu(i,j,0)) * ( dvdx(i,j,kk) - dudy(i,j,kk) );
            });
        } else {
            if (CS.use_circulation) {
                ParallelFor(bxVort, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                {
                    const int k = kstart + kk - 1;
                    vort_xy(i,j,kk) = mask2dBu(i,j,0) * IareaBu(i,j,0) * (
                        ((v(i+1,j,k)*dyCv(i+1,j,0)) - (v(i,j,k)*dyCv(i,j,0)))
                      - ((u(i,j+1,k)*dxCu(i,j+1,0)) - (u(i,j,k)*dxCu(i,j,0)))
                    );
                });
            } else {
                ParallelFor(bxVort, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                {
                    vort_xy(i,j,kk) = mask2dBu(i,j,0) * ( dvdx(i,j,kk) - dudy(i,j,kk) );
                });
            }
        }
    }

    if (CS.use_Leithy) {
        if (CS.no_slip) {
            ParallelFor(bxQvsmooth, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                vort_xy_smooth(i,j,kk) = (2.0_rt-mask2dBu(i,j,0)) * ( dvdx_smooth(i,j,kk) - dudy_smooth(i,j,kk) );
            });
        } else {
            ParallelFor(bxQvsmooth, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                vort_xy_smooth(i,j,kk) = mask2dBu(i,j,0) * ( dvdx_smooth(i,j,kk) - dudy_smooth(i,j,kk) );
            });
        }
    }
}

void hor_visc_h_coeffs_a(const Box& bxH,
                         const Box& bxKh,
                         const Box& bxHstr,
                         const hv_opts_CS_C& CS,
                         Array4<const Real> const& MEKE_Ku,
                         Array4<const Real> const& Res_fn_h,
                         Array4<const Real> const& BS_struct,
                         Array4<const Real> const& Kh_bg_xx,
                         Array4<const Real> const& Ah_bg_xx,
                         Array4<const Real> const& Kh_Max_xx,
                         Array4<const Real> const& n1n2_h,
                         Array4<const Real> const& n1n1_m_n2n2_h,
                         Array4<const Real> const& grid_sp_h2,
                         Array4<const Real> const& Laplac2_const_xx,
                         Array4<const Real> const& Biharm6_const_xx,
                         Array4<const Real> const& Laplac3_const_xx,
                         Array4<const Real> const& Biharm_const_xx,
                         Array4<const Real> const& Biharm_const2_xx,
                         const int kstart,
                         const Real h_neglect,
                         const Real inv_PI3,
                         const Real inv_PI6,
                         const bool use_kh_struct,
                         const bool use_MEKE_Ku,
                         const bool rescale_Kh,
                         Array4<const Real> const& u,
                         Array4<const Real> const& v,
                         Array4<const Real> const& h,
                         Array4<Real>       const& Kh_h,
                         Array4<Real>       const& grid_Re_Kh,
                         Array4<Real>       const& div_xx_h,
                         Array4<Real>       const& sh_xx_h,
                         Array4<const Real> const& h_u,
                         Array4<const Real> const& h_v,
                         Array4<const Real> const& sh_xx,
                         Array4<Real>       const& str_xx,
                         Array4<const Real> const& grad_vort_mag_h,
                         Array4<const Real> const& grad_vort_mag_h_2d,
                         Array4<const Real> const& grad_div_mag_h,
                         Array4<const Real> const& dudx,
                         Array4<const Real> const& dvdy,
                         Array4<const Real> const& sh_xy,
                         Array4<const Real> const& Del2vort_q,
                         Array4<Real>       const& Ah,
                         Array4<Real>       const& Kh,
                         Array4<Real>       const& Shear_mag,
                         Array4<Real>       const& vert_vort_mag,
                         Array4<Real>       const& hrat_min,
                         Array4<Real>       const& visc_bound_rem)
{
    BL_PROFILE("hor_visc_h_coeffs_a");

    if ((CS.Smagorinsky_Kh) || (CS.Smagorinsky_Ah)) {
        ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            /// Square of tension (sh_xx) [T-2 ~> s-2]
            Real const sh_xx_sq = sh_xx(i,j,kk)*sh_xx(i,j,kk);
            /// Square of shearing strain (sh_xy) [T-2 ~> s-2]
            Real const sh_xy_sq = 0.25_rt * ( ((sh_xy(i-1,j-1,kk)*sh_xy(i-1,j-1,kk)) + (sh_xy(i,j,kk)*sh_xy(i,j,kk)))
                                         + ((sh_xy(i-1,j,kk)*sh_xy(i-1,j,kk)) + (sh_xy(i,j-1,kk)*sh_xy(i,j-1,kk))) );
            Shear_mag(i,j,kk) = std::sqrt(sh_xx_sq + sh_xy_sq);
        });
    }

    if (CS.bound_Ah || CS.bound_Kh) {
        ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            const int k = kstart + kk - 1;
            /// Minimum h at the 4 neighboring velocity points [H ~> m]
            Real const h_min = amrex::min(h_u(i,j,kk), h_u(i-1,j,kk), h_v(i,j,kk), h_v(i,j-1,kk));
            hrat_min(i,j,kk) = amrex::min(1.0_rt, h_min / (h(i,j,k) + h_neglect));
        });
    }

    if (CS.Laplacian) {
        /// Determine the Laplacian viscosity at h points, using the
        /// largest value from several parameterizations. Also get
        /// the Laplacian component of str_xx.

        if ((CS.Leith_Kh) || (CS.Leith_Ah) || (CS.use_Leithy)) {
            if (CS.use_QG_Leith_visc) {
                ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                {
                    /// Vorticity gradient magnitude [L-1 T-1 ~> m-1 s-1]
                    Real const grad_vort = grad_vort_mag_h(i,j,kk) + grad_div_mag_h(i,j,kk);
                    /// QG-based vorticity gradient magnitude [L-1 T-1 ~> m-1 s-1]
                    Real const grad_vort_qg = 3._rt * grad_vort_mag_h_2d(i,j,kk);
                    vert_vort_mag(i,j,kk) = amrex::min(grad_vort, grad_vort_qg);
                });
            } else {
                ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                {
                    vert_vort_mag(i,j,kk) = grad_vort_mag_h(i,j,kk) + grad_div_mag_h(i,j,kk);
                });
            }
        }

        /// Static (pre-computed) background viscosity
        ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            Kh(i,j,kk) = Kh_bg_xx(i,j,0);
        });

        /// NOTE: The following do-block can be decomposed and vectorized after the
        ///   stack size has been reduced.
        if (CS.add_LES_viscosity) {
            if (CS.Smagorinsky_Kh) {
                ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                {
                    Kh(i,j,kk) = Kh(i,j,kk) + Laplac2_const_xx(i,j,0) * Shear_mag(i,j,kk);
                });
            }
            if (CS.Leith_Kh) {
                ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                {
                    Kh(i,j,kk) = Kh(i,j,kk) + Laplac3_const_xx(i,j,0) * vert_vort_mag(i,j,kk) * inv_PI3;
                });
            }
        } else {
            if (CS.Smagorinsky_Kh) {
                ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                {
                    Kh(i,j,kk) = amrex::max(Kh(i,j,kk), Laplac2_const_xx(i,j,0) * Shear_mag(i,j,kk));
                });
            }
            if (CS.Leith_Kh) {
                ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                {
                    Kh(i,j,kk) = amrex::max(Kh(i,j,kk), Laplac3_const_xx(i,j,0) * vert_vort_mag(i,j,kk) * inv_PI3);
                });
            }
        }

        /// All viscosity contributions above are subject to resolution scaling

        if (rescale_Kh) {
            ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                Kh(i,j,kk) = Res_fn_h(i,j,0) * Kh(i,j,kk);
            });
        }

        /// Place a floor on the viscosity, if desired.
        const Real Kh_bg_min = CS.Kh_bg_min;
        ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            Kh(i,j,kk) = amrex::max(Kh(i,j,kk), Kh_bg_min);
        });

        if (use_MEKE_Ku && !CS.EY24_EBT_BS) {
            /// *Add* the MEKE contribution (which might be negative)
            if (use_kh_struct) {
                if (CS.res_scale_MEKE) {
                    ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                    {
                        const int k = kstart + kk - 1;
                        Kh(i,j,kk) = Kh(i,j,kk) + MEKE_Ku(i,j,0) * Res_fn_h(i,j,0) * BS_struct(i,j,k);
                    });
                } else {
                    ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                    {
                        const int k = kstart + kk - 1;
                        Kh(i,j,kk) = Kh(i,j,kk) + MEKE_Ku(i,j,0) * BS_struct(i,j,k);
                    });
                }
            } else {
                if (CS.res_scale_MEKE) {
                    ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                    {
                        Kh(i,j,kk) = Kh(i,j,kk) + MEKE_Ku(i,j,0) * Res_fn_h(i,j,0);
                    });
                } else {
                    ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                    {
                        Kh(i,j,kk) = Kh(i,j,kk) + MEKE_Ku(i,j,0);
                    });
                }
            }
        }

        if (CS.anisotropic) {
            const Real Kh_aniso = CS.Kh_aniso;
            ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                /// *Add* the tension component of anisotropic viscosity
                Kh(i,j,kk) = Kh(i,j,kk) + Kh_aniso * (1._rt - n1n2_h(i,j,0)*n1n2_h(i,j,0));
            });
        }


        /// Newer method of bounding for stability
        if ((CS.bound_Kh) && (CS.bound_Ah)) {
            const bool backscatter_underbound = CS.backscatter_underbound;
            ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                visc_bound_rem(i,j,kk) = 1.0_rt;
                /// The local maximum Laplacian viscosity for stability [L2 T-1 ~> m2 s-1]
                Real const Kh_max_here = hrat_min(i,j,kk) * Kh_Max_xx(i,j,0);
                if (Kh(i,j,kk) >= Kh_max_here) {
                    visc_bound_rem(i,j,kk) = 0.0_rt;
                    Kh(i,j,kk) = Kh_max_here;
                } else if ((Kh(i,j,kk) > 0.0_rt) || (backscatter_underbound && (Kh_max_here > 0.0_rt))) {
                    visc_bound_rem(i,j,kk) = 1.0_rt - Kh(i,j,kk) / Kh_max_here;
                }
            });
        } else if (CS.bound_Kh) {
            ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                Kh(i,j,kk) = amrex::min(Kh(i,j,kk), hrat_min(i,j,kk) * Kh_Max_xx(i,j,0));
            });
        }

        /// In Leith+E parameterization Kh is computed after Ah in the biharmonic loop.
        /// The harmonic component of str_xx is added in the biharmonic loop.
        if (CS.use_Leithy) {
            FArrayBox(Kh).setVal<RunOn::Device>(0._rt, bxKh);
        }

        if (CS.do_Kh_h || CS.debug) {
            ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                Kh_h(i,j,k) = Kh(i,j,kk);
            });
        }

        if (CS.do_grid_Re_Kh) {
            const Real min_grid_Kh = CS.min_grid_Kh;
            ParallelFor(bxH, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                Real const u_sum = u(i,j,k)+u(i-1,j,k);
                Real const v_sum = v(i,j,k)+v(i,j-1,k);
                /// Local kinetic energy [L2 T-2 ~> m2 s-2]
                Real const KE = 0.125_rt*((u_sum*u_sum) + (v_sum*v_sum));
                /// Laplacian viscosity bound by grid [L2 T-1 ~> m2 s-1]
                Real const grid_Kh = amrex::max(Kh(i,j,kk), min_grid_Kh);
                grid_Re_Kh(i,j,k) = (std::sqrt(KE) * std::sqrt(grid_sp_h2(i,j,0))) / grid_Kh;
            });
        }

        if (CS.do_div_xx_h) {
            ParallelFor(bxH, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                div_xx_h(i,j,k) = dudx(i,j,kk) + dvdy(i,j,kk);
            });
        }

        if (CS.do_sh_xx_h) {
            ParallelFor(bxH, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                sh_xx_h(i,j,k) = sh_xx(i,j,kk);
            });
        }

        ParallelFor(bxHstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            str_xx(i,j,kk) = -Kh(i,j,kk) * sh_xx(i,j,kk);
        });
    } else {
        FArrayBox(str_xx).setVal<RunOn::Device>(0.0_rt, bxHstr);
    } /// Get Kh at h points and get Laplacian component of str_xx

    if (CS.anisotropic) {
        const Real Kh_aniso = CS.Kh_aniso;
        ParallelFor(bxHstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            /// Shearing-strain averaged to h-points
            /// Local variable for interpolating computed strain rates [T-1 ~> s-1].
            Real const local_strain = 0.25_rt * ( (sh_xy(i,j,kk) + sh_xy(i-1,j-1,kk)) +
                                               (sh_xy(i-1,j,kk) + sh_xy(i,j-1,kk)) );
            /// *Add* the shear-strain contribution to the xx-component of stress
            str_xx(i,j,kk) = str_xx(i,j,kk) - Kh_aniso * n1n2_h(i,j,0) * n1n1_m_n2n2_h(i,j,0) * local_strain;
        });
    }

    if (CS.biharmonic) {
        /// Determine the biharmonic viscosity at h points, using the
        /// largest value from several parameterizations. Also get the
        /// biharmonic component of str_xx.
        ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            Ah(i,j,kk) = Ah_bg_xx(i,j,0);
        });

        if ((CS.Smagorinsky_Ah) || (CS.Leith_Ah) || (CS.use_Leithy)) {
            if (CS.Smagorinsky_Ah) {
                if (CS.bound_Coriolis) {
                    ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                    {
                        /// Smagorinsky biharmonic viscosity [L4 T-1 ~> m4 s-1]
                        Real const AhSm = Shear_mag(i,j,kk) * (Biharm_const_xx(i,j,0)
                            + Biharm_const2_xx(i,j,0) * Shear_mag(i,j,kk));
                        Ah(i,j,kk) = amrex::max(Ah(i,j,kk), AhSm);
                    });
                } else {
                    ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                    {
                        /// Smagorinsky biharmonic viscosity [L4 T-1 ~> m4 s-1]
                        Real const AhSm = Biharm_const_xx(i,j,0) * Shear_mag(i,j,kk);
                        Ah(i,j,kk) = amrex::max(Ah(i,j,kk), AhSm);
                    });
                }
            }

            if (CS.Leith_Ah) {
                ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                {
                    /// Laplacian of vorticity at h-points [L-2 T-1 ~> m-2 s-1]
                    Real const Del2vort_h = 0.25_rt * ((Del2vort_q(i,j,kk) + Del2vort_q(i-1,j-1,kk)) +
                                                    (Del2vort_q(i-1,j,kk) + Del2vort_q(i,j-1,kk)));
                    /// 2D Leith biharmonic viscosity [L4 T-1 ~> m4 s-1]
                    Real const AhLth = Biharm6_const_xx(i,j,0) * amrex::Math::abs(Del2vort_h) * inv_PI6;
                    Ah(i,j,kk) = amrex::max(Ah(i,j,kk), AhLth);
                });
            }
        } /// Smagorinsky_Ah or Leith_Ah or Leith+E
    }
}

///
/// @brief Diagnoses str_xx*d_x u - str_yy*d_y v + str_xy*(d_y u + d_x v) at one h-point, using the old
///        formulation that includes energy diffusion.
///
/// @param i, j, k, kk The point, with k the absolute layer and kk its index in the k-block
/// @param H_to_RZ     A constant that translates thickness units to the units of mass per unit area
///                    [R Z H-1 ~> kg m-3 or 1]
/// @param str_xx      The diagonal term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2]
/// @param str_xy      The cross term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2]
/// @param u           The zonal velocity [L T-1 ~> m s-1]
/// @param v           The meridional velocity [L T-1 ~> m s-1]
/// @param IdxT        1/dxT [L-1 ~> m-1]
/// @param IdyT        IdyT is 1/dyT [L-1 ~> m-1]
/// @param IdyBu       1/dyBu [L-1 ~> m-1]
/// @param IdxBu       1/dxBu [L-1 ~> m-1]
/// @return            The work done by the stresses [R Z L2 T-3 ~> W m-2]
///
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
static Real hor_visc_frict_work_bug(const int i, const int j, const int k, const int kk,
                                    Real const H_to_RZ,
                                    Array4<const Real> const& str_xx,
                                    Array4<const Real> const& str_xy,
                                    Array4<const Real> const& u,
                                    Array4<const Real> const& v,
                                    Array4<const Real> const& IdxT,
                                    Array4<const Real> const& IdyT,
                                    Array4<const Real> const& IdyBu,
                                    Array4<const Real> const& IdxBu) noexcept
{
    return H_to_RZ * (
        ((str_xx(i,j,kk) * (u(i,j,k)-u(i-1,j,k))*IdxT(i,j,0))
        - (str_xx(i,j,kk) * (v(i,j,k)-v(i,j-1,k))*IdyT(i,j,0)))
        + 0.25_rt*(( (str_xy(i,j,kk) *
               (((u(i,j+1,k)-u(i,j,k))*IdyBu(i,j,0))
              + ((v(i+1,j,k)-v(i,j,k))*IdxBu(i,j,0))))
            + (str_xy(i-1,j-1,kk) *
               (((u(i-1,j,k)-u(i-1,j-1,k))*IdyBu(i-1,j-1,0))
              + ((v(i,j-1,k)-v(i-1,j-1,k))*IdxBu(i-1,j-1,0)))) )
          + ( (str_xy(i-1,j,kk) *
               (((u(i-1,j+1,k)-u(i-1,j,k))*IdyBu(i-1,j,0))
              + ((v(i,j,k)-v(i-1,j,k))*IdxBu(i-1,j,0))))
            + (str_xy(i,j-1,kk) *
               (((u(i,j,k)-u(i,j-1,k))*IdyBu(i,j-1,0))
              + ((v(i+1,j-1,k)-v(i,j-1,k))*IdxBu(i,j-1,0)))) ) ) );
}

///
/// @brief Diagnoses str_xx*d_x u - str_yy*d_y v + str_xy*(d_y u + d_x v) at one h-point from the volume
///        transports.
///
/// @param i, j, k, kk The point, with k the absolute layer and kk its index in the k-block
/// @param H_to_RZ     A constant that translates thickness units to the units of mass per unit area
///                    [R Z H-1 ~> kg m-3 or 1]
/// @param h_neglect   thickness so small it can be lost in roundoff and so neglected [H ~> m or kg m-2]
/// @param str_xx      The diagonal term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2]
/// @param str_xy      The cross term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2]
/// @param uh          The zonal volume transport [H L2 T-1 ~> m3 s-1]
/// @param vh          The meridional volume transport [H L2 T-1 ~> m3 s-1]
/// @param h_u         Thickness interpolated to u points [H ~> m or kg m-2]
/// @param h_v         Thickness interpolated to v points [H ~> m or kg m-2]
/// @param IareaT      1/areaT [L-2 ~> m-2]
/// @param dy2h        Pre-calculated dy^2 at h points [L2 ~> m2]
/// @param dx2h        Pre-calculated dx^2 at h points [L2 ~> m2]
/// @param dxCu        dxCu is delta x at u points [L ~> m]
/// @param IdyCu       1/dyCu [L-1 ~> m-1]
/// @param IareaCu     The masked inverse areas of u-grid cells [L-2 ~> m-2]
/// @param dyCv        dyCv is delta y at v points [L ~> m]
/// @param IdxCv       1/dxCv [L-1 ~> m-1]
/// @param IareaCv     The masked inverse areas of v-grid cells [L-2 ~> m-2]
/// @param dx2q        Pre-calculated dx^2 at q points [L2 ~> m2]
/// @param dy2q        Pre-calculated dy^2 at q points [L2 ~> m2]
/// @return            The work done by the stresses [R Z L2 T-3 ~> W m-2]
///
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
static Real hor_visc_frict_work(const int i, const int j, const int k, const int kk,
                                Real const H_to_RZ,
                                Real const h_neglect,
                                Array4<const Real> const& str_xx,
                                Array4<const Real> const& str_xy,
                                Array4<const Real> const& uh,
                                Array4<const Real> const& vh,
                                Array4<const Real> const& h_u,
                                Array4<const Real> const& h_v,
                                Array4<const Real> const& IareaT,
                                Array4<const Real> const& dy2h,
                                Array4<const Real> const& dx2h,
                                Array4<const Real> const& dxCu,
                                Array4<const Real> const& IdyCu,
                                Array4<const Real> const& IareaCu,
                                Array4<const Real> const& dyCv,
                                Array4<const Real> const& IdxCv,
                                Array4<const Real> const& IareaCv,
                                Array4<const Real> const& dx2q,
                                Array4<const Real> const& dy2q) noexcept
{
    return H_to_RZ * IareaT(i,j,0) * (
        ((str_xx(i,j,kk)*dy2h(i,j,0) * (
              (uh(i,j,k)*dxCu(i,j,0)*IdyCu(i,j,0)*IareaCu(i,j,0)/(h_u(i,j,kk)+h_neglect))
            - (uh(i-1,j,k)*dxCu(i-1,j,0)*IdyCu(i-1,j,0)*IareaCu(i-1,j,0)/(h_u(i-1,j,kk)+h_neglect)) ) )
        - (str_xx(i,j,kk)*dx2h(i,j,0) * (
              (vh(i,j,k)*dyCv(i,j,0)*IdxCv(i,j,0)*IareaCv(i,j,0)/(h_v(i,j,kk)+h_neglect))
            - (vh(i,j-1,k)*dyCv(i,j-1,0)*IdxCv(i,j-1,0)*IareaCv(i,j-1,0)/(h_v(i,j-1,kk)+h_neglect)) ) ))
        + (0.25_rt*(((str_xy(i,j,kk)*(
                 (dx2q(i,j,0)*((uh(i,j+1,k)*IareaCu(i,j+1,0)/(h_u(i,j+1,kk)+h_neglect))
                              - (uh(i,j,k)*IareaCu(i,j,0)/(h_u(i,j,kk)+h_neglect))))
               + (dy2q(i,j,0)*((vh(i+1,j,k)*IareaCv(i+1,j,0)/(h_v(i+1,j,kk)+h_neglect))
                              - (vh(i,j,k)*IareaCv(i,j,0)/(h_v(i,j,kk)+h_neglect)))) ))
            +(str_xy(i-1,j-1,kk)*(
                 (dx2q(i-1,j-1,0)*((uh(i-1,j,k)*IareaCu(i-1,j,0)/(h_u(i-1,j,kk)+h_neglect))
                                  - (uh(i-1,j-1,k)*IareaCu(i-1,j-1,0)/(h_u(i-1,j-1,kk)+h_neglect))))
               + (dy2q(i-1,j-1,0)*((vh(i,j-1,k)*IareaCv(i,j-1,0)/(h_v(i,j-1,kk)+h_neglect))
                                  - (vh(i-1,j-1,k)*IareaCv(i-1,j-1,0)/(h_v(i-1,j-1,kk)+h_neglect)))) )) )
           +((str_xy(i-1,j,kk)*(
                 (dx2q(i-1,j,0)*((uh(i-1,j+1,k)*IareaCu(i-1,j+1,0)/(h_u(i-1,j+1,kk)+h_neglect))
                                - (uh(i-1,j,k)*IareaCu(i-1,j,0)/(h_u(i-1,j,kk)+h_neglect))))
               + (dy2q(i-1,j,0)*((vh(i,j,k)*IareaCv(i,j,0)/(h_v(i,j,kk)+h_neglect))
                                - (vh(i-1,j,k)*IareaCv(i-1,j,0)/(h_v(i-1,j,kk)+h_neglect)))) ))
            +(str_xy(i,j-1,kk)*(
                 (dx2q(i,j-1,0)*((uh(i,j,k)*IareaCu(i,j,0)/(h_u(i,j,kk)+h_neglect))
                                - (uh(i,j-1,k)*IareaCu(i,j-1,0)/(h_u(i,j-1,kk)+h_neglect))))
               + (dy2q(i,j-1,0)*((vh(i+1,j-1,k)*IareaCv(i+1,j-1,0)/(h_v(i+1,j-1,kk)+h_neglect))
                                - (vh(i,j-1,k)*IareaCv(i,j-1,0)/(h_v(i,j-1,kk)+h_neglect)))) )) ) )) );
}

///
/// @brief Accumulates the work done by lateral friction into the vertically integrated MEKE source
///        terms, adjusting for any energy loss seen as a reduction in the (biharmonic) frictional
///        source term.
///
/// @tparam old_RoScl          If true, use the Rossby number function of answer dates 20190102-20241130
/// @tparam has_mom_src_bh     If true, MEKE_mom_src_bh is allocated
/// @param bxH_2d              The h-point iteration box of this tile, restricted to a single level
/// @param kmax                The number of layers in the k-block
/// @param kstart              The first absolute k-layer of the current k-block
/// @param s_to_T              A constant that translates time intervals in seconds to the units of time
///                            [T s-1 ~> 1]
/// @param backscat_subround   The ratio of f over Shear_mag that is so small that the backscatter
///                            calculation gives the same value as if f were 0 [nondim].
/// @param backscatter_Ro_pow  Power in Rossby number function for backscatter [nondim]
/// @param backscatter_Ro_c    Coefficient in Rossby number function for backscatter [nondim]
/// @param CoriolisBu          The Coriolis parameter at corner points [T-1 ~> s-1]
/// @param sh_xx               horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
/// @param sh_xy               horizontal shearing strain (du/dy + dv/dx) including metric terms [T-1 ~> s-1]
/// @param FrictWork           work done by MKE dissipation mechanisms [R Z L2 T-3 ~> W m-2]
/// @param FrictWork_bh        work done by the biharmonic MKE dissipation mechanisms [R Z L2 T-3 ~> W m-2]
/// @param MEKE_mom_src        MEKE source from lateral friction in the momentum equations
///                            [R Z L2 T-3 ~> W m-2]
/// @param MEKE_mom_src_bh     MEKE source from the biharmonic part of the lateral friction in the momentum
///                            equations [R Z L2 T-3 ~> W m-2]
///
template <bool old_RoScl, bool has_mom_src_bh>
void hor_visc_meke_mom_src_loop(const Box& bxH_2d,
                                const int kmax,
                                const int kstart,
                                const Real s_to_T,
                                const Real backscat_subround,
                                const Real backscatter_Ro_pow,
                                const Real backscatter_Ro_c,
                                Array4<const Real> const& CoriolisBu,
                                Array4<const Real> const& sh_xx,
                                Array4<const Real> const& sh_xy,
                                Array4<const Real> const& FrictWork,
                                Array4<const Real> const& FrictWork_bh,
                                Array4<Real>       const& MEKE_mom_src,
                                Array4<Real>       const& MEKE_mom_src_bh)
{
    ParallelFor(bxH_2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        for (int kk = 0; kk < kmax; ++kk) {
            const int k = kstart + kk - 1;
            /// abs(f) at h-point for MEKE source term [T-1 ~> s-1]
            Real FatH = 0.25_rt*( (amrex::Math::abs(CoriolisBu(i-1,j-1,0)) + amrex::Math::abs(CoriolisBu(i,j,0))) +
                               (amrex::Math::abs(CoriolisBu(i-1,j,0)) + amrex::Math::abs(CoriolisBu(i,j-1,0))) );
            /// Shear_mag value in backscatter [T-1 ~> s-1]
            Real Shear_mag_bc = std::sqrt(sh_xx(i,j,kk) * sh_xx(i,j,kk) +
                0.25_rt*(((sh_xy(i-1,j-1,kk)*sh_xy(i-1,j-1,kk)) + (sh_xy(i,j,kk)*sh_xy(i,j,kk))) +
                      ((sh_xy(i-1,j,kk)*sh_xy(i-1,j,kk)) + (sh_xy(i,j-1,kk)*sh_xy(i,j-1,kk)))));
            /// The scaling function for MEKE source term [nondim]
            Real RoScl;
            if constexpr (old_RoScl) {
                FatH = std::pow(s_to_T*FatH, backscatter_Ro_pow); /// f^n
                /// Note the hard-coded dimensional constant in the following line that can not
                /// be rescaled for dimensional consistency.
                Shear_mag_bc = ((std::pow(s_to_T * Shear_mag_bc, backscatter_Ro_pow)) + 1.e-30_rt)
                            * backscatter_Ro_c; /// c * D^n
                /// The Rossby number function is g(Ro) = 1/(1+c.Ro^n)
                /// RoScl = 1 - g(Ro)
                RoScl = Shear_mag_bc / (FatH + Shear_mag_bc); /// = 1 - f^n/(f^n+c*D^n)
            } else {
                if (FatH <= backscat_subround*Shear_mag_bc) {
                    RoScl = 1.0_rt;
                } else {
                    /// The ratio of shear over the absolute value of f raised to some power and rescaled [nondim]
                    Real const Sh_F_pow = backscatter_Ro_c * std::pow(Shear_mag_bc / FatH, backscatter_Ro_pow);
                    RoScl = Sh_F_pow / (1.0_rt + Sh_F_pow); /// = 1 - f^n/(f^n+c*D^n)
                }
            }

            MEKE_mom_src(i,j,0) = MEKE_mom_src(i,j,0) + (FrictWork(i,j,k) - RoScl*FrictWork_bh(i,j,k));

            if constexpr (has_mom_src_bh) {
                MEKE_mom_src_bh(i,j,0) = MEKE_mom_src_bh(i,j,0)
                    + (FrictWork_bh(i,j,k) - RoScl * FrictWork_bh(i,j,k));
            }
        }
    });
}

void hor_visc_stress_div(const Box& bxH,
                         const Box& bxU,
                         const Box& bxV,
                         const Box& bxHstr,
                         const Box& bxQstr,
                         const hv_opts_CS_C& CS,
                         Array4<Real>       const& MEKE_mom_src,
                         Array4<Real>       const& MEKE_mom_src_bh,
                         Array4<Real>       const& MEKE_GME_snk,
                         const Real backscatter_Ro_pow,
                         const Real backscatter_Ro_c,
                         Array4<Real>       const& skeb_diss,
                         const Real skeb_frict_coef,
                         Array4<const Real> const& mask2dBu,
                         Array4<const Real> const& IdxCu,
                         Array4<const Real> const& IdyCu,
                         Array4<const Real> const& IareaCu,
                         Array4<const Real> const& IdyCv,
                         Array4<const Real> const& IdxCv,
                         Array4<const Real> const& IareaCv,
                         Array4<const Real> const& IdxT,
                         Array4<const Real> const& IdyT,
                         Array4<const Real> const& IdyBu,
                         Array4<const Real> const& IdxBu,
                         Array4<const Real> const& IareaT,
                         Array4<const Real> const& dxCu,
                         Array4<const Real> const& dyCv,
                         Array4<const Real> const& CoriolisBu,
                         Array4<const Real> const& reduction_xx,
                         Array4<const Real> const& reduction_xy,
                         Array4<const Real> const& dx2h,
                         Array4<const Real> const& dy2h,
                         Array4<const Real> const& dx2q,
                         Array4<const Real> const& dy2q,
                         OceanOBC* OBC,
                         const int kstart,
                         const Real h_neglect,
                         const bool apply_OBC,
                         const bool find_FrictWork,
                         const bool find_FrictWork_bh,
                         const Real backscat_subround,
                         const bool skeb_use_frict,
                         Array4<const Real> const& u,
                         Array4<const Real> const& v,
                         Array4<const Real> const& h,
                         Array4<const Real> const& uh,
                         Array4<const Real> const& vh,
                         Array4<Real>       const& diffu,
                         Array4<Real>       const& diffv,
                         Array4<const Real> const& str_xx_GME,
                         Array4<const Real> const& str_xy_GME,
                         Array4<Real>       const& FrictWork,
                         Array4<Real>       const& FrictWork_bh,
                         Array4<Real>       const& FrictWork_GME,
                         Array4<const Real> const& visc_limit_h_flag,
                         Array4<const Real> const& h_u,
                         Array4<const Real> const& h_v,
                         Array4<const Real> const& sh_xx,
                         Array4<Real>       const& str_xx,
                         Array4<const Real> const& bhstr_xx,
                         Array4<const Real> const& sh_xy,
                         Array4<Real>       const& str_xy,
                         Array4<const Real> const& bhstr_xy,
                         Array4<const Real> const& hq,
                         const Real H_to_RZ,
                         const Real s_to_T)
{
    BL_PROFILE("hor_visc_stress_div");

    if (OBC != nullptr) {
        AMREX_ABORT_LOC("OBC pointer provided but not yet implemented");
    }

    const int kmax = bxH.length(2);
    Box bxH_2d = bxH;
    bxH_2d.setSmall(2, 0);
    bxH_2d.setBig(2, 0);

    if (CS.use_GME) {

        /// This changes the units of str_xx from [L2 T-2 ~> m2 s-2] to [H L2 T-2 ~> m3 s-2 or kg s-2].
        ParallelFor(bxHstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            const int k = kstart + kk - 1;
            str_xx(i,j,kk) = (str_xx(i,j,kk) + str_xx_GME(i,j,kk)) * (h(i,j,k) * reduction_xx(i,j,0));
        });

        /// This adds in GME and changes the units of str_xx from [L2 T-2 ~> m2 s-2] to [H L2 T-2 ~> m3 s-2 or kg s-2].
        if (CS.no_slip) {
            ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                str_xy(i,j,kk) = (str_xy(i,j,kk) + str_xy_GME(i,j,kk)) * (hq(i,j,kk) * reduction_xy(i,j,0));
            });
        } else {
            ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                str_xy(i,j,kk) = (str_xy(i,j,kk) + str_xy_GME(i,j,kk)) *
                                 (hq(i,j,kk) * mask2dBu(i,j,0) * reduction_xy(i,j,0));
            });
        }
    } else { /// .not. use_GME
        /// This changes the units of str_xx from [L2 T-2 ~> m2 s-2] to [H L2 T-2 ~> m3 s-2 or kg s-2].
        ParallelFor(bxHstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            const int k = kstart + kk - 1;
            str_xx(i,j,kk) = str_xx(i,j,kk) * (h(i,j,k) * reduction_xx(i,j,0));
        });

        /// This changes the units of str_xy from [L2 T-2 ~> m2 s-2] to [H L2 T-2 ~> m3 s-2 or kg s-2].
        if (CS.no_slip) {
            ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                str_xy(i,j,kk) = str_xy(i,j,kk) * (hq(i,j,kk) * reduction_xy(i,j,0));
            });
        } else {
            ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                str_xy(i,j,kk) = str_xy(i,j,kk) * (hq(i,j,kk) * mask2dBu(i,j,0) * reduction_xy(i,j,0));
            });
        }
    } /// use_GME

    /// Evaluate 1/h x.Div(h Grad u) or the biharmonic equivalent.
    ParallelFor(bxU, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
    {
        const int k = kstart + kk - 1;
        diffu(i,j,k) = ((IdxCu(i,j,0)*((dx2q(i,j-1,0)*str_xy(i,j-1,kk)) - (dx2q(i,j,0)*str_xy(i,j,kk))) +
                         IdyCu(i,j,0)*((dy2h(i,j,0)*str_xx(i,j,kk)) - (dy2h(i+1,j,0)*str_xx(i+1,j,kk)))) *
                       IareaCu(i,j,0)) / (h_u(i,j,kk) + h_neglect);
    });

    /*
    if (apply_OBC) {
        /// This is not the right boundary condition. If all the masking of tendencies are done
        /// correctly later then eliminating this block should not change answers.
        for (int n = 0; n < OBC->number_of_segments; ++n) {
            auto const& seg = OBC->segment[n];
            if (seg.is_E_or_W) {
                const int I = seg.HI.IsdB;
                const Box seg_bx(IntVect(I, seg.HI.jsd, 0), IntVect(I, seg.HI.jed, kmax-1));
                ParallelFor(seg_bx, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                {
                    const int k = kstart + kk - 1;
                    diffu(i,j,k) = 0._rt;
                });
            }
        }
    }
    */

    /// Evaluate 1/h y.Div(h Grad u) or the biharmonic equivalent.
    ParallelFor(bxV, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
    {
        const int k = kstart + kk - 1;
        diffv(i,j,k) = ((IdyCv(i,j,0)*((dy2q(i-1,j,0)*str_xy(i-1,j,kk)) - (dy2q(i,j,0)*str_xy(i,j,kk))) -
                         IdxCv(i,j,0)*((dx2h(i,j,0)*str_xx(i,j,kk)) - (dx2h(i,j+1,0)*str_xx(i,j+1,kk)))) *
                       IareaCv(i,j,0)) / (h_v(i,j,kk) + h_neglect);
    });

    /*
    if (apply_OBC) {
        /// This is not the right boundary condition. If all the masking of tendencies are done
        /// correctly later then eliminating this block should not change answers.
        for (int n = 0; n < OBC->number_of_segments; ++n) {
            auto const& seg = OBC->segment[n];
            if (seg.is_N_or_S) {
                const int J = seg.HI.JsdB;
                const Box seg_bx(IntVect(seg.HI.isd, J, 0), IntVect(seg.HI.ied, J, kmax-1));
                ParallelFor(seg_bx, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                {
                    const int k = kstart + kk - 1;
                    diffv(i,j,k) = 0._rt;
                });
            }
        }
    }
    */


    if (find_FrictWork) {

        if (CS.FrictWork_bug) {
            /// Diagnose   str_xx*d_x u - str_yy*d_y v + str_xy*(d_y u + d_x v)
            /// This is the old formulation that includes energy diffusion
            ParallelFor(bxH, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                FrictWork(i,j,k) = hor_visc_frict_work_bug(i, j, k, kk, H_to_RZ, str_xx, str_xy, u, v,
                                                           IdxT, IdyT, IdyBu, IdxBu);
            });
        } else {
            ParallelFor(bxH, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                FrictWork(i,j,k) = hor_visc_frict_work(i, j, k, kk, H_to_RZ, h_neglect, str_xx, str_xy, uh, vh,
                                                       h_u, h_v, IareaT, dy2h, dx2h, dxCu, IdyCu, IareaCu,
                                                       dyCv, IdxCv, IareaCv, dx2q, dy2q);
            });
        }

        if (CS.EY24_EBT_BS) {
            ParallelFor(bxH, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                FrictWork(i,j,k) = (1._rt - visc_limit_h_flag(i,j,k)) * FrictWork(i,j,k);
            });
        }
    }

    if (find_FrictWork_bh) {
        if (CS.FrictWork_bug) {
            /// Diagnose   bhstr_xx*d_x u - bhstr_yy*d_y v + bhstr_xy*(d_y u + d_x v)
            /// This is the old formulation that includes energy diffusion !cyc
            ParallelFor(bxH, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                FrictWork_bh(i,j,k) = hor_visc_frict_work_bug(i, j, k, kk, H_to_RZ, bhstr_xx, bhstr_xy, u, v,
                                                              IdxT, IdyT, IdyBu, IdxBu);
            });
        } else {
            ParallelFor(bxH, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                /// Diagnose   bhstr_xx*d_x u - bhstr_yy*d_y v + bhstr_xy*(d_y u + d_x v)
                FrictWork_bh(i,j,k) = hor_visc_frict_work(i, j, k, kk, H_to_RZ, h_neglect, bhstr_xx, bhstr_xy,
                                                          uh, vh, h_u, h_v, IareaT, dy2h, dx2h, dxCu, IdyCu,
                                                          IareaCu, dyCv, IdxCv, IareaCv, dx2q, dy2q);
            });
        }

        if (CS.EY24_EBT_BS) {
            ParallelFor(bxH, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                FrictWork_bh(i,j,k) = (1._rt - visc_limit_h_flag(i,j,k)) * FrictWork_bh(i,j,k);
            });
        }
    }

    if (CS.use_GME) {
        if (CS.FrictWork_bug) {
            ParallelFor(bxH, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                /// Diagnose   str_xx_GME*d_x u - str_yy_GME*d_y v + str_xy_GME*(d_y u + d_x v)
                /// This is the old formulation that includes energy diffusion
                FrictWork_GME(i,j,k) = hor_visc_frict_work_bug(i, j, k, kk, H_to_RZ, str_xx_GME, str_xy_GME, u, v,
                                                               IdxT, IdyT, IdyBu, IdxBu);
            });
        } else {
            ParallelFor(bxH, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                FrictWork_GME(i,j,k) = hor_visc_frict_work(i, j, k, kk, H_to_RZ, h_neglect, str_xx_GME, str_xy_GME,
                                                           uh, vh, h_u, h_v, IareaT, dy2h, dx2h, dxCu, IdyCu,
                                                           IareaCu, dyCv, IdxCv, IareaCv, dx2q, dy2q);
            });
        }
        if (find_FrictWork && (MEKE_mom_src.p != nullptr) && (MEKE_GME_snk.p != nullptr)) {
            ParallelFor(bxH_2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                for (int kk = 0; kk < kmax; ++kk) {
                    const int k = kstart + kk - 1;
                    MEKE_GME_snk(i,j,0) = MEKE_GME_snk(i,j,0) + FrictWork_GME(i,j,k);
                }
            });
        }
    }

    if (skeb_use_frict) {
        ParallelFor(bxH, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            const int k = kstart + kk - 1;
            /// Note that the sign convention is FrictWork < 0 means energy dissipation.
            skeb_diss(i,j,k) = skeb_diss(i,j,k) - skeb_frict_coef *
                                     FrictWork(i,j,k) / (H_to_RZ * (h(i,j,k) + h_neglect));
        });
    }

    /// Make a similar calculation as for FrictWork above but accumulating into
    /// the vertically integrated MEKE source term, and adjusting for any
    /// energy loss seen as a reduction in the (biharmonic) frictional source term.
    if (find_FrictWork && (MEKE_mom_src.p != nullptr)) {
        if (backscatter_Ro_c != 0._rt) {
            const bool old_RoScl = (CS.answer_date > 20190101) && (CS.answer_date < 20241201);
            const bool has_mom_src_bh = (MEKE_mom_src_bh.p != nullptr);
            if (old_RoScl) {
                if (has_mom_src_bh) {
                    hor_visc_meke_mom_src_loop<true, true>(bxH_2d, kmax, kstart, s_to_T, backscat_subround,
                                                           backscatter_Ro_pow, backscatter_Ro_c, CoriolisBu,
                                                           sh_xx, sh_xy, FrictWork, FrictWork_bh,
                                                           MEKE_mom_src, MEKE_mom_src_bh);
                } else {
                    hor_visc_meke_mom_src_loop<true, false>(bxH_2d, kmax, kstart, s_to_T, backscat_subround,
                                                            backscatter_Ro_pow, backscatter_Ro_c, CoriolisBu,
                                                            sh_xx, sh_xy, FrictWork, FrictWork_bh,
                                                            MEKE_mom_src, MEKE_mom_src_bh);
                }
            } else {
                if (has_mom_src_bh) {
                    hor_visc_meke_mom_src_loop<false, true>(bxH_2d, kmax, kstart, s_to_T, backscat_subround,
                                                            backscatter_Ro_pow, backscatter_Ro_c, CoriolisBu,
                                                            sh_xx, sh_xy, FrictWork, FrictWork_bh,
                                                            MEKE_mom_src, MEKE_mom_src_bh);
                } else {
                    hor_visc_meke_mom_src_loop<false, false>(bxH_2d, kmax, kstart, s_to_T, backscat_subround,
                                                             backscatter_Ro_pow, backscatter_Ro_c, CoriolisBu,
                                                             sh_xx, sh_xy, FrictWork, FrictWork_bh,
                                                             MEKE_mom_src, MEKE_mom_src_bh);
                }
            }
        } else {
            ParallelFor(bxH_2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                for (int kk = 0; kk < kmax; ++kk) {
                    const int k = kstart + kk - 1;
                    MEKE_mom_src(i,j,0) = MEKE_mom_src(i,j,0) + FrictWork(i,j,k);
                }
            });

            if (MEKE_mom_src_bh.p != nullptr) {
                ParallelFor(bxH_2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
                {
                    for (int kk = 0; kk < kmax; ++kk) {
                        const int k = kstart + kk - 1;
                        MEKE_mom_src_bh(i,j,0) = MEKE_mom_src_bh(i,j,0) + FrictWork_bh(i,j,k);
                    }
                });
            }
        } /// MEKE_a%backscatter_Ro_c

    } /// find_FrictWork and associated(mom_src)
}

void hor_visc_h_coeffs_b(const Box& bxH,
                         const Box& bxQ,
                         const Box& bxKh,
                         const Box& bxHstr,
                         const Box& bxQstr,
                         const hv_opts_CS_C& CS,
                         Array4<const Real> const& MEKE_Au,
                         Array4<const Real> const& MEKE_Ku,
                         Array4<const Real> const& BS_struct,
                         Array4<const Real> const& IdyCu,
                         Array4<const Real> const& IdxCv,
                         Array4<const Real> const& IdyCv,
                         Array4<const Real> const& IdxCu,
                         Array4<const Real> const& mask2dBu,
                         Array4<const Real> const& mask2dCu,
                         Array4<const Real> const& mask2dCv,
                         Array4<const Real> const& reduction_xx,
                         Array4<const Real> const& Ah_Max_xx,
                         Array4<const Real> const& Ah_Max_xx_KS,
                         Array4<const Real> const& grid_sp_h3,
                         Array4<const Real> const& dx_dyT,
                         Array4<const Real> const& dy_dxT,
                         Array4<const Real> const& dx_dyBu,
                         Array4<const Real> const& dy_dxBu,
                         Array4<const Real> const& Re_Ah_const_xx,
                         OceanOBC* OBC,
                         const int kstart,
                         const Real h_neglect,
                         const Real h_neglect3,
                         const bool apply_OBC,
                         const bool use_kh_struct,
                         const bool use_MEKE_Au,
                         Array4<const Real> const& u,
                         Array4<const Real> const& v,
                         Array4<const Real> const& h,
                         Array4<Real>       const& Ah_h,
                         Array4<Real>       const& Kh_h,
                         Array4<Real>       const& grid_Re_Ah,
                         Array4<Real>       const& visc_limit_h,
                         Array4<Real>       const& visc_limit_h_flag,
                         Array4<Real>       const& visc_limit_h_frac,
                         Array4<Real>       const& BS_coeff_h,
                         Array4<const Real> const& Del2u,
                         Array4<const Real> const& h_u,
                         Array4<const Real> const& Del2v,
                         Array4<const Real> const& h_v,
                         Array4<const Real> const& sh_xx,
                         Array4<const Real> const& sh_xx_smooth,
                         Array4<Real>       const& str_xx,
                         Array4<Real>       const& bhstr_xx,
                         Array4<const Real> const& m_leithy,
                         Array4<Real>       const& dDel2vdx,
                         Array4<Real>       const& dDel2udy,
                         Array4<const Real> const& sh_xy,
                         Array4<Real>       const& hq,
                         Array4<Real>       const& Ah,
                         Array4<Real>       const& Kh,
                         Array4<Real>       const& Shear_mag,
                         Array4<Real>       const& hrat_min,
                         Array4<const Real> const& visc_bound_rem)
{
    BL_PROFILE("hor_visc_h_coeffs_b");

    if (OBC != nullptr) {
        AMREX_ABORT_LOC("OBC pointer provided but not yet implemented");
    }

    if (CS.biharmonic) {

        if (use_MEKE_Au) {
            /// *Add* the MEKE contribution
            ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                Ah(i,j,kk) = Ah(i,j,kk) + MEKE_Au(i,j,0);
            });
        }

        if (CS.Re_Ah > 0.0_rt) {
            ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                Real const u_sum = u(i,j,k)+u(i-1,j,k);
                Real const v_sum = v(i,j,k)+v(i,j-1,k);
                /// Local kinetic energy [L2 T-2 ~> m2 s-2]
                Real const KE = 0.125_rt*((u_sum*u_sum) + (v_sum*v_sum));
                Ah(i,j,kk) = std::sqrt(KE) * Re_Ah_const_xx(i,j,0);
            });
        }

        if (CS.bound_Ah) {
            if (CS.bound_Kh) {
                ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                {
                    Ah(i,j,kk) = amrex::min(Ah(i,j,kk), visc_bound_rem(i,j,kk) * hrat_min(i,j,kk) * Ah_Max_xx(i,j,0));
                });
            } else {
                ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                {
                    Ah(i,j,kk) = amrex::min(Ah(i,j,kk), hrat_min(i,j,kk) * Ah_Max_xx(i,j,0));
                });
            }
        }

        if (CS.EY24_EBT_BS) {
            const Real KS_coef = CS.KS_coef;
            ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                Real const tmp = KS_coef * hrat_min(i,j,kk) * Ah_Max_xx_KS(i,j,0);
                visc_limit_h(i,j,k) = tmp;
                visc_limit_h_frac(i,j,k) = Ah(i,j,kk) / (KS_coef * hrat_min(i,j,kk) * Ah_Max_xx_KS(i,j,0));
                if (Ah(i,j,kk) >= tmp) {
                    visc_limit_h_flag(i,j,k) = 1._rt;
                }
            });
        }

        if ((CS.do_Ah_h) || CS.debug || CS.use_Leithy) {
            ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                Ah_h(i,j,k) = Ah(i,j,kk);
            });
        }

        if (CS.use_Leithy) {
            /// Compute Leith+E Kh after bounds have been applied to Ah
            /// and after it has been smoothed. Kh = -m_leithy * Ah
            ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                Kh(i,j,kk) = -m_leithy(i,j,kk) * Ah(i,j,kk);
                Kh_h(i,j,k) = Kh(i,j,kk);
            });
        }

        if (CS.do_grid_Re_Ah) {
            const Real min_grid_Ah = CS.min_grid_Ah;
            ParallelFor(bxH, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                Real const u_sum = u(i,j,k) + u(i-1,j,k);
                Real const v_sum = v(i,j,k) + v(i,j-1,k);
                /// Local kinetic energy [L2 T-2 ~> m2 s-2]
                Real const KE = 0.125_rt * ((u_sum*u_sum) + (v_sum*v_sum));
                /// Biharmonic viscosity bound by grid [L4 T-1 ~> m4 s-1]
                Real const grid_Ah = amrex::max(Ah(i,j,kk), min_grid_Ah);
                grid_Re_Ah(i,j,k) = (std::sqrt(KE) * grid_sp_h3(i,j,0)) / grid_Ah;
            });
        }

        if (CS.use_Leithy) {
            ParallelFor(bxHstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                /// dy-weighted Laplacian(u) diff in x [L-2 T-1 ~> m-2 s-1]
                Real const d_del2u = (IdyCu(i,j,0) * Del2u(i,j,kk)) - (IdyCu(i-1,j,0) * Del2u(i-1,j,kk));
                /// dx-weighted Laplacian(v) diff in y [L-2 T-1 ~> m-2 s-1]
                Real const d_del2v = (IdxCv(i,j,0) * Del2v(i,j,kk)) - (IdxCv(i,j-1,0) * Del2v(i,j-1,kk));
                /// Stress tensor update [L2 T-2 ~> m2 s-2]
                Real const d_str = Ah(i,j,kk) * ((dy_dxT(i,j,0) * d_del2u) - (dx_dyT(i,j,0) * d_del2v));

                str_xx(i,j,kk) = str_xx(i,j,kk) + d_str;

                str_xx(i,j,kk) = str_xx(i,j,kk) - Kh(i,j,kk) * sh_xx_smooth(i,j,kk);

                /// Keep a copy of the biharmonic contribution for backscatter parameterization
                bhstr_xx(i,j,kk) = d_str * (h(i,j,k) * reduction_xx(i,j,0));
            });
        } else {
            ParallelFor(bxHstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                /// dy-weighted Laplacian(u) diff in x [L-2 T-1 ~> m-2 s-1]
                Real const d_del2u = (IdyCu(i,j,0) * Del2u(i,j,kk)) - (IdyCu(i-1,j,0) * Del2u(i-1,j,kk));
                /// dx-weighted Laplacian(v) diff in y [L-2 T-1 ~> m-2 s-1]
                Real const d_del2v = (IdxCv(i,j,0) * Del2v(i,j,kk)) - (IdxCv(i,j-1,0) * Del2v(i,j-1,kk));
                /// Stress tensor update [L2 T-2 ~> m2 s-2]
                Real const d_str = Ah(i,j,kk) * ((dy_dxT(i,j,0) * d_del2u) - (dx_dyT(i,j,0) * d_del2v));

                str_xx(i,j,kk) = str_xx(i,j,kk) + d_str;

                /// Keep a copy of the biharmonic contribution for backscatter parameterization
                bhstr_xx(i,j,kk) = d_str * (h(i,j,k) * reduction_xx(i,j,0));
            });
        }
    } /// Get biharmonic coefficient at h points and biharmonic part of str_xx

    /// Backscatter using MEKE
    if (CS.EY24_EBT_BS) {
        hor_visc_backscatter_h(bxHstr, CS, MEKE_Ku, BS_struct, use_kh_struct, kstart,
                               visc_limit_h_flag, sh_xx, str_xx, BS_coeff_h);
    } /// Backscatter

    if (CS.biharmonic) {
        /// Gradient of Laplacian, for use in bi-harmonic term
        ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            dDel2vdx(i,j,kk) = dy_dxBu(i,j,0)*((Del2v(i+1,j,kk)*IdyCv(i+1,j,0)) - (Del2v(i,j,kk)*IdyCv(i,j,0)));
            dDel2udy(i,j,kk) = dx_dyBu(i,j,0)*((Del2u(i,j+1,kk)*IdxCu(i,j+1,0)) - (Del2u(i,j,kk)*IdxCu(i,j,0)));
        });
        /// Adjust contributions to shearing strain on open boundaries.
        /*
        if (apply_OBC) {
            if ((OBC->strain_config == OBC_STRAIN_ZERO) ||
                (OBC->strain_config == OBC_STRAIN_FREESLIP)) {
                const int is   = bxH.smallEnd(0);
                const int js   = bxH.smallEnd(1);
                const int Ieq  = bxQ.bigEnd(0);
                const int Jeq  = bxQ.bigEnd(1);
                const int kmax = bxH.length(2);
                const int strain_config = OBC->strain_config;
                for (int n = 0; n < OBC->number_of_segments; ++n) {
                    auto const& seg = OBC->segment[n];
                    const int J = seg.HI.JsdB;
                    const int I = seg.HI.IsdB;
                    if (seg.is_N_or_S && (J >= js-1) && (J <= Jeq)) {
                        const Box seg_bx(IntVect(seg.HI.IsdB, J, 0), IntVect(seg.HI.IedB, J, kmax-1));
                        ParallelFor(seg_bx, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                        {
                            if (strain_config == OBC_STRAIN_ZERO) {
                                dDel2vdx(i,j,kk) = 0._rt; dDel2udy(i,j,kk) = 0._rt;
                            } else if (strain_config == OBC_STRAIN_FREESLIP) {
                                dDel2udy(i,j,kk) = 0._rt;
                            }
                        });
                    } else if (seg.is_E_or_W && (I >= is-1) && (I <= Ieq)) {
                        const Box seg_bx(IntVect(I, seg.HI.JsdB, 0), IntVect(I, seg.HI.JedB, kmax-1));
                        ParallelFor(seg_bx, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                        {
                            if (strain_config == OBC_STRAIN_ZERO) {
                                dDel2vdx(i,j,kk) = 0._rt; dDel2udy(i,j,kk) = 0._rt;
                            } else if (strain_config == OBC_STRAIN_FREESLIP) {
                                dDel2vdx(i,j,kk) = 0._rt;
                            }
                        });
                    }
                }
            }
        }
        */
    }

    if ((CS.Smagorinsky_Kh) || (CS.Smagorinsky_Ah)) {
        ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            /// Square of shearing strain (sh_xy) [T-2 ~> s-2]
            Real const sh_xy_sq = sh_xy(i,j,kk)*sh_xy(i,j,kk);
            /// Square of tension (sh_xx) [T-2 ~> s-2]
            Real const sh_xx_sq = 0.25_rt * ( ((sh_xx(i,j,kk)*sh_xx(i,j,kk)) + (sh_xx(i+1,j+1,kk)*sh_xx(i+1,j+1,kk)))
                                         + ((sh_xx(i,j+1,kk)*sh_xx(i,j+1,kk)) + (sh_xx(i+1,j,kk)*sh_xx(i+1,j,kk))) );
            Shear_mag(i,j,kk) = std::sqrt(sh_xy_sq + sh_xx_sq);
        });
    }

    ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
    {
        /// temporary variables [H2 ~> m2 or kg2 m-4].
        Real const h2uq = 4.0_rt * (h_u(i,j,kk) * h_u(i,j+1,kk));
        Real const h2vq = 4.0_rt * (h_v(i,j,kk) * h_v(i+1,j,kk));
        hq(i,j,kk) = (2.0_rt * (h2uq * h2vq))
            / (h_neglect3 + (h2uq + h2vq) * ((h_u(i,j,kk) + h_u(i,j+1,kk)) + (h_v(i,j,kk) + h_v(i+1,j,kk))));
    });

    if (CS.bound_Ah || CS.bound_Kh) {
        ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            /// Minimum h at the 4 neighboring velocity points [H ~> m]
            Real const h_min = amrex::min(h_u(i,j,kk), h_u(i,j+1,kk), h_v(i,j,kk), h_v(i+1,j,kk));
            hrat_min(i,j,kk) = amrex::min(1.0_rt, h_min / (hq(i,j,kk) + h_neglect));
        });

    }

    /// TODO: GPU??  Are h_[uv] on CPU?  update to hrat_min?
    if (CS.no_slip) {
        const bool no_slip = CS.no_slip;
        ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            if (no_slip && (mask2dBu(i,j,0) < 0.5_rt)) {
                if ((mask2dCu(i,j,0) + mask2dCu(i,j+1,0)) +
                    (mask2dCv(i,j,0) + mask2dCv(i+1,j,0)) > 0.0_rt) {
                    /// This is a coastal vorticity point, so modify hq and hrat_min.

                    /// Thicknesses interpolated by arithmetic means to corner
                    /// points; these are first interpolated to u or v velocity
                    /// points where masks are applied [H ~> m or kg m-2].
                    Real const hu = mask2dCu(i,j,0) * h_u(i,j,kk) + mask2dCu(i,j+1,0) * h_u(i,j+1,kk);
                    Real const hv = mask2dCv(i,j,0) * h_v(i,j,kk) + mask2dCv(i+1,j,0) * h_v(i+1,j,kk);
                    if ((mask2dCu(i,j,0) + mask2dCu(i,j+1,0)) *
                        (mask2dCv(i,j,0) + mask2dCv(i+1,j,0)) == 0.0_rt) {
                        /// Only one of hu and hv is nonzero, so just add them.
                        hq(i,j,kk) = hu + hv;
                        hrat_min(i,j,kk) = 1.0_rt;
                    } else {
                        /// Both hu and hv are nonzero, so take the harmonic mean.
                        hq(i,j,kk) = 2.0_rt * (hu * hv) / ((hu + hv) + h_neglect);
                        hrat_min(i,j,kk) = amrex::min(1.0_rt, amrex::min(hu, hv) / (hq(i,j,kk) + h_neglect) );
                    }
                }
            }
        });
    }
}

///
/// @brief Adds the MEKE contribution (which might be negative) to the Laplacian viscosity at q points.
///
/// @tparam use_kh_struct   If true, the MEKE viscosity contributions are shaped with VarMix%BS_struct
/// @tparam res_scale_MEKE  If true, the MEKE viscosity is scaled by the resolution function
/// @param bxQstr           The q points of this tile where the shear stresses are computed,
///                         [isc-1:iec, jsc-1:jec, 1:kmax]
/// @param kstart           The first absolute k-layer of the current k-block
/// @param MEKE_Ku          The MEKE-derived lateral viscosity coefficient [L2 T-1 ~> m2 s-1]
/// @param Res_fn_q         Non-dimensional function of the ratio the first baroclinic deformation radius to the
///                         grid spacing at q points [nondim]
/// @param BS_struct        Vertical structure function used in backscatter [nondim]
/// @param Kh               Laplacian viscosity (h or q) [L2 T-1 ~> m2 s-1]
///
template <bool use_kh_struct, bool res_scale_MEKE>
void hor_visc_meke_ku_q_loop(const Box& bxQstr,
                             const int kstart,
                             Array4<const Real> const& MEKE_Ku,
                             Array4<const Real> const& Res_fn_q,
                             Array4<const Real> const& BS_struct,
                             Array4<Real>       const& Kh)
{
    ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
    {
        /// A copy of the resolution scaling factor if being applied to MEKE [nondim]. Otherwise = 1.
        Real meke_res_fn = 1._rt;
        if constexpr (res_scale_MEKE) {
            meke_res_fn = Res_fn_q(i,j,0);
        }

        if constexpr (use_kh_struct) {
            const int k = kstart + kk - 1;
            Kh(i,j,kk) = Kh(i,j,kk) + 0.25_rt*( ((MEKE_Ku(i,j,0)*BS_struct(i,j,k)) +
                                       (MEKE_Ku(i+1,j+1,0)*BS_struct(i+1,j+1,k))) +
                                       ((MEKE_Ku(i+1,j,0)*BS_struct(i+1,j,k)) +
                                       (MEKE_Ku(i,j+1,0)*BS_struct(i,j+1,k))) ) * meke_res_fn;
        } else {
            Kh(i,j,kk) = Kh(i,j,kk) + 0.25_rt * (
                (MEKE_Ku(i,j,0) + MEKE_Ku(i+1,j+1,0)) +
                                       (MEKE_Ku(i+1,j,0) +
                                        MEKE_Ku(i,j+1,0)) ) * meke_res_fn;
        }
    });
}

void hor_visc_q_coeffs(const Box& bxQstr,
                       const Box& bxHgrad,
                       const Box& bxQshear,
                       const hv_opts_CS_C& CS,
                       Array4<const Real> const& MEKE_Ku,
                       Array4<const Real> const& MEKE_Au,
                       Array4<const Real> const& Res_fn_q,
                       Array4<const Real> const& BS_struct,
                       Array4<const Real> const& mask2dBu,
                       Array4<const Real> const& Kh_bg_xy,
                       Array4<const Real> const& Ah_bg_xy,
                       Array4<const Real> const& reduction_xy,
                       Array4<const Real> const& Kh_Max_xy,
                       Array4<const Real> const& Ah_Max_xy,
                       Array4<const Real> const& Ah_Max_xy_KS,
                       Array4<const Real> const& n1n2_q,
                       Array4<const Real> const& n1n1_m_n2n2_q,
                       Array4<const Real> const& Laplac2_const_xy,
                       Array4<const Real> const& Biharm6_const_xy,
                       Array4<const Real> const& Laplac3_const_xy,
                       Array4<const Real> const& Biharm_const_xy,
                       Array4<const Real> const& Biharm_const2_xy,
                       Array4<const Real> const& Re_Ah_const_xy,
                       const int kstart,
                       const Real inv_PI3,
                       const Real inv_PI6,
                       const bool find_FrictWork,
                       const bool use_kh_struct,
                       const bool use_MEKE_Ku,
                       const bool use_MEKE_Au,
                       const bool rescale_Kh,
                       Array4<const Real> const& u,
                       Array4<const Real> const& v,
                       Array4<const Real> const& sh_xx_bt,
                       Array4<const Real> const& sh_xy_bt,
                       Array4<const Real> const& GME_effic_h,
                       Array4<const Real> const& GME_effic_q,
                       Array4<const Real> const& KH_u_GME,
                       Array4<const Real> const& KH_v_GME,
                       Array4<Real>       const& GME_coeff_h,
                       Array4<Real>       const& GME_coeff_q,
                       Array4<Real>       const& str_xx_GME,
                       Array4<Real>       const& str_xy_GME,
                       Array4<const Real> const& Ah_h,
                       Array4<Real>       const& Ah_q,
                       Array4<const Real> const& Kh_h,
                       Array4<Real>       const& Kh_q,
                       Array4<Real>       const& vort_xy_q,
                       Array4<Real>       const& sh_xy_q,
                       Array4<Real>       const& visc_limit_q,
                       Array4<Real>       const& visc_limit_q_flag,
                       Array4<Real>       const& visc_limit_q_frac,
                       Array4<Real>       const& BS_coeff_q,
                       Array4<const Real> const& sh_xx,
                       Array4<const Real> const& dDel2vdx,
                       Array4<const Real> const& dDel2udy,
                       Array4<const Real> const& sh_xy,
                       Array4<const Real> const& sh_xy_smooth,
                       Array4<Real>       const& str_xy,
                       Array4<Real>       const& bhstr_xy,
                       Array4<const Real> const& vort_xy,
                       Array4<const Real> const& grad_vort_mag_q,
                       Array4<const Real> const& grad_vort_mag_q_2d,
                       Array4<const Real> const& Del2vort_q,
                       Array4<const Real> const& grad_div_mag_q,
                       Array4<const Real> const& hq,
                       Array4<Real>       const& Ah,
                       Array4<Real>       const& Kh,
                       Array4<const Real> const& Shear_mag,
                       Array4<Real>       const& vert_vort_mag,
                       Array4<const Real> const& hrat_min,
                       Array4<Real>       const& visc_bound_rem)
{
    BL_PROFILE("hor_visc_q_coeffs");

    if (CS.Laplacian) {
        /// Determine the Laplacian viscosity at q points, using the
        /// largest value from several parameterizations. Also get the
        /// Laplacian component of str_xy.

        if ((CS.Leith_Kh) || (CS.Leith_Ah)) {
            if (CS.use_QG_Leith_visc) {
                ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                {
                    /// Vorticity gradient magnitude [L-1 T-1 ~> m-1 s-1]
                    Real const grad_vort = grad_vort_mag_q(i,j,kk) + grad_div_mag_q(i,j,kk);
                    /// QG-based vorticity gradient magnitude [L-1 T-1 ~> m-1 s-1]
                    Real const grad_vort_qg = 3._rt * grad_vort_mag_q_2d(i,j,kk);
                    vert_vort_mag(i,j,kk) = amrex::min(grad_vort, grad_vort_qg);
                });
            } else {
                ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                {
                    vert_vort_mag(i,j,kk) = grad_vort_mag_q(i,j,kk) + grad_div_mag_q(i,j,kk);
                });
            }
        }

        /// Static (pre-computed) background viscosity
        ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            Kh(i,j,kk) = Kh_bg_xy(i,j,0);
        });

        if (CS.Smagorinsky_Kh) {
            if (CS.add_LES_viscosity) {
                ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                {
                    Kh(i,j,kk) = Kh(i,j,kk) + Laplac2_const_xy(i,j,0) * Shear_mag(i,j,kk);
                });
            } else {
                ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                {
                    Kh(i,j,kk) = amrex::max(Kh(i,j,kk), Laplac2_const_xy(i,j,0) * Shear_mag(i,j,kk) );
                });
            }
        }

        if (CS.Leith_Kh) {
            if (CS.add_LES_viscosity) {
                ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                {
                    /// Is this right? -AJA
                    Kh(i,j,kk) = Kh(i,j,kk) + Laplac3_const_xy(i,j,0) * vert_vort_mag(i,j,kk) * inv_PI3;
                });
            } else {
                ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                {
                    Kh(i,j,kk) = amrex::max(Kh(i,j,kk), Laplac3_const_xy(i,j,0) * vert_vort_mag(i,j,kk) * inv_PI3);
                });
            }
        }

        /// All viscosity contributions above are subject to resolution scaling

        if (rescale_Kh) {
            ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                Kh(i,j,kk) = Res_fn_q(i,j,0) * Kh(i,j,kk);
            });
        }

        const Real Kh_bg_min = CS.Kh_bg_min;
        ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            Kh(i,j,kk) = amrex::max(Kh(i,j,kk), Kh_bg_min); /// Place a floor on the viscosity, if desired.
        });

        if (use_MEKE_Ku && !CS.EY24_EBT_BS) {
            if (use_kh_struct) {
                if (CS.res_scale_MEKE) {
                    hor_visc_meke_ku_q_loop<true, true>(bxQstr, kstart, MEKE_Ku, Res_fn_q, BS_struct, Kh);
                } else {
                    hor_visc_meke_ku_q_loop<true, false>(bxQstr, kstart, MEKE_Ku, Res_fn_q, BS_struct, Kh);
                }
            } else {
                if (CS.res_scale_MEKE) {
                    hor_visc_meke_ku_q_loop<false, true>(bxQstr, kstart, MEKE_Ku, Res_fn_q, BS_struct, Kh);
                } else {
                    hor_visc_meke_ku_q_loop<false, false>(bxQstr, kstart, MEKE_Ku, Res_fn_q, BS_struct, Kh);
                }
            }
        }

        if (CS.anisotropic) {
            const Real Kh_aniso = CS.Kh_aniso;
            /// *Add* the shear component of anisotropic viscosity
            ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                Kh(i,j,kk) = Kh(i,j,kk) + Kh_aniso * (n1n2_q(i,j,0)*n1n2_q(i,j,0));
            });
        }

        /// Newer method of bounding for stability
        if ((CS.bound_Kh) && (CS.bound_Ah)) {
            const bool backscatter_underbound = CS.backscatter_underbound;
            ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                visc_bound_rem(i,j,kk) = 1.0_rt;
                /// The local maximum Laplacian viscosity for stability [L2 T-1 ~> m2 s-1]
                Real const Kh_max_here = hrat_min(i,j,kk) * Kh_Max_xy(i,j,0);
                if (Kh(i,j,kk) >= Kh_max_here) {
                    visc_bound_rem(i,j,kk) = 0.0_rt;
                    Kh(i,j,kk) = Kh_max_here;
                } else if ((Kh(i,j,kk) > 0.0_rt) || (backscatter_underbound && (Kh_max_here > 0.0_rt))) {
                    visc_bound_rem(i,j,kk) = 1.0_rt - Kh(i,j,kk) / Kh_max_here;
                }
            });
        } else if (CS.bound_Kh) {
            ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                Kh(i,j,kk) = amrex::min(Kh(i,j,kk), hrat_min(i,j,kk) * Kh_Max_xy(i,j,0));
            });
        }

        if (CS.use_Leithy) {
            /// Leith+E doesn't recompute Kh at q points, it just interpolates it from h to q points
            ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                Kh(i,j,kk) = 0.25_rt * ((Kh_h(i,j,k) + Kh_h(i+1,j+1,k)) + (Kh_h(i,j+1,k) + Kh_h(i+1,j,k)));
            });
        }

        if (CS.do_Kh_q || CS.debug) {
            ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                Kh_q(i,j,k) = Kh(i,j,kk);
            });
        }

        if (CS.do_vort_xy_q) {
            ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                vort_xy_q(i,j,k) = vort_xy(i,j,kk);
            });
        }

        if (CS.do_sh_xy_q) {
            ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                sh_xy_q(i,j,k) = sh_xy(i,j,kk);
            });
        }

        if (!CS.use_Leithy) {
            ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                str_xy(i,j,kk) = -Kh(i,j,kk) * sh_xy(i,j,kk);
            });
        } else {
            ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                str_xy(i,j,kk) = -Kh(i,j,kk) * sh_xy_smooth(i,j,kk);
            });
        }
    } else {
        FArrayBox(str_xy).setVal<RunOn::Device>(0._rt, bxQstr);
    } /// get harmonic coefficient Kh at q points and harmonic part of str_xy

    if (CS.anisotropic) {
        const Real Kh_aniso = CS.Kh_aniso;
        ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            /// Horizontal-tension averaged to q-points
            /// Local variable for interpolating computed strain rates [T-1 ~> s-1].
            Real const local_strain = 0.25_rt * ( (sh_xx(i,j,kk) + sh_xx(i+1,j+1,kk)) +
                                               (sh_xx(i+1,j,kk) + sh_xx(i,j+1,kk)) );
            /// *Add* the tension contribution to the xy-component of stress
            str_xy(i,j,kk) = str_xy(i,j,kk) - Kh_aniso * n1n2_q(i,j,0) * n1n1_m_n2n2_q(i,j,0) * local_strain;
        });
    }

    if (CS.biharmonic) {
        /// Determine the biharmonic viscosity at q points, using the
        /// largest value from several parameterizations. Also get the
        /// biharmonic component of str_xy.
        ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            Ah(i,j,kk) = Ah_bg_xy(i,j,0);
        });

        if (CS.Smagorinsky_Ah || CS.Leith_Ah) {
            if (CS.Smagorinsky_Ah) {
                if (CS.bound_Coriolis) {
                    ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                    {
                        /// Smagorinsky biharmonic viscosity [L4 T-1 ~> m4 s-1]
                        Real const AhSm = Shear_mag(i,j,kk) * (Biharm_const_xy(i,j,0)
                            + Biharm_const2_xy(i,j,0) * Shear_mag(i,j,kk));
                        Ah(i,j,kk) = amrex::max(Ah(i,j,kk), AhSm);
                    });
                } else {
                    ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                    {
                        /// Smagorinsky biharmonic viscosity [L4 T-1 ~> m4 s-1]
                        Real const AhSm = Biharm_const_xy(i,j,0) * Shear_mag(i,j,kk);
                        Ah(i,j,kk) = amrex::max(Ah(i,j,kk), AhSm);
                    });
                }
            }

            if (CS.Leith_Ah) {
                ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                {
                    /// 2D Leith biharmonic viscosity [L4 T-1 ~> m4 s-1]
                    Real const AhLth = Biharm6_const_xy(i,j,0) * amrex::Math::abs(Del2vort_q(i,j,kk)) * inv_PI6;
                    Ah(i,j,kk) = amrex::max(Ah(i,j,kk), AhLth);
                });
            }
        } /// Smagorinsky_Ah or Leith_Ah

        if (use_MEKE_Au) {
            /// *Add* the MEKE contribution
            ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                Ah(i,j,kk) = Ah(i,j,kk) + 0.25_rt * (
                    (MEKE_Au(i,j,0) + MEKE_Au(i+1,j+1,0)) + (MEKE_Au(i+1,j,0) + MEKE_Au(i,j+1,0)) );
            });
        }

        /// XXX: It is just overwrites the values!
        if (CS.Re_Ah > 0.0_rt) {
            ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                Real const u_sum = u(i,j,k) + u(i,j+1,k);
                Real const v_sum = v(i,j,k) + v(i+1,j,k);
                /// Local kinetic energy [L2 T-2 ~> m2 s-2]
                Real const KE = 0.125_rt * ((u_sum*u_sum) + (v_sum*v_sum));
                Ah(i,j,kk) = std::sqrt(KE) * Re_Ah_const_xy(i,j,0);
            });
        }

        if (CS.bound_Ah) {
            if (CS.bound_Kh) {
                ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                {
                    Ah(i,j,kk) = amrex::min(Ah(i,j,kk), visc_bound_rem(i,j,kk) * hrat_min(i,j,kk) * Ah_Max_xy(i,j,0));
                });
            } else {
                ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                {
                    Ah(i,j,kk) = amrex::min(Ah(i,j,kk), hrat_min(i,j,kk) * Ah_Max_xy(i,j,0));
                });
            }
        }

        if (CS.EY24_EBT_BS) {
            const Real KS_coef = CS.KS_coef;
            ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                Real const tmp = KS_coef *hrat_min(i,j,kk) * Ah_Max_xy_KS(i,j,0);
                visc_limit_q(i,j,k) = tmp;
                visc_limit_q_frac(i,j,k) = Ah(i,j,kk) / (KS_coef * hrat_min(i,j,kk) * Ah_Max_xy_KS(i,j,0));
                if (Ah(i,j,kk) >= tmp) {
                    visc_limit_q_flag(i,j,k) = 1._rt;
                }
            });
        }

        /// Leith+E doesn't recompute Ah at q points, it just interpolates it from h to q points
        if (CS.use_Leithy) {
            ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                Ah(i,j,kk) = 0.25_rt * ((Ah_h(i,j,k) + Ah_h(i+1,j+1,k)) + (Ah_h(i,j+1,k) + Ah_h(i+1,j,k)));
            });
        }

        if (CS.do_Ah_q || CS.debug) {
            ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                Ah_q(i,j,k) = Ah(i,j,kk);
            });
        }

        /// Again, need to initialize str_xy as if its biharmonic
        ParallelFor(bxQstr, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            /// Stress tensor update [L2 T-2 ~> m2 s-2]
            Real const d_str = Ah(i,j,kk) * (dDel2vdx(i,j,kk) + dDel2udy(i,j,kk));

            str_xy(i,j,kk) = str_xy(i,j,kk) + d_str;

            /// Keep a copy of the biharmonic contribution for backscatter parameterization
            bhstr_xy(i,j,kk) = d_str * (hq(i,j,kk) * mask2dBu(i,j,0) * reduction_xy(i,j,0));
        });
    } /// Get Ah at q points and biharmonic part of str_xy

    /// Backscatter using MEKE
    if (CS.EY24_EBT_BS) {
        hor_visc_backscatter_q(bxQstr, CS, MEKE_Ku, BS_struct, use_kh_struct, kstart,
                               visc_limit_q_flag, sh_xy, str_xy, BS_coeff_q);
    } /// Backscatter

    if (CS.use_GME) {
        const Real GME_limiter = CS.GME_limiter;

        /// The wider halo here is to permit one pass of smoothing without a halo update.
        if ((CS.do_GME_coeff_h) || find_FrictWork) {
            ParallelFor(bxHgrad, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                /// The GME (negative) viscosity coefficient [L2 T-1 ~> m2 s-1]
                Real GME_coeff = GME_effic_h(i,j,0) * 0.25_rt *
                    ((KH_u_GME(i,j,k)+KH_u_GME(i-1,j,k)) + (KH_v_GME(i,j,k)+KH_v_GME(i,j-1,k)));
                GME_coeff = amrex::min(GME_coeff, GME_limiter);

                GME_coeff_h(i,j,k) = GME_coeff;
                str_xx_GME(i,j,kk) = GME_coeff * sh_xx_bt(i,j,0);
            });
        } else {
            ParallelFor(bxHgrad, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                /// The GME (negative) viscosity coefficient [L2 T-1 ~> m2 s-1]
                Real GME_coeff = GME_effic_h(i,j,0) * 0.25_rt *
                    ((KH_u_GME(i,j,k)+KH_u_GME(i-1,j,k)) + (KH_v_GME(i,j,k)+KH_v_GME(i,j-1,k)));
                GME_coeff = amrex::min(GME_coeff, GME_limiter);

                str_xx_GME(i,j,kk) = GME_coeff * sh_xx_bt(i,j,0);
            });
        }

        /// The wider halo here is to permit one pass of smoothing without a halo update.
        if (CS.do_GME_coeff_q) {
            ParallelFor(bxQshear, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                /// The GME (negative) viscosity coefficient [L2 T-1 ~> m2 s-1]
                Real GME_coeff = GME_effic_q(i,j,0) * 0.25_rt *
                    ((KH_u_GME(i,j,k)+KH_u_GME(i,j+1,k)) + (KH_v_GME(i,j,k)+KH_v_GME(i+1,j,k)));
                GME_coeff = amrex::min(GME_coeff, GME_limiter);

                GME_coeff_q(i,j,k) = GME_coeff;
                str_xy_GME(i,j,kk) = GME_coeff * sh_xy_bt(i,j,0);
            });
        } else {
            ParallelFor(bxQshear, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
            {
                const int k = kstart + kk - 1;
                /// The GME (negative) viscosity coefficient [L2 T-1 ~> m2 s-1]
                Real GME_coeff = GME_effic_q(i,j,0) * 0.25_rt *
                    ((KH_u_GME(i,j,k)+KH_u_GME(i,j+1,k)) + (KH_v_GME(i,j,k)+KH_v_GME(i+1,j,k)));
                GME_coeff = amrex::min(GME_coeff, GME_limiter);

                str_xy_GME(i,j,kk) = GME_coeff * sh_xy_bt(i,j,0);
            });
        }
    }
}

void horizontal_viscosity_tr(const Box& bxH0,
                             const Box& bxQ0,
                             const hv_opts_CS_C& CS,
                             const hv_grid_dev& G,
                             const hv_metrics_dev& M,
                             const hv_mixing_dev& X,
                             OceanOBC* OBC,
                             const int nkblock,
                             const Real H_subroundoff,
                             const Real H_to_RZ,
                             const Real s_to_T,
                             Array4<const Real> const& u,
                             Array4<const Real> const& v,
                             Array4<const Real> const& h,
                             Array4<const Real> const& uh,
                             Array4<const Real> const& vh,
                             Array4<Real>       const& diffu,
                             Array4<Real>       const& diffv,
                             Array4<const Real> const& hu_cont,
                             Array4<const Real> const& hv_cont,
                             Array4<const Real> const& u_smooth,
                             Array4<const Real> const& v_smooth,
                             Array4<const Real> const& sh_xx_bt,
                             Array4<const Real> const& sh_xy_bt,
                             Array4<const Real> const& GME_effic_h,
                             Array4<const Real> const& GME_effic_q,
                             Array4<const Real> const& KH_u_GME,
                             Array4<const Real> const& KH_v_GME,
                             Array4<Real>       const& GME_coeff_h,
                             Array4<Real>       const& GME_coeff_q,
                             Array4<Real>       const& str_xx_GME,
                             Array4<Real>       const& str_xy_GME,
                             Array4<Real>       const& NoSt,
                             Array4<Real>       const& ShSt,
                             Array4<Real>       const& FrictWork,
                             Array4<Real>       const& FrictWork_bh,
                             Array4<Real>       const& FrictWork_GME,
                             Array4<Real>       const& Ah_h,
                             Array4<Real>       const& Ah_q,
                             Array4<Real>       const& Kh_h,
                             Array4<Real>       const& Kh_q,
                             Array4<Real>       const& grid_Re_Kh,
                             Array4<Real>       const& grid_Re_Ah,
                             Array4<Real>       const& div_xx_h,
                             Array4<Real>       const& sh_xx_h,
                             Array4<Real>       const& vort_xy_q,
                             Array4<Real>       const& sh_xy_q,
                             Array4<Real>       const& visc_limit_h,
                             Array4<Real>       const& visc_limit_q,
                             Array4<Real>       const& visc_limit_h_flag,
                             Array4<Real>       const& visc_limit_q_flag,
                             Array4<Real>       const& visc_limit_h_frac,
                             Array4<Real>       const& visc_limit_q_frac,
                             Array4<Real>       const& BS_coeff_h,
                             Array4<Real>       const& BS_coeff_q)
{
    BL_PROFILE("horizontal_viscosity_tr");

    if (OBC != nullptr) {
        AMREX_ABORT_LOC("OBC pointer provided but not yet implemented");
    }

    /// True if any Leith parameterizations are enabled
    const bool use_Leith = CS.Leith_Kh || CS.Leith_Ah || CS.use_Leithy;
    /// True if vort_xy must be computed
    const bool use_vort_xy = use_Leith || CS.do_vort_xy_q || CS.use_ZB2020;

    const int is  = bxH0.smallEnd(0);
    const int ie  = bxH0.bigEnd(0);
    const int js  = bxH0.smallEnd(1);
    const int je  = bxH0.bigEnd(1);
    const int nz  = bxH0.bigEnd(2) + 1;
    const int Isq = bxQ0.smallEnd(0);
    const int Ieq = bxQ0.bigEnd(0);
    const int Jsq = bxQ0.smallEnd(1);
    const int Jeq = bxQ0.bigEnd(1);

    /// thickness so small it can be lost in roundoff and so neglected [H ~> m or kg m-2]
    const Real h_neglect  = H_subroundoff;
    /// h_neglect3 = h_neglect**3
    /// h_neglect^3 [H3 ~> m3 or kg3 m-6]
    const Real h_neglect3 = h_neglect*h_neglect*h_neglect;
    const Real pi = 4.0_rt*std::atan(1.0_rt);
    /// Powers of the inverse of pi [nondim]
    const Real inv_PI3 = 1.0_rt/(pi*pi*pi);
    /// Powers of the inverse of pi [nondim]
    const Real inv_PI6 = inv_PI3 * inv_PI3;

    if (CS.EY24_EBT_BS) {
        FArrayBox(visc_limit_h).setVal<RunOn::Device>(0._rt);
        FArrayBox(visc_limit_q).setVal<RunOn::Device>(0._rt);
        FArrayBox(visc_limit_h_flag).setVal<RunOn::Device>(0._rt);
        FArrayBox(visc_limit_q_flag).setVal<RunOn::Device>(0._rt);
        FArrayBox(visc_limit_h_frac).setVal<RunOn::Device>(0._rt);
        FArrayBox(visc_limit_q_frac).setVal<RunOn::Device>(0._rt);
    }

    const bool skeb_use_frict = X.skeb_use_frict;

    /// Fields evaluated on active layers, used for constructing 3D stress fields
    /// NOTE: The position of these declarations can impact performance, due to the
    ///   very large number of stack arrays in this function.  Move with caution!
    /// NOTE: Several of these are declared with the memory extent of q-points, but the
    ///   same arrays are also used at h-points to reduce the memory footprint of this
    ///   module, so they should never be used in halo point or checksum calls.
    /// [u_a%lb(1):u_a%ub(1), h_a%lb(2):h_a%ub(2), 1:nkblock]
    const Box bxUH_kb(IntVect(Box(u).smallEnd(0), Box(h).smallEnd(1), 0),
                       IntVect(Box(u).bigEnd(0), Box(h).bigEnd(1), nkblock-1));
    /// [h_a%lb(1):h_a%ub(1), v_a%lb(2):v_a%ub(2), 1:nkblock]
    const Box bxHV_kb(IntVect(Box(h).smallEnd(0), Box(v).smallEnd(1), 0),
                       IntVect(Box(h).bigEnd(0), Box(v).bigEnd(1), nkblock-1));
    /// [h_a%lb(1):h_a%ub(1), h_a%lb(2):h_a%ub(2), 1:nkblock]
    const Box bxHH_kb(IntVect(Box(h).smallEnd(0), Box(h).smallEnd(1), 0),
                       IntVect(Box(h).bigEnd(0), Box(h).bigEnd(1), nkblock-1));
    /// [u_a%lb(1):u_a%ub(1), v_a%lb(2):v_a%ub(2), 1:nkblock]
    const Box bxUV_kb(IntVect(Box(u).smallEnd(0), Box(v).smallEnd(1), 0),
                       IntVect(Box(u).bigEnd(0), Box(v).bigEnd(1), nkblock-1));

    /// Thickness interpolated to u points [H ~> m or kg m-2].
    FArrayBox h_u_fab(bxUH_kb, 1, The_Async_Arena());
    Array4<Real> const h_u = h_u_fab.array();
    /// The u-component of the Laplacian of velocity [L-1 T-1 ~> m-1 s-1]
    FArrayBox Del2u_fab(bxUH_kb, 1, The_Async_Arena());
    Array4<Real> const Del2u = Del2u_fab.array();
    /// Thickness interpolated to v points [H ~> m or kg m-2].
    FArrayBox h_v_fab(bxHV_kb, 1, The_Async_Arena());
    Array4<Real> const h_v = h_v_fab.array();
    /// The v-component of the Laplacian of velocity [L-1 T-1 ~> m-1 s-1]
    FArrayBox Del2v_fab(bxHV_kb, 1, The_Async_Arena());
    Array4<Real> const Del2v = Del2v_fab.array();
    /// horizontal tension from smoothed velocity including metric terms [T-1 ~> s-1]
    FArrayBox sh_xx_smooth_fab(bxHH_kb, 1, The_Async_Arena());
    Array4<Real> const sh_xx_smooth = sh_xx_smooth_fab.array();
    /// horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
    FArrayBox sh_xx_fab(bxHH_kb, 1, The_Async_Arena());
    Array4<Real> const sh_xx = sh_xx_fab.array();
    /// components in the horizontal tension [T-1 ~> s-1]
    FArrayBox dvdy_fab(bxHH_kb, 1, The_Async_Arena());
    Array4<Real> const dvdy = dvdy_fab.array();
    /// components in the horizontal tension [T-1 ~> s-1]
    FArrayBox dudx_fab(bxHH_kb, 1, The_Async_Arena());
    Array4<Real> const dudx = dudx_fab.array();
    /// Magnitude of divergence gradient at h-points [L-1 T-1 ~> m-1 s-1]
    FArrayBox grad_div_mag_h_fab(bxHH_kb, 1, The_Async_Arena());
    Array4<Real> const grad_div_mag_h = grad_div_mag_h_fab.array();
    /// Magnitude of 2d vorticity gradient at h-points [L-1 T-1 ~> m-1 s-1]
    FArrayBox grad_vort_mag_h_2d_fab(bxHH_kb, 1, The_Async_Arena());
    Array4<Real> const grad_vort_mag_h_2d = grad_vort_mag_h_2d_fab.array();
    /// Magnitude of vorticity gradient at h-points [L-1 T-1 ~> m-1 s-1]
    FArrayBox grad_vort_mag_h_fab(bxHH_kb, 1, The_Async_Arena());
    Array4<Real> const grad_vort_mag_h = grad_vort_mag_h_fab.array();
    /// A copy of str_xx that only contains the biharmonic contribution [H L2 T-2 ~> m3 s-2 or kg s-2]
    FArrayBox bhstr_xx_fab(bxHH_kb, 1, The_Async_Arena());
    Array4<Real> const bhstr_xx = bhstr_xx_fab.array();
    /// str_xx is the diagonal term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2], but at some points in the code
///     it is not yet layer integrated, so is in [L2 T-2 ~> m2 s-2].
    FArrayBox str_xx_fab(bxHH_kb, 1, The_Async_Arena());
    Array4<Real> const str_xx = str_xx_fab.array();
    /// Kh=m_leithy*Ah in Leith+E parameterization [L-2 ~> m-2]
    FArrayBox m_leithy_fab(bxHH_kb, 1, The_Async_Arena());
    Array4<Real> const m_leithy = m_leithy_fab.array();
    /// Components in the biharmonic equivalent of the shearing strain [L-2 T-1 ~> m-2 s-1]
    FArrayBox dDel2udy_fab(bxUV_kb, 1, The_Async_Arena());
    Array4<Real> const dDel2udy = dDel2udy_fab.array();
    /// Components in the biharmonic equivalent of the shearing strain [L-2 T-1 ~> m-2 s-1]
    FArrayBox dDel2vdx_fab(bxUV_kb, 1, The_Async_Arena());
    Array4<Real> const dDel2vdx = dDel2vdx_fab.array();
    /// components in the shearing strain from smoothed velocity [T-1 ~> s-1]
    FArrayBox dudy_smooth_fab(bxUV_kb, 1, The_Async_Arena());
    Array4<Real> const dudy_smooth = dudy_smooth_fab.array();
    /// components in the shearing strain from smoothed velocity [T-1 ~> s-1]
    FArrayBox dvdx_smooth_fab(bxUV_kb, 1, The_Async_Arena());
    Array4<Real> const dvdx_smooth = dvdx_smooth_fab.array();
    /// components in the shearing strain [T-1 ~> s-1]
    FArrayBox dudy_fab(bxUV_kb, 1, The_Async_Arena());
    Array4<Real> const dudy = dudy_fab.array();
    /// components in the shearing strain [T-1 ~> s-1]
    FArrayBox dvdx_fab(bxUV_kb, 1, The_Async_Arena());
    Array4<Real> const dvdx = dvdx_fab.array();
    /// horizontal shearing strain from smoothed velocity including metric terms [T-1 ~> s-1]
    FArrayBox sh_xy_smooth_fab(bxUV_kb, 1, The_Async_Arena());
    Array4<Real> const sh_xy_smooth = sh_xy_smooth_fab.array();
    /// horizontal shearing strain (du/dy + dv/dx) including metric terms [T-1 ~> s-1]
    FArrayBox sh_xy_fab(bxUV_kb, 1, The_Async_Arena());
    Array4<Real> const sh_xy = sh_xy_fab.array();
    /// harmonic mean of the harmonic means of the u- & v point thicknesses [H ~> m or kg m-2] This form guarantees that
///     hq/hu < 4.
    FArrayBox hq_fab(bxUV_kb, 1, The_Async_Arena());
    Array4<Real> const hq = hq_fab.array();
    /// Magnitude of divergence gradient at q-points [L-1 T-1 ~> m-1 s-1]
    FArrayBox grad_div_mag_q_fab(bxUV_kb, 1, The_Async_Arena());
    Array4<Real> const grad_div_mag_q = grad_div_mag_q_fab.array();
    /// Laplacian of vorticity at q-points [L-2 T-1 ~> m-2 s-1]
    FArrayBox Del2vort_q_fab(bxUV_kb, 1, The_Async_Arena());
    Array4<Real> const Del2vort_q = Del2vort_q_fab.array();
    /// Magnitude of 2d vorticity gradient at q-points [L-1 T-1 ~> m-1 s-1]
    FArrayBox grad_vort_mag_q_2d_fab(bxUV_kb, 1, The_Async_Arena());
    Array4<Real> const grad_vort_mag_q_2d = grad_vort_mag_q_2d_fab.array();
    /// Magnitude of vorticity gradient at q-points [L-1 T-1 ~> m-1 s-1]
    FArrayBox grad_vort_mag_q_fab(bxUV_kb, 1, The_Async_Arena());
    Array4<Real> const grad_vort_mag_q = grad_vort_mag_q_fab.array();
    /// Vertical vorticity including metric terms, smoothed [T-1 ~> s-1]
    FArrayBox vort_xy_smooth_fab(bxUV_kb, 1, The_Async_Arena());
    Array4<Real> const vort_xy_smooth = vort_xy_smooth_fab.array();
    /// Vertical vorticity (dv/dx - du/dy) including metric terms [T-1 ~> s-1]
    FArrayBox vort_xy_fab(bxUV_kb, 1, The_Async_Arena());
    Array4<Real> const vort_xy = vort_xy_fab.array();
    /// A copy of str_xy that only contains the biharmonic contribution [H L2 T-2 ~> m3 s-2 or kg s-2]
    FArrayBox bhstr_xy_fab(bxUV_kb, 1, The_Async_Arena());
    Array4<Real> const bhstr_xy = bhstr_xy_fab.array();
    /// str_xy is the cross term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2], but at some points in the code it
///     is not yet layer integrated, so is in [L2 T-2 ~> m2 s-2].
    FArrayBox str_xy_fab(bxUV_kb, 1, The_Async_Arena());
    Array4<Real> const str_xy = str_xy_fab.array();
    /// fraction of overall viscous bounds that remain to be applied (h or q) [nondim]
    FArrayBox visc_bound_rem_fab(bxUV_kb, 1, The_Async_Arena());
    Array4<Real> const visc_bound_rem = visc_bound_rem_fab.array();
    /// h_min divided by the thickness at the stress point (h or q) [nondim]
    FArrayBox hrat_min_fab(bxUV_kb, 1, The_Async_Arena());
    Array4<Real> const hrat_min = hrat_min_fab.array();
    /// magnitude of gradient of smoothed vertical vorticity (h or q) [L-1 T-1 ~> m-1 s-1]
    FArrayBox vert_vort_mag_smooth_fab(bxUV_kb, 1, The_Async_Arena());
    Array4<Real> const vert_vort_mag_smooth = vert_vort_mag_smooth_fab.array();
    /// magnitude of the vertical vorticity gradient (h or q) [L-1 T-1 ~> m-1 s-1]
    FArrayBox vert_vort_mag_fab(bxUV_kb, 1, The_Async_Arena());
    Array4<Real> const vert_vort_mag = vert_vort_mag_fab.array();
    /// magnitude of the shear (h or q) [T-1 ~> s-1]
    FArrayBox Shear_mag_fab(bxUV_kb, 1, The_Async_Arena());
    Array4<Real> const Shear_mag = Shear_mag_fab.array();
    /// Laplacian  viscosity (h or q) [L2 T-1 ~> m2 s-1]
    FArrayBox Kh_fab(bxUV_kb, 1, The_Async_Arena());
    Array4<Real> const Kh = Kh_fab.array();
    /// biharmonic viscosity (h or q) [L4 T-1 ~> m4 s-1]
    FArrayBox Ah_fab(bxUV_kb, 1, The_Async_Arena());
    Array4<Real> const Ah = Ah_fab.array();
    /// y-derivative of horizontal divergence (d/dy(du/dx + dv/dy)) [L-1 T-1 ~> m-1 s-1]
    FArrayBox div_xx_dy_fab(bxHV_kb, 1, The_Async_Arena());
    Array4<Real> const div_xx_dy = div_xx_dy_fab.array();
    /// x-derivative of smoothed vertical vorticity [L-1 T-1 ~> m-1 s-1]
    FArrayBox vort_xy_dx_smooth_fab(bxHV_kb, 1, The_Async_Arena());
    Array4<Real> const vort_xy_dx_smooth = vort_xy_dx_smooth_fab.array();
    /// x-derivative of vertical vorticity (d/dx(dv/dx - du/dy)) [L-1 T-1 ~> m-1 s-1]
    FArrayBox vort_xy_dx_fab(bxHV_kb, 1, The_Async_Arena());
    Array4<Real> const vort_xy_dx = vort_xy_dx_fab.array();
    /// x-derivative of horizontal divergence (d/dx(du/dx + dv/dy)) [L-1 T-1 ~> m-1 s-1]
    FArrayBox div_xx_dx_fab(bxUH_kb, 1, The_Async_Arena());
    Array4<Real> const div_xx_dx = div_xx_dx_fab.array();
    /// y-derivative of smoothed vertical vorticity [L-1 T-1 ~> m-1 s-1]
    FArrayBox vort_xy_dy_smooth_fab(bxUH_kb, 1, The_Async_Arena());
    Array4<Real> const vort_xy_dy_smooth = vort_xy_dy_smooth_fab.array();
    /// y-derivative of vertical vorticity (d/dy(dv/dx - du/dy)) [L-1 T-1 ~> m-1 s-1]
    FArrayBox vort_xy_dy_fab(bxUH_kb, 1, The_Async_Arena());
    Array4<Real> const vort_xy_dy = vort_xy_dy_fab.array();
    /// The square of the biharmonic viscosity [L8 T-2 ~> m8 s-2]
    FArrayBox Ah_sq_fab(bxHH_kb, 1, The_Async_Arena());
    Array4<Real> const Ah_sq = Ah_sq_fab.array();

    if (CS.use_Leithy) {
        /// Initialize
        FArrayBox(m_leithy).setVal<RunOn::Device>(0.0_rt);
    }

    bool apply_OBC = false;
    /*
    if (OBC != nullptr) {
        if (OBC->OBC_pe) {
            apply_OBC = OBC->Flather_u_BCs_exist_globally || OBC->Flather_v_BCs_exist_globally;
            apply_OBC = true;
        }
    }
    */

    bool apply_OBC_strain = false;
    /*
    if (OBC != nullptr) {
        apply_OBC_strain = (OBC->strain_config != OBC_STRAIN_NONE)
          && ((!CS.OBC_strain_bug) || (OBC->strain_config != OBC_STRAIN_SPECIFIED));
    }
    */


    const bool find_FrictWork = CS.do_FrictWork || CS.do_FrictWorkIntz
        || (X.MEKE_mom_src.p != nullptr);
    const bool find_FrictWork_bh = CS.do_FrictWork_bh || CS.do_FrictWorkIntz_bh
        || (X.MEKE_mom_src_bh.p != nullptr);

    const bool use_kh_struct = (X.BS_struct.p != nullptr);
    /// The ratio of f over Shear_mag that is so small that the backscatter
    /// calculation gives the same value as if f were 0 [nondim].
    Real backscat_subround = 0.0_rt;
    if (find_FrictWork && (X.MEKE_mom_src.p != nullptr) && (X.backscatter_Ro_c > 0.0_rt) &&
        (X.backscatter_Ro_pow != 0.0_rt)) {
        backscat_subround = std::pow(1.0e-16_rt/X.backscatter_Ro_c, 1.0_rt/X.backscatter_Ro_pow);
    }

    /// Toggle whether to use a Laplacian viscosity derived from MEKE
    const bool use_MEKE_Ku = (X.MEKE_Ku.p != nullptr);
    const bool use_MEKE_Au = (X.MEKE_Au.p != nullptr);

    const bool use_cont_huv = CS.use_cont_thick && (hu_cont.p != nullptr) && (hv_cont.p != nullptr);

    bool rescale_Kh = false;
    if (X.use_variable_mixing) {
        rescale_Kh = X.Resoln_scaled_Kh;
    }

    /// The number of points below the h points, in i and j, over which the
    /// thickness-point viscosities are computed [nondim]
    int kh_lo[2];
    /// Set the halo sizes used for the thickness-point viscosities.
    /// They always extend one point above the h points, for the normal stresses over bxHstr.
    if (CS.use_Leithy || CS.debug) {
        /// One point below: Leith+E interpolates the h-point Ah to q points from isc-1,
        /// and the debug checksums cover it
        kh_lo[0] = 1;
        kh_lo[1] = 1;
    } else {
        /// Down to the first u (v) point: 1 in symmetric memory, else 0
        kh_lo[0] = is-Isq;
        kh_lo[1] = js-Jsq;
    }

    Box bxH0_2d = bxH0;
    bxH0_2d.setSmall(2, 0);
    bxH0_2d.setBig(2, 0);
    if (find_FrictWork && (X.MEKE_mom_src.p != nullptr)) {
        FArrayBox(X.MEKE_mom_src).setVal<RunOn::Device>(0._rt, bxH0_2d);

        if (X.MEKE_mom_src_bh.p != nullptr) {
            FArrayBox(X.MEKE_mom_src_bh).setVal<RunOn::Device>(0._rt, bxH0_2d);
        }

        if (X.MEKE_GME_snk.p != nullptr) {
            FArrayBox(X.MEKE_GME_snk).setVal<RunOn::Device>(0._rt, bxH0_2d);
        }
    }

    /// TODO: Only needed if FrictWork_bh is true, and currently only used on CPU,
    ///   but I do not yet see any benefit to breaking up the calculation.

    for (int kstart = 1; kstart <= nz; kstart += nkblock) {
        const int kend = amrex::min(kstart+nkblock-1, nz);
        const int kmax = kend-kstart+1;
        /// The h-point iteration box of one tile, [isc:iec, jsc:jec, 1:kmax]
        const Box bxH(IntVect(is, js, 0), IntVect(ie, je, kmax-1));
        /// The B-grid-index iteration box of one tile, [IscB:IecB, JscB:JecB, 1:kmax]
        const Box bxQ(IntVect(Isq, Jsq, 0), IntVect(Ieq, Jeq, kmax-1));
        /// The u-point iteration box of one tile, [IscB:IecB, jsc:jec, 1:kmax]
        const Box bxU(IntVect(Isq, js, 0), IntVect(Ieq, je, kmax-1));
        /// The v-point iteration box of one tile, [isc:iec, JscB:JecB, 1:kmax]
        const Box bxV(IntVect(is, Jsq, 0), IntVect(ie, Jeq, kmax-1));
        /// The thickness-point viscosity box of one tile, [is_Kh:ie_Kh, js_Kh:je_Kh, 1:kmax]
        const Box bxKh = grow_by(bxH, kh_lo[0], kh_lo[1], 1, 1);
        /// The thickness-point viscosity range in i only, [is_Kh:ie_Kh, jsc:jec, 1:kmax]
        const Box bxKh_i = grow_by(bxH, kh_lo[0], 0, 1, 0);
        /// The thickness-point viscosity range in j only, [isc:iec, js_Kh:je_Kh, 1:kmax]
        const Box bxKh_j = grow_by(bxH, 0, kh_lo[1], 0, 1);
        /// h points of one tile where the normal stresses are computed,
        /// [IscB:IecB+1, JscB:JecB+1, 1:kmax]: the u- and v-point divergence needs
        /// them on both sides of each velocity point; built from bxQ because the
        /// range starts at IscB. Used by strain, h_coeffs_a/b, backscatter_h, stress_div
        const Box bxHstr = grow_by(bxQ, 0, 0, 1, 1);
        /// q points of one tile where the shear stresses and the q-point viscosities
        /// are computed, [isc-1:iec, jsc-1:jec, 1:kmax]: the u- and v-point divergence
        /// needs them on both sides of each velocity point; built from bxH because the
        /// range starts at isc-1 in both memory modes. Used by strain, Leith_grad_a/b,
        /// h_coeffs_b, q_coeffs, backscatter_q, stress_div
        const Box bxQstr = grow_by(bxH, 1, 1, 0, 0);
        /// h points of one tile where the velocity gradients dudx, dvdy and sh_xx (and
        /// str_xx_GME) are computed, [IscB-1:IecB+2, JscB-1:JecB+2, 1:kmax]: Del2u and
        /// Del2v, one point beyond the velocity points, read them on both sides; built
        /// from bxQ because the range is set in q indices. Used by strain, q_coeffs
        const Box bxHgrad = grow_by(bxQ, 1, 1, 2, 2);
        /// q points of one tile where the shearing strain sh_xy (and str_xy_GME) is
        /// computed, [isc-2:iec+1, jsc-2:jec+1, 1:kmax]: Del2u and Del2v, one point
        /// beyond the velocity points, read it on both sides; built from bxH because
        /// the range starts at isc-2 in both memory modes. Used by strain, q_coeffs
        const Box bxQshear = grow_by(bxH, 2, 2, 1, 1);
        /// q points of one tile where the smoothed vorticity terms and Del2vort_q
        /// are computed, [is_Kh-1:ie_Kh, js_Kh-1:je_Kh, 1:kmax]: the h-point
        /// quantities over bxKh average the four surrounding q points. Used by
        /// strain, Leith_grad_a
        const Box bxQvsmooth = grow_by(bxKh, 1, 1, 0, 0);

        hor_visc_strain(bxH, bxQ, bxU, bxV, bxKh, bxHstr, bxQstr, bxHgrad, bxQshear, bxQvsmooth, CS, G.IdyCu, G.IdxCv,
                        G.IdyCv, G.IdxCu, G.mask2dT, G.dxBu, G.mask2dBu, G.IareaBu, G.dyCv, G.dxCu, M.dx2h, M.dy2h,
                        M.dx_dyT, M.dy_dxT, M.dx2q, M.dy2q, M.dx_dyBu, M.dy_dxBu, M.Idx2dyCu, M.Idxdy2u, M.Idx2dyCv,
                        M.Idxdy2v, OBC, kstart, use_vort_xy, apply_OBC, apply_OBC_strain, use_cont_huv, u, v, h,
                        hu_cont, hv_cont, u_smooth, v_smooth, NoSt, ShSt, Del2u, h_u, Del2v, h_v, sh_xx, sh_xx_smooth,
                        dudx, dvdy, dvdx, dudy, dvdx_smooth, dudy_smooth, sh_xy, sh_xy_smooth, vort_xy, vort_xy_smooth);

        if (use_Leith) {
            hor_visc_leith_grad_a(bxH, bxKh, bxKh_i, bxKh_j, bxQstr, bxQvsmooth, CS, G.dyBu, G.IdxBu, G.IdyCu, G.dxBu,
                                  G.IdyBu, G.IdxCv, G.IdyCv, G.IdxCu, G.df_dx, G.df_dy, dudx, dvdy, vort_xy,
                                  vort_xy_smooth, grad_vort_mag_h_2d, grad_div_mag_h, grad_vort_mag_q_2d,
                                  grad_div_mag_q, Del2vort_q, vort_xy_dx, vort_xy_dx_smooth, div_xx_dy, vort_xy_dy,
                                  vort_xy_dy_smooth, div_xx_dx);
            if (CS.use_QG_Leith_visc) {
                AMREX_ABORT_LOC("horizontal_viscosity_tr: QG Leith viscosity not yet implemented");
                /// call hor_visc_QG_Leith_host(host, nkblock, kstart, kmax, h_a, dz_a, slope_x_a, slope_y_a, &
                ///                             div_xx_dx_a, div_xx_dy_a, vort_xy_dx_a, vort_xy_dy_a)
            }
            hor_visc_leith_grad_b(bxKh, bxQstr, CS, grad_vort_mag_h, grad_vort_mag_q, vert_vort_mag_smooth, vort_xy_dx,
                                  vort_xy_dx_smooth, vort_xy_dy, vort_xy_dy_smooth);
        } /// CS%Leith_Kh

        hor_visc_h_coeffs_a(bxH, bxKh, bxHstr, CS, X.MEKE_Ku, X.Res_fn_h, X.BS_struct, M.Kh_bg_xx, M.Ah_bg_xx,
                            M.Kh_Max_xx, M.n1n2_h, M.n1n1_m_n2n2_h, M.grid_sp_h2, M.Laplac2_const_xx,
                            M.Biharm6_const_xx, M.Laplac3_const_xx, M.Biharm_const_xx, M.Biharm_const2_xx, kstart,
                            h_neglect, inv_PI3, inv_PI6, use_kh_struct, use_MEKE_Ku, rescale_Kh, u, v, h, Kh_h,
                            grid_Re_Kh, div_xx_h, sh_xx_h, h_u, h_v, sh_xx, str_xx, grad_vort_mag_h, grad_vort_mag_h_2d,
                            grad_div_mag_h, dudx, dvdy, sh_xy, Del2vort_q, Ah, Kh, Shear_mag, vert_vort_mag, hrat_min,
                            visc_bound_rem);

        if (CS.biharmonic) {
            if ((CS.Smagorinsky_Ah) || (CS.Leith_Ah) || (CS.use_Leithy)) {
                if (CS.use_Leithy) {
                    hor_visc_leithy_m(bxKh, CS, G.mask2dBu, M.Ah_bg_xx, M.m_const_leithy, M.m_leithy_max,
                                      M.Biharm6_const_xx, inv_PI6, Del2vort_q, vert_vort_mag, vort_xy_smooth,
                                      m_leithy);
                    if (CS.smooth_Ah) {
                        AMREX_ABORT_LOC("horizontal_viscosity_tr: m_leithy halo update not yet implemented");
                        /// One layer of the field being smoothed [arbitrary]
                        Box layer_bx = Box(m_leithy);
                        layer_bx.setSmall(2, 0);
                        layer_bx.setBig(2, 0);
                        FArrayBox layer_fab(layer_bx, 1, The_Async_Arena());
                        Array4<Real> const layer = layer_fab.array();
                        for (int kk = 0; kk < kmax; ++kk) {
                            /// call pass_var(m_leithy(:,:,kk), G%Domain, halo=2)
                            ParallelFor(layer_bx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
                            {
                                layer(i,j,0) = m_leithy(i,j,kk);
                            });
                            smooth_x9_h(bxH0, G.mask2dT, layer, true);
                            ParallelFor(layer_bx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
                            {
                                m_leithy(i,j,kk) = layer(i,j,0);
                            });
                            /// call pass_var(m_leithy(:,:,kk), G%Domain)
                        }
                    }
                    hor_visc_leithy_ah(bxKh, CS, M.Ah_bg_xx, M.Biharm6_const_xx, inv_PI6, Del2vort_q,
                                       vert_vort_mag_smooth, Ah, m_leithy, Ah_sq);
                    if (CS.smooth_Ah) {
                        AMREX_ABORT_LOC("horizontal_viscosity_tr: Ah_sq halo update not yet implemented");
                        /// One layer of the field being smoothed [arbitrary]
                        Box layer_bx = Box(Ah_sq);
                        layer_bx.setSmall(2, 0);
                        layer_bx.setBig(2, 0);
                        FArrayBox layer_fab(layer_bx, 1, The_Async_Arena());
                        Array4<Real> const layer = layer_fab.array();
                        for (int kk = 0; kk < kmax; ++kk) {
                            /// call pass_var(Ah_sq(:,:,kk), G%Domain, halo=2)
                            ParallelFor(layer_bx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
                            {
                                layer(i,j,0) = Ah_sq(i,j,kk);
                            });
                            smooth_x9_h(bxH0, G.mask2dT, layer, false);
                            ParallelFor(layer_bx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
                            {
                                Ah_sq(i,j,kk) = layer(i,j,0);
                            });
                            /// call pass_var(Ah_sq(:,:,kk), G%Domain)
                        }
                    }
                    hor_visc_leithy_ah_h(bxKh, CS, M.Ah_bg_xx, kstart, Ah, Ah_h, Ah_sq);
                }
            } /// Smagorinsky_Ah or Leith_Ah or Leith+E
        }

        hor_visc_h_coeffs_b(bxH, bxQ, bxKh, bxHstr, bxQstr, CS, X.MEKE_Au, X.MEKE_Ku, X.BS_struct, G.IdyCu, G.IdxCv,
                            G.IdyCv, G.IdxCu, G.mask2dBu, G.mask2dCu, G.mask2dCv, M.reduction_xx, M.Ah_Max_xx,
                            M.Ah_Max_xx_KS, M.grid_sp_h3, M.dx_dyT, M.dy_dxT, M.dx_dyBu, M.dy_dxBu, M.Re_Ah_const_xx,
                            OBC, kstart, h_neglect, h_neglect3, apply_OBC, use_kh_struct, use_MEKE_Au, u, v, h, Ah_h,
                            Kh_h, grid_Re_Ah, visc_limit_h, visc_limit_h_flag, visc_limit_h_frac, BS_coeff_h, Del2u,
                            h_u, Del2v, h_v, sh_xx, sh_xx_smooth, str_xx, bhstr_xx, m_leithy, dDel2vdx, dDel2udy, sh_xy,
                            hq, Ah, Kh, Shear_mag, hrat_min, visc_bound_rem);

        hor_visc_q_coeffs(bxQstr, bxHgrad, bxQshear, CS, X.MEKE_Ku, X.MEKE_Au, X.Res_fn_q, X.BS_struct, G.mask2dBu,
                          M.Kh_bg_xy, M.Ah_bg_xy, M.reduction_xy, M.Kh_Max_xy, M.Ah_Max_xy, M.Ah_Max_xy_KS, M.n1n2_q,
                          M.n1n1_m_n2n2_q, M.Laplac2_const_xy, M.Biharm6_const_xy, M.Laplac3_const_xy,
                          M.Biharm_const_xy, M.Biharm_const2_xy, M.Re_Ah_const_xy, kstart, inv_PI3, inv_PI6,
                          find_FrictWork, use_kh_struct, use_MEKE_Ku, use_MEKE_Au, rescale_Kh, u, v, sh_xx_bt, sh_xy_bt,
                          GME_effic_h, GME_effic_q, KH_u_GME, KH_v_GME, GME_coeff_h, GME_coeff_q, str_xx_GME,
                          str_xy_GME, Ah_h, Ah_q, Kh_h, Kh_q, vort_xy_q, sh_xy_q, visc_limit_q, visc_limit_q_flag,
                          visc_limit_q_frac, BS_coeff_q, sh_xx, dDel2vdx, dDel2udy, sh_xy, sh_xy_smooth, str_xy,
                          bhstr_xy, vort_xy, grad_vort_mag_q, grad_vort_mag_q_2d, Del2vort_q, grad_div_mag_q, hq, Ah,
                          Kh, Shear_mag, vert_vort_mag, hrat_min, visc_bound_rem);

        if (CS.use_GME) {
            AMREX_ABORT_LOC("horizontal_viscosity_tr: GME smoothing not yet implemented");
            /// Applying GME diagonal term.  This is linear and the arguments can be rescaled.
            /// smooth_GME is an impure host routine, so str_xx_GME/str_xy_GME must be brought back
            /// from the device before this loop and pushed back afterward.
            /// call hor_visc_GME_smooth_host(host, Gcore, bxH0, bxQ0, nkblock, kmax, str_xx_GME_a, &
            ///                               str_xy_GME_a)
        }

        hor_visc_stress_div(bxH, bxU, bxV, bxHstr, bxQstr, CS, X.MEKE_mom_src, X.MEKE_mom_src_bh, X.MEKE_GME_snk,
                            X.backscatter_Ro_pow, X.backscatter_Ro_c, X.skeb_diss, X.skeb_frict_coef, G.mask2dBu,
                            G.IdxCu, G.IdyCu, G.IareaCu, G.IdyCv, G.IdxCv, G.IareaCv, G.IdxT, G.IdyT, G.IdyBu, G.IdxBu,
                            G.IareaT, G.dxCu, G.dyCv, G.CoriolisBu, M.reduction_xx, M.reduction_xy, M.dx2h, M.dy2h,
                            M.dx2q, M.dy2q, OBC, kstart, h_neglect, apply_OBC, find_FrictWork, find_FrictWork_bh,
                            backscat_subround, skeb_use_frict, u, v, h, uh, vh, diffu, diffv, str_xx_GME, str_xy_GME,
                            FrictWork, FrictWork_bh, FrictWork_GME, visc_limit_h_flag, h_u, h_v, sh_xx, str_xx,
                            bhstr_xx, sh_xy, str_xy, bhstr_xy, hq, H_to_RZ, s_to_T);

        /// Pass the velocity gradients and thickness to ZB2020
        if (CS.use_ZB2020) {
            AMREX_ABORT_LOC("horizontal_viscosity_tr: ZB2020 gradient copy not yet implemented");
            /// call hor_visc_ZB2020_copy_host(host, nkblock, kstart, kmax, sh_xx_a, sh_xy_a, vort_xy_a, hq_a)
        }
    } /// end of kstart block loop
}

} // namespace MOM
