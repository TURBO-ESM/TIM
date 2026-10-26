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

using namespace amrex;

namespace MOM {

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
        if (mask(i,j,0) > 0.0) {
            Real Iwts = 0.0625;  ///< The inverse of the sum of the weights [nondim]
            if constexpr (!zero_land) {
                Iwts = 1.0 / ( (4.0*mask(i,j,0) +
                                ( 2.0*((mask(i-1,j,0) + mask(i+1,j,0)) +
                                       (mask(i,j-1,0) + mask(i,j+1,0))) +
                                 ((mask(i-1,j-1,0) + mask(i+1,j+1,0)) +
                                  (mask(i-1,j+1,0) + mask(i+1,j-1,0))) ) ) + 1.0e-16 );
            }
            field(i,j,0) = Iwts * ( 4.0*mask(i,j,0) * f_prev(i,j,0)
                                  + (2.0*((mask(i-1,j,0) * f_prev(i-1,j,0) + mask(i+1,j,0) * f_prev(i+1,j,0)) +
                                          (mask(i,j-1,0) * f_prev(i,j-1,0) + mask(i,j+1,0) * f_prev(i,j+1,0)))
                                    + ((mask(i-1,j-1,0) * f_prev(i-1,j-1,0) + mask(i+1,j+1,0) * f_prev(i+1,j+1,0)) +
                                       (mask(i-1,j+1,0) * f_prev(i-1,j+1,0) + mask(i+1,j-1,0) * f_prev(i-1,j-1,0))) ));
        }
    });
}

///
/// @brief Copy a field into the buffer holding its previous value.
///
/// @param bx     The box covering the whole field
/// @param field  The field being filtered [arbitrary]
/// @param f_prev The previous value of the field [arbitrary]
///
static void copy_prev_field(const Box& bx,
                            Array4<const Real> const& field,
                            Array4<Real>       const& f_prev)
{
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        f_prev(i,j,k) = field(i,j,k);
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
        /// bxH0 grown by s cells for pass s, [is-s:ie+s, js-s:je+s]
        const Box bxH0_gs = amrex::grow(bxH0_2d, IntVect(s, s, 0));

        copy_prev_field(fh_bx, field_h, fh_prev);
        /// apply smoothing on field_h using rotationally symmetric expressions.
        if (zero_land) {
            smooth_x9_pass<true>(bxH0_gs, mask2dT, fh_prev, field_h);
        } else {
            smooth_x9_pass<false>(bxH0_gs, mask2dT, fh_prev, field_h);
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
        /// bxU0 grown by s cells for pass s, [Isq-s:Ieq+s, js-s:je+s]
        const Box bxU0_gs = amrex::grow(bxU0_2d, IntVect(s, s, 0));
        /// bxV0 grown by s cells for pass s, [is-s:ie+s, Jsq-s:Jeq+s]
        const Box bxV0_gs = amrex::grow(bxV0_2d, IntVect(s, s, 0));

        copy_prev_field(fu_bx, field_u, fu_prev);
        /// apply smoothing on field_u using the original non-rotationally symmetric expressions.
        if (zero_land) {
            smooth_x9_pass<true>(bxU0_gs, mask2dCu, fu_prev, field_u);
        } else {
            smooth_x9_pass<false>(bxU0_gs, mask2dCu, fu_prev, field_u);
        }

        copy_prev_field(fv_bx, field_v, fv_prev);
        /// apply smoothing on field_v using the original non-rotationally symmetric expressions.
        if (zero_land) {
            smooth_x9_pass<true>(bxV0_gs, mask2dCv, fv_prev, field_v);
        } else {
            smooth_x9_pass<false>(bxV0_gs, mask2dCv, fv_prev, field_v);
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
        if (mask(i,j,0) == 0.0) return;
        /// compute weights
        Real const ww = 0.125 * mask(i-1,j,0);
        Real const we = 0.125 * mask(i+1,j,0);
        Real const ws = 0.125 * mask(i,j-1,0);
        Real const wn = 0.125 * mask(i,j+1,0);
        Real const wc = 1.0 - ((ww+we)+(wn+ws));
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

    copy_prev_field(orig_bx, GME_flux_h, GME_flux_h_original);
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

    copy_prev_field(orig_bx, GME_flux_q, GME_flux_q_original);
    smooth_gme_pass(bxPass_2d, mask2dBu, GME_flux_q_original, GME_flux_q);
}

///
/// @brief Adds the MEKE-based (EY24_EBT_BS) backscatter contribution to the diagonal term of
///        the stress tensor at h-points.
///
/// @tparam use_kh_struct    If true, shape the backscatter coefficient with VarMix%BS_struct
/// @tparam do_BS_coeff_h    If true, the BS_coeff_h diagnostic is requested.
/// @param bxQ_pij           bxQ grown to [Isq:Ieq+1, Jsq:Jeq+1]
/// @param MEKE_Ku           The MEKE-derived lateral viscosity coefficient [L2 T-1 ~> m2 s-1]
/// @param BS_struct         Vertical structure function used in backscatter [nondim]
/// @param kstart            The first absolute k-layer of the current k-block
/// @param visc_limit_h_flag determines whether backscatter is shut off [nondim]
/// @param sh_xx             horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
/// @param str_xx            The diagonal term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2]
/// @param BS_coeff_h        A diagnostic array of the backscatter coefficient [L2 T-1 ~> m2 s-1]
///
template <bool use_kh_struct, bool do_BS_coeff_h>
void hor_visc_backscatter_h_loop(const Box& bxQ_pij,
                                 Array4<const Real> const& MEKE_Ku,
                                 Array4<const Real> const& BS_struct,
                                 const int kstart,
                                 Array4<const Real> const& visc_limit_h_flag,
                                 Array4<const Real> const& sh_xx,
                                 Array4<Real>       const& str_xx,
                                 Array4<Real>       const& BS_coeff_h)
{
    ParallelFor(bxQ_pij, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
    {
        const int k = kstart + kk - 1;

        Real BS_struct_k = 0.0;
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

void hor_visc_backscatter_h(const Box& bxQ,
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

    /// bxQ grown to [Isq:Ieq+1, Jsq:Jeq+1]
    Box bxQ_pij = bxQ;
    bxQ_pij.growHi(0, 1);
    bxQ_pij.growHi(1, 1);

    if (use_kh_struct) {
        if (CS.do_BS_coeff_h) {
            hor_visc_backscatter_h_loop<true, true>(bxQ_pij, MEKE_Ku, BS_struct, kstart,
                                                    visc_limit_h_flag, sh_xx, str_xx, BS_coeff_h);
        } else {
            hor_visc_backscatter_h_loop<true, false>(bxQ_pij, MEKE_Ku, BS_struct, kstart,
                                                     visc_limit_h_flag, sh_xx, str_xx, BS_coeff_h);
        }
    } else {
        if (CS.do_BS_coeff_h) {
            hor_visc_backscatter_h_loop<false, true>(bxQ_pij, MEKE_Ku, BS_struct, kstart,
                                                     visc_limit_h_flag, sh_xx, str_xx, BS_coeff_h);
        } else {
            hor_visc_backscatter_h_loop<false, false>(bxQ_pij, MEKE_Ku, BS_struct, kstart,
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
/// @param bxH_mij           bxH grown to [is-1:ie, js-1:je]
/// @param MEKE_Ku           The MEKE-derived lateral viscosity coefficient [L2 T-1 ~> m2 s-1]
/// @param BS_struct         Vertical structure function used in backscatter [nondim]
/// @param kstart            The first absolute k-layer of the current k-block
/// @param visc_limit_q_flag determines whether backscatter is shut off [nondim]
/// @param sh_xy             horizontal shearing strain (du/dy + dv/dx) including metric terms [T-1 ~> s-1]
/// @param str_xy            The cross term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2]
/// @param BS_coeff_q        A diagnostic array of the backscatter coefficient [L2 T-1 ~> m2 s-1]
///
template <bool use_kh_struct, bool do_BS_coeff_q>
void hor_visc_backscatter_q_loop(const Box& bxH_mij,
                                 Array4<const Real> const& MEKE_Ku,
                                 Array4<const Real> const& BS_struct,
                                 const int kstart,
                                 Array4<const Real> const& visc_limit_q_flag,
                                 Array4<const Real> const& sh_xy,
                                 Array4<Real>       const& str_xy,
                                 Array4<Real>       const& BS_coeff_q)
{
    ParallelFor(bxH_mij, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
    {
        const int k = kstart + kk - 1;

        Real Kh_BS;  ///< Laplacian antiviscosity [L2 T-1 ~> m2 s-1]
        if (visc_limit_q_flag(i,j,k) > 0.0) {
            Kh_BS = 0.;
        } else {
            if constexpr (use_kh_struct) {
                Kh_BS = 0.25*( ((MEKE_Ku(i,j,0)*BS_struct(i,j,k)) +
                                (MEKE_Ku(i+1,j+1,0)*BS_struct(i+1,j+1,k))) +
                               ((MEKE_Ku(i+1,j,0)*BS_struct(i+1,j,k)) +
                                (MEKE_Ku(i,j+1,0)*BS_struct(i,j+1,k))) );
            } else {
                Kh_BS = 0.25*( (MEKE_Ku(i,j,0) + MEKE_Ku(i+1,j+1,0)) +
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

void hor_visc_backscatter_q(const Box& bxH,
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

    /// bxH grown to [is-1:ie, js-1:je]
    Box bxH_mij = bxH;
    bxH_mij.growLo(0, 1);
    bxH_mij.growLo(1, 1);

    if (use_kh_struct) {
        if (CS.do_BS_coeff_q) {
            hor_visc_backscatter_q_loop<true, true>(bxH_mij, MEKE_Ku, BS_struct, kstart,
                                                    visc_limit_q_flag, sh_xy, str_xy, BS_coeff_q);
        } else {
            hor_visc_backscatter_q_loop<true, false>(bxH_mij, MEKE_Ku, BS_struct, kstart,
                                                     visc_limit_q_flag, sh_xy, str_xy, BS_coeff_q);
        }
    } else {
        if (CS.do_BS_coeff_q) {
            hor_visc_backscatter_q_loop<false, true>(bxH_mij, MEKE_Ku, BS_struct, kstart,
                                                     visc_limit_q_flag, sh_xy, str_xy, BS_coeff_q);
        } else {
            hor_visc_backscatter_q_loop<false, false>(bxH_mij, MEKE_Ku, BS_struct, kstart,
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

    /// bxQ0 grown to [Isq-2:Ieq+2, Jsq-2:Jeq+2]
    const Box bxQ0_g2ij = amrex::grow(bxQ0_2d, IntVect(2, 2, 0));
    /// bxQ0 grown to [Isq-1:Ieq+2, Jsq-1:Jeq+2]
    Box bxQ0_mij_p2ij = bxQ0_2d;
    bxQ0_mij_p2ij.growLo(0, 1).growLo(1, 1).growHi(0, 2).growHi(1, 2);
    /// bxH0 grown to [is-2:ie+2, js-2:je+2]
    const Box bxH0_g2ij = amrex::grow(bxH0_2d, IntVect(2, 2, 0));
    /// bxH0 grown to [is-2:ie+1, js-2:je+1]
    Box bxH0_m2ij_pij = bxH0_2d;
    bxH0_m2ij_pij.growLo(0, 2).growLo(1, 2).growHi(0, 1).growHi(1, 1);

    /// The total thickness of all layers [H ~> m or kg m-2]
    Box htot_bx = Box(h);
    htot_bx.setSmall(2, 0);
    htot_bx.setBig(2, 0);
    FArrayBox htot_fab(htot_bx, 1, The_Async_Arena());
    Array4<Real> const htot = htot_fab.array();

    /// Calculate the barotropic horizontal tension
    ParallelFor(bxH0_g2ij, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        dudx_bt(i,j,0) = dy_dxT(i,j,0)*((IdyCu(i,j,0) * ubtav(i,j,0)) -
                                        (IdyCu(i-1,j,0) * ubtav(i-1,j,0)));
        dvdy_bt(i,j,0) = dx_dyT(i,j,0)*((IdxCv(i,j,0) * vbtav(i,j,0)) -
                                        (IdxCv(i,j-1,0) * vbtav(i,j-1,0)));
    });
    ParallelFor(bxQ0_mij_p2ij, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        sh_xx_bt(i,j,0) = dudx_bt(i,j,0) - dvdy_bt(i,j,0);
    });

    /// Components for the barotropic shearing strain
    ParallelFor(bxQ0_g2ij, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        dvdx_bt(i,j,0) = dy_dxBu(i,j,0)*((vbtav(i+1,j,0)*IdyCv(i+1,j,0))
                                         - (vbtav(i,j,0)*IdyCv(i,j,0)));
        dudy_bt(i,j,0) = dx_dyBu(i,j,0)*((ubtav(i,j+1,0)*IdxCu(i,j+1,0))
                                         - (ubtav(i,j,0)*IdxCu(i,j,0)));
    });

    if (CS.no_slip) {
        ParallelFor(bxH0_m2ij_pij, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            sh_xy_bt(i,j,0) = (2.0-mask2dBu(i,j,0)) * ( dvdx_bt(i,j,0) + dudy_bt(i,j,0) );
        });
    } else {
        ParallelFor(bxH0_m2ij_pij, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            sh_xy_bt(i,j,0) = mask2dBu(i,j,0) * ( dvdx_bt(i,j,0) + dudy_bt(i,j,0) );
        });
    }

    ParallelFor(bxH0_g2ij, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        htot(i,j,0) = 0.0;
        for (int k = 0; k < nz; ++k) {
            htot(i,j,0) = htot(i,j,0) + h(i,j,k);
        }
    });

    /// The inverse of GME tapering scale [H-1 ~> m-1 or m2 kg-1]
    const Real I_GME_h0 = 1.0 / CS.GME_h0;
    /// The nondimensional prefactor multiplying the GME coefficient [nondim]
    const Real GME_efficiency = CS.GME_efficiency;

    ParallelFor(bxQ0_mij_p2ij, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        /// A mask that zeroes out cells with at least one land edge [nondim]
        Real const boundary_mask_h = (mask2dCu(i,j,0) * mask2dCu(i-1,j,0)) * (mask2dCv(i,j,0) * mask2dCv(i,j-1,0));
        Real const dvdx_bt_h = 0.25*((dvdx_bt(i,j,0)+dvdx_bt(i-1,j-1,0)) + (dvdx_bt(i,j-1,0)+dvdx_bt(i-1,j,0)));
        Real const dudy_bt_h = 0.25*((dudy_bt(i,j,0)+dudy_bt(i-1,j-1,0)) + (dudy_bt(i,j-1,0)+dudy_bt(i-1,j,0)));
        /// Magnitude of the barotropic velocity gradient tensor squared at h-points [T-2 ~> s-2]
        Real const grad_vel_mag_bt_h = mask2dT(i,j,0) * boundary_mask_h * (dudx_bt(i,j,0)*dudx_bt(i,j,0) +
              dvdy_bt(i,j,0)*dvdy_bt(i,j,0) + dvdx_bt_h*dvdx_bt_h + dudy_bt_h*dudy_bt_h);
        /// Probably the following test could be simplified to
        /// if (boundary_mask_h * G%mask2dT(I,J) > 0.0) then
        if (grad_vel_mag_bt_h > 0.0) {
            Real const htot_frac = amrex::min(htot(i,j,0) * I_GME_h0, 1.0);
            GME_effic_h(i,j,0) = GME_efficiency * mask2dT(i,j,0) * (htot_frac*htot_frac);
        } else {
            GME_effic_h(i,j,0) = 0.0;
        }
    });

    ParallelFor(bxH0_m2ij_pij, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        /// A mask that zeroes out cells with at least one land edge [nondim]
        Real const boundary_mask_q = (mask2dCv(i,j,0) * mask2dCv(i+1,j,0)) * (mask2dCu(i,j,0) * mask2dCu(i,j+1,0));
        Real const dudx_bt_q = 0.25*((dudx_bt(i,j,0)+dudx_bt(i+1,j+1,0)) + (dudx_bt(i,j+1,0)+dudx_bt(i+1,j,0)));
        Real const dvdy_bt_q = 0.25*((dvdy_bt(i,j,0)+dvdy_bt(i+1,j+1,0)) + (dvdy_bt(i,j+1,0)+dvdy_bt(i+1,j,0)));
        /// Magnitude of the barotropic velocity gradient tensor squared at q-points [T-2 ~> s-2]
        Real const grad_vel_mag_bt_q = mask2dBu(i,j,0) * boundary_mask_q * (dvdx_bt(i,j,0)*dvdx_bt(i,j,0) +
              dudy_bt(i,j,0)*dudy_bt(i,j,0) + dudx_bt_q*dudx_bt_q + dvdy_bt_q*dvdy_bt_q);
        /// Probably the following test could be simplified to
        /// if (boundary_mask_q * G%mask2dBu(I,J) > 0.0) then
        if (grad_vel_mag_bt_q > 0.0) {
            /// The arithmetic mean total thickness at q points [H ~> m or kg m-2]
            Real const h_arith_q = 0.25 * ((htot(i,j,0) + htot(i+1,j+1,0)) + (htot(i+1,j,0) + htot(i,j+1,0)));
            Real const h_arith_frac = amrex::min(h_arith_q * I_GME_h0, 1.0);
            GME_effic_q(i,j,0) = GME_efficiency * mask2dBu(i,j,0) * (h_arith_frac*h_arith_frac);
        } else {
            GME_effic_q(i,j,0) = 0.0;
        }
    });
}

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

void hor_visc_leith_grad_a(const Box& bxH,
                           const Box& bxKh,
                           const Box& bxKh_i,
                           const Box& bxKh_j,
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

    const Box bxKh_gij          = grow_by(bxKh,   1, 1, 1, 1);  ///< bxKh grown to [is_Kh-1:ie_Kh+1, js_Kh-1:je_Kh+1]
    const Box bxKh_i_mi_m2j_pij = grow_by(bxKh_i, 1, 2, 1, 1);  ///< bxKh_i grown to [is_Kh-1:ie_Kh+1, js-2:je+1]
    const Box bxKh_i_mij_pj     = grow_by(bxKh_i, 1, 1, 0, 1);  ///< bxKh_i grown to [is_Kh-1:ie_Kh, js-1:je+1]
    const Box bxKh_j_m2i_mj_pij = grow_by(bxKh_j, 2, 1, 1, 1);  ///< bxKh_j grown to [is-2:ie+1, js_Kh-1:je_Kh+1]
    const Box bxKh_j_mij_pi     = grow_by(bxKh_j, 1, 1, 1, 0);  ///< bxKh_j grown to [is-1:ie+1, js_Kh-1:je_Kh]
    const Box bxKh_mi           = grow_by(bxKh,   1, 0, 0, 0);  ///< bxKh grown to [is_Kh-1:ie_Kh, js_Kh:je_Kh]
    const Box bxKh_mij          = grow_by(bxKh,   1, 1, 0, 0);  ///< bxKh grown to [is_Kh-1:ie_Kh, js_Kh-1:je_Kh]
    const Box bxKh_mj           = grow_by(bxKh,   0, 1, 0, 0);  ///< bxKh grown to [is_Kh:ie_Kh, js_Kh-1:je_Kh]
    const Box bxH_m2i_mj_pij    = grow_by(bxH,    2, 1, 1, 1);  ///< bxH grown to [is-2:ie+1, js-1:je+1]
    const Box bxH_mi_m2j_pij    = grow_by(bxH,    1, 2, 1, 1);  ///< bxH grown to [is-1:ie+1, js-2:je+1]
    const Box bxH_mij           = grow_by(bxH,    1, 1, 0, 0);  ///< bxH grown to [is-1:ie, js-1:je]

    /// Vorticity gradient
    ParallelFor(bxKh_i_mi_m2j_pij, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
    {
        /// Ratio of meridional over zonal grid spacing at vertices [nondim]
        Real const DY_dxBu = dyBu(i,j,0) * IdxBu(i,j,0);
        vort_xy_dx(i,j,kk) = DY_dxBu * ((vort_xy(i,j,kk) * IdyCu(i,j,0)) - (vort_xy(i-1,j,kk) * IdyCu(i-1,j,0)));
    });

    ParallelFor(bxKh_j_m2i_mj_pij, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
    {
        /// Ratio of zonal over meridional grid spacing at vertices [nondim]
        Real const DX_dyBu = dxBu(i,j,0) * IdyBu(i,j,0);
        vort_xy_dy(i,j,kk) = DX_dyBu * ((vort_xy(i,j,kk) * IdxCv(i,j,0)) - (vort_xy(i,j-1,kk) * IdxCv(i,j-1,0)));
    });

    if (CS.use_Leithy) {
        /// Gradient of smoothed vorticity
        ParallelFor(bxKh_mj, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            Real const DY_dxBu = dyBu(i,j,0) * IdxBu(i,j,0);
            vort_xy_dx_smooth(i,j,kk) = DY_dxBu *
                        ((vort_xy_smooth(i,j,kk) * IdyCu(i,j,0)) - (vort_xy_smooth(i-1,j,kk) * IdyCu(i-1,j,0)));
        });

        ParallelFor(bxKh_mi, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            Real const DX_dyBu = dxBu(i,j,0) * IdyBu(i,j,0);
            vort_xy_dy_smooth(i,j,kk) = DX_dyBu *
                        ((vort_xy_smooth(i,j,kk) * IdxCv(i,j,0)) - (vort_xy_smooth(i,j-1,kk) * IdxCv(i,j-1,0)));
        });
    }

    /// Laplacian of vorticity
    /// if (CS%Leith_Ah .or. CS%use_Leithy) then
    ParallelFor(bxKh_mij, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
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
        ParallelFor(bxKh_gij, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            div_xx(i,j,kk) = dudx(i,j,kk) + dvdy(i,j,kk);
        });

        /// Divergence gradient
        ParallelFor(bxKh_i_mij_pj, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            div_xx_dx(i,j,kk) = IdxCu(i,j,0)*(div_xx(i+1,j,kk) - div_xx(i,j,kk));
        });
        ParallelFor(bxKh_j_mij_pi, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            div_xx_dy(i,j,kk) = IdyCv(i,j,0)*(div_xx(i,j+1,kk) - div_xx(i,j,kk));
        });

        /// Magnitude of divergence gradient
        ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            Real const div_xx_dx_h = 0.5*(div_xx_dx(i,j,kk) + div_xx_dx(i-1,j,kk));
            Real const div_xx_dy_h = 0.5*(div_xx_dy(i,j,kk) + div_xx_dy(i,j-1,kk));
            grad_div_mag_h(i,j,kk) = std::sqrt((div_xx_dx_h*div_xx_dx_h) + (div_xx_dy_h*div_xx_dy_h));
        });
        ParallelFor(bxH_mij, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            Real const div_xx_dx_q = 0.5*(div_xx_dx(i,j,kk) + div_xx_dx(i,j+1,kk));
            Real const div_xx_dy_q = 0.5*(div_xx_dy(i,j,kk) + div_xx_dy(i+1,j,kk));
            grad_div_mag_q(i,j,kk) = std::sqrt((div_xx_dx_q*div_xx_dx_q) + (div_xx_dy_q*div_xx_dy_q));
        });

    } else {

        ParallelFor(bxKh_i_mij_pj, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            div_xx_dx(i,j,kk) = 0.0;
        });
        ParallelFor(bxKh_j_mij_pi, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            div_xx_dy(i,j,kk) = 0.0;
        });
        ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            grad_div_mag_h(i,j,kk) = 0.0;
        });
        ParallelFor(bxH_mij, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            grad_div_mag_q(i,j,kk) = 0.0;
        });

    }

    /// Add in beta for the Leith viscosity
    if (CS.use_beta_in_Leith) {
        ParallelFor(bxH_mi_m2j_pij, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            vort_xy_dx(i,j,kk) = vort_xy_dx(i,j,kk) + 0.5 * ( dF_dx(i,j,0) + dF_dx(i,j+1,0));
        });
        ParallelFor(bxH_m2i_mj_pij, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            vort_xy_dy(i,j,kk) = vort_xy_dy(i,j,kk) + 0.5 * ( dF_dy(i,j,0) + dF_dy(i+1,j,0));
        });
    }

    if (CS.use_QG_Leith_visc) {

        ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            Real const vort_xy_dx_h = 0.5*(vort_xy_dx(i,j,kk) + vort_xy_dx(i,j-1,kk));
            Real const vort_xy_dy_h = 0.5*(vort_xy_dy(i,j,kk) + vort_xy_dy(i-1,j,kk));
            grad_vort_mag_h_2d(i,j,kk) = std::sqrt((vort_xy_dx_h*vort_xy_dx_h) + (vort_xy_dy_h*vort_xy_dy_h));
        });
        ParallelFor(bxH_mij, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            Real const vort_xy_dx_q = 0.5*(vort_xy_dx(i,j,kk) + vort_xy_dx(i+1,j,kk));
            Real const vort_xy_dy_q = 0.5*(vort_xy_dy(i,j,kk) + vort_xy_dy(i,j+1,kk));
            grad_vort_mag_q_2d(i,j,kk) = std::sqrt((vort_xy_dx_q*vort_xy_dx_q) + (vort_xy_dy_q*vort_xy_dy_q));
        });
    }
}

void hor_visc_leith_grad_b(const Box& bxH,
                           const Box& bxKh,
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

    const Box bxH_mij = grow_by(bxH, 1, 1, 0, 0);  ///< bxH grown to [is-1:ie, js-1:je]

    ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
    {
        Real const vort_xy_dx_h = 0.5*(vort_xy_dx(i,j,kk) + vort_xy_dx(i,j-1,kk));
        Real const vort_xy_dy_h = 0.5*(vort_xy_dy(i,j,kk) + vort_xy_dy(i-1,j,kk));
        grad_vort_mag_h(i,j,kk) = std::sqrt((vort_xy_dx_h*vort_xy_dx_h) + (vort_xy_dy_h*vort_xy_dy_h));
    });
    ParallelFor(bxH_mij, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
    {
        Real const vort_xy_dx_q = 0.5*(vort_xy_dx(i,j,kk) + vort_xy_dx(i+1,j,kk));
        Real const vort_xy_dy_q = 0.5*(vort_xy_dy(i,j,kk) + vort_xy_dy(i,j+1,kk));
        grad_vort_mag_q(i,j,kk) = std::sqrt((vort_xy_dx_q*vort_xy_dx_q) + (vort_xy_dy_q*vort_xy_dy_q));
    });

    if (CS.use_Leithy) {
        ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            Real const vort_xy_dx_smooth_h = 0.5*(vort_xy_dx_smooth(i,j,kk) + vort_xy_dx_smooth(i,j-1,kk));
            Real const vort_xy_dy_smooth_h = 0.5*(vort_xy_dy_smooth(i,j,kk) + vort_xy_dy_smooth(i-1,j,kk));
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
        ParallelFor(m_leithy_bx, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            m_leithy(i,j,kk) = 0.0;  ///< This is here to initialize domain edge halo values.
        });
    }
    ParallelFor(bxKh, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
    {
        /// Laplacian of vorticity at h-points [L-2 T-1 ~> m-2 s-1]
        Real const Del2vort_h = 0.25 * ((Del2vort_q(i,j,kk) + Del2vort_q(i-1,j-1,kk)) +
                                        (Del2vort_q(i-1,j,kk) + Del2vort_q(i,j-1,kk)));
        /// 2D Leith biharmonic viscosity [L4 T-1 ~> m4 s-1]
        Real const AhLth = Biharm6_const_xx(i,j,0) * inv_PI6 * amrex::Math::abs(Del2vort_h);
        if (AhLth <= Ah_bg_xx(i,j,0)) {
            m_leithy(i,j,kk) = 0.0;
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
        Real const Del2vort_h = 0.25 * ((Del2vort_q(i,j,kk) + Del2vort_q(i-1,j-1,kk)) +
                                        (Del2vort_q(i-1,j,kk) + Del2vort_q(i,j-1,kk)));
        /// 2D Leith+E biharmonic viscosity [L4 T-1 ~> m4 s-1]
        Real const AhLthy = Biharm6_const_xx(i,j,0) * inv_PI6 *
                std::sqrt(amrex::max(0., Del2vort_h*Del2vort_h -
                                         m_leithy(i,j,kk)*(vert_vort_mag_smooth(i,j,kk)*vert_vort_mag_smooth(i,j,kk))));
        Ah(i,j,kk) = amrex::max(Ah_bg_xx(i,j,0), AhLthy);
    });
    if (CS.smooth_Ah) {
        /// Smooth Ah before applying upper bound.  Square Ah, then smooth, then take its square root.
        Box Ah_sq_bx = Box(Ah_sq);
        Ah_sq_bx.setSmall(2, bxKh.smallEnd(2));
        Ah_sq_bx.setBig(2, bxKh.bigEnd(2));
        ParallelFor(Ah_sq_bx, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            Ah_sq(i,j,kk) = 0.0;  ///< This is here to initialize domain edge halo values.
        });
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

} // namespace MOM
