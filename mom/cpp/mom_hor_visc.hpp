// SKILLS: 0.3.3
#pragma once
///
/// @file mom_hor_visc.hpp
/// @brief Box-level AMReX kernel declarations for the MOM6 horizontal viscosity terms.
///
#include <AMReX_Box.H>
#include <AMReX_Array4.H>

#include "mom_hor_visc_kernel.hpp"
#include "turbotmp_bridge_c_types.h"

struct OceanOBC;  // Undefined at the moment

/// @brief Options used by the horizontal viscosity kernels.
/// Field-for-field mirror of the Fortran `bind(C)` type `hv_opts_CS_C` -- order and
/// types must not change.
struct hv_opts_CS_C {
    bool   Laplacian;              ///< Use a Laplacian horizontal viscosity if true.
    bool   biharmonic;             ///< Use a biharmonic horizontal viscosity if true.
    bool   debug;                  ///< If true, write verbose checksums for debugging purposes.
    bool   no_slip;                ///< If true, no slip boundary conditions are used. Otherwise free slip
                                   ///< boundary conditions are assumed. The implementation of the free slip
                                   ///< boundary conditions on a C-grid is much cleaner than the no slip boundary
                                   ///< conditions. The use of free slip b.c.s is strongly encouraged. The no
                                   ///< slip b.c.s are not implemented with the biharmonic viscosity.
    bool   bound_Kh;               ///< If true, the Laplacian coefficient is locally limited to guarantee
                                   ///< stability.
    bool   EY24_EBT_BS;            ///< If true, use an equivalent barotropic backscatter with a stabilizing
                                   ///< kill switch in MEKE, developed by Yankovsky et al. 2024
    bool   bound_Ah;               ///< If true, the biharmonic coefficient is locally limited to guarantee
                                   ///< stability.
    double Re_Ah;                  ///< If nonzero, the biharmonic coefficient is scaled so that the
                                   ///< biharmonic Reynolds number is equal to this [nondim].
    double KS_coef;                ///< A nondimensional coefficient on the biharmonic viscosity that sets the
                                   ///< kill switch for backscatter. Default is 1.0 [nondim].
    bool   backscatter_underbound; ///< If true, the bounds on the biharmonic viscosity are allowed to
                                   ///< increase where the Laplacian viscosity is negative (due to backscatter
                                   ///< parameterizations) beyond the largest timestep-dependent stable values of
                                   ///< biharmonic viscosity when no Laplacian viscosity is applied.  The default
                                   ///< is true for historical reasons, but this option probably should not be
                                   ///< used as it can lead to numerical instabilities.
    bool   Smagorinsky_Kh;         ///< If true, use Smagorinsky nonlinear eddy viscosity. KH is the
                                   ///< background value.
    bool   Smagorinsky_Ah;         ///< If true, use a biharmonic form of Smagorinsky nonlinear eddy
                                   ///< viscosity. AH is the background.
    bool   Leith_Kh;               ///< If true, use 2D Leith nonlinear eddy viscosity. KH is the background
                                   ///< value.
    bool   Modified_Leith;         ///< If true, use extra component of Leith viscosity to damp divergent
                                   ///< flow. To use, still set Leith_Kh=.TRUE.
    bool   use_beta_in_Leith;      ///< If true, includes the beta term in the Leith viscosity
    bool   Leith_Ah;               ///< If true, use a biharmonic form of 2D Leith nonlinear eddy viscosity.
                                   ///< AH is the background.
    bool   use_Leithy;             ///< If true, use a biharmonic form of 2D Leith nonlinear eddy viscosity
                                   ///< with harmonic backscatter. Ah is the background. Leithy = Leith+E
    double c_K;                    ///< Fraction of energy dissipated by the biharmonic term that gets
                                   ///< backscattered in the Leith+E scheme. [nondim]
    bool   smooth_Ah;              ///< If true (default), then Ah and m_leithy are smoothed. This smoothing
                                   ///< requires a lot of blocking communication.
    bool   use_QG_Leith_visc;      ///< If true, use QG Leith nonlinear eddy viscosity. KH is the background
                                   ///< value.
    bool   bound_Coriolis;         ///< If true & SMAGORINSKY_AH is used, the biharmonic viscosity is modified
                                   ///< to include a term that scales quadratically with the velocity shears.
    double Kh_bg_min;              ///< The minimum value allowed for Laplacian horizontal viscosity
                                   ///< [L2 T-1 ~> m2 s-1]. The default is 0.0.
    bool   FrictWork_bug;          ///< If true, retain an answer-changing bug in calculating FrictWork, which
                                   ///< cancels the h in thickness flux and the h at velocity point.
    bool   OBC_strain_bug;         ///< If true, recover a bug that specified shear strain option at open
                                   ///< boundaries cannot be applied.
    bool   use_land_mask;          ///< Use the land mask for the computation of thicknesses at velocity
                                   ///< locations. This eliminates the dependence on arbitrary values over land
                                   ///< or outside of the domain. Default is False to maintain answers with
                                   ///< legacy experiments but should be changed to True for new experiments.
    bool   anisotropic;            ///< If true, allow anisotropic component to the viscosity.
    bool   add_LES_viscosity;      ///< If true, adds the viscosity from Smagorinsky and Leith to the
                                   ///< background viscosity instead of taking the maximum.
    double Kh_aniso;               ///< The anisotropic viscosity [L2 T-1 ~> m2 s-1].
    bool   res_scale_MEKE;         ///< If true, the viscosity contribution from MEKE is scaled by the
                                   ///< resolution function.
    bool   use_GME;                ///< If true, use GME backscatter scheme.
    int    answer_date;            ///< The vintage of the order of arithmetic and expressions in the
                                   ///< horizontal viscosity calculations.  Values below 20190101 recover the
                                   ///< answers from the end of 2018, while higher values use updated and more
                                   ///< robust forms of the same expressions.
    double GME_h0;                 ///< The strength of GME tapers quadratically to zero when the bathymetric
                                   ///< total water column thickness is less than GME_H0 [H ~> m or kg m-2]
    double GME_efficiency;         ///< The nondimensional prefactor multiplying the GME coefficient [nondim]
    double GME_limiter;            ///< The absolute maximum value the GME coefficient is allowed to take
                                   ///< [L2 T-1 ~> m2 s-1].
    double min_grid_Kh;            ///< Minimum horizontal Laplacian viscosity used to limit the grid Reynolds
                                   ///< number [L2 T-1 ~> m2 s-1]
    double min_grid_Ah;            ///< Minimun horizontal biharmonic viscosity used to limit grid Reynolds
                                   ///< number [L4 T-1 ~> m4 s-1]
    bool   use_cont_thick;         ///< If true, thickness at velocity points adopts h[uv] in BT_cont from
                                   ///< continuity solver.
    bool   use_ZB2020;             ///< If true, use Zanna-Bolton 2020 parameterization.
    bool   use_circulation;        ///< If true, use circulation theorem to compute vorticity (for ZB20 or
                                   ///< Leith)
    bool   do_grid_Re_Ah;          ///< If true, the grid_Re_Ah diagnostic is requested.
    bool   do_grid_Re_Kh;          ///< If true, the grid_Re_Kh diagnostic is requested.
    bool   do_Ah_h;                ///< If true, the Ah_h diagnostic is requested.
    bool   do_Ah_q;                ///< If true, the Ah_q diagnostic is requested.
    bool   do_Kh_h;                ///< If true, the Kh_h diagnostic is requested.
    bool   do_Kh_q;                ///< If true, the Kh_q diagnostic is requested.
    bool   do_GME_coeff_h;         ///< If true, the GME_coeff_h diagnostic is requested.
    bool   do_GME_coeff_q;         ///< If true, the GME_coeff_q diagnostic is requested.
    bool   do_vort_xy_q;           ///< If true, the vort_xy_q diagnostic is requested.
    bool   do_div_xx_h;            ///< If true, the div_xx_h diagnostic is requested.
    bool   do_sh_xy_q;             ///< If true, the sh_xy_q diagnostic is requested.
    bool   do_sh_xx_h;             ///< If true, the sh_xx_h diagnostic is requested.
    bool   do_FrictWork;           ///< If true, the FrictWork diagnostic is requested.
    bool   do_FrictWorkIntz;       ///< If true, the FrictWorkIntz diagnostic is requested.
    bool   do_FrictWork_bh;        ///< If true, the FrictWork_bh diagnostic is requested.
    bool   do_FrictWorkIntz_bh;    ///< If true, the FrictWorkIntz_bh diagnostic is requested.
    bool   do_normstress;          ///< If true, the normstress diagnostic is requested.
    bool   do_shearstress;         ///< If true, the shearstress diagnostic is requested.
    bool   do_BS_coeff_h;          ///< If true, the BS_coeff_h diagnostic is requested.
    bool   do_BS_coeff_q;          ///< If true, the BS_coeff_q diagnostic is requested.
};

namespace MOM {
using amrex::Box;
using amrex::Array4;
using amrex::Real;

///
/// @brief Apply a 9-point smoothing filter twice to a field staggered at a thickness point to reduce
///        horizontal two-grid-point noise.
///        Note that this subroutine does not conserve mass, so don't use it in situations where you
///        need conservation.  Also note that it assumes that the input field has valid values in the
///        first two halo points upon entry.
///
/// @param bxH0      The h-point box of the whole column, [isc:iec, jsc:jec, 1:nz]
/// @param mask2dT   0 for land points and 1 for ocean points on the h-grid [nondim]
/// @param field_h   h-point field to be smoothed [arbitrary]
/// @param zero_land If false, return the average of the surrounding ocean points when
///                  smoothing, otherwise use a value of 0 for land points and include them
///                  in the averages.
///
void smooth_x9_h(const Box& bxH0,
                 Array4<const Real> const& mask2dT,
                 Array4<Real>       const& field_h,
                 const bool zero_land);

///
/// @brief Apply a 9-point smoothing filter twice to a pair of velocity components to reduce
///        horizontal two-grid-point noise.
///        Note that this subroutine does not conserve angular momentum, so don't use it
///        in situations where you need conservation.  Also note that it assumes that the
///        input fields have valid values in the first two halo points upon entry.
///
/// @param bxU0      The u-point box of the whole column, [IscB:IecB, jsc:jec, 1:nz]
/// @param bxV0      The v-point box of the whole column, [isc:iec, JscB:JecB, 1:nz]
/// @param mask2dCu  0 for boundary points and 1 for ocean points on the u grid [nondim]
/// @param mask2dCv  0 for boundary points and 1 for ocean points on the v grid [nondim]
/// @param field_u   u-point field to be smoothed [arbitrary]
/// @param field_v   v-point field to be smoothed [arbitrary]
/// @param zero_land If false, return the average of the surrounding ocean points when
///                  smoothing, otherwise use a value of 0 for land points and include them
///                  in the averages.
///
void smooth_x9_uv(const Box& bxU0,
                  const Box& bxV0,
                  Array4<const Real> const& mask2dCu,
                  Array4<const Real> const& mask2dCv,
                  Array4<Real>       const& field_u,
                  Array4<Real>       const& field_v,
                  const bool zero_land);

///
/// @brief Applies one pass of the 1-1-4-1-1 Laplacian filter to a GME diffusive flux at h points,
///        over the computational domain extended by xh points.
///
/// @param bxPass     The h-point range of this smoothing pass, the computational domain extended
///                   by the valid halo points
/// @param mask2dT    0 for land points and 1 for ocean points on the h-grid [nondim]
/// @param GME_flux_h GME diffusive flux at h points [L2 T-2 ~> m2 s-2]
///
void smooth_gme_h_pass(const Box& bxPass,
                       Array4<const Real> const& mask2dT,
                       Array4<Real>       const& GME_flux_h);

///
/// @brief Applies one pass of the 1-1-4-1-1 Laplacian filter to a GME diffusive flux at q points,
///        over the computational domain extended by xq points.
///
/// @param bxPass     The q-point range of this smoothing pass, the computational domain extended
///                   by the valid halo points
/// @param mask2dBu   0 for boundary points and 1 for ocean points on the q grid [nondim]
/// @param GME_flux_q GME diffusive flux at q points [L2 T-2 ~> m2 s-2]
///
void smooth_gme_q_pass(const Box& bxPass,
                       Array4<const Real> const& mask2dBu,
                       Array4<Real>       const& GME_flux_q);

///
/// @brief Adds the MEKE-based (EY24_EBT_BS) backscatter contribution to the diagonal term of
///        the stress tensor at h-points. The caller must only invoke this routine when
///        CS%EY24_EBT_BS is true.
///
/// @param bxQ               The B-grid-index iteration box of this tile, [IscB:IecB, JscB:JecB, 1:kmax]
/// @param CS                Options used by the horizontal viscosity kernels
/// @param MEKE_Ku           The MEKE-derived lateral viscosity coefficient [L2 T-1 ~> m2 s-1];
///                          may be unallocated; pass a default-constructed Array4 when unallocated
/// @param BS_struct         Vertical structure function used in backscatter [nondim];
///                          may be unallocated; pass a default-constructed Array4 when unallocated
/// @param use_kh_struct     If true, shape the backscatter coefficient with VarMix%BS_struct
/// @param kstart            The first absolute k-layer of the current k-block
/// @param visc_limit_h_flag determines whether backscatter is shut off [nondim]
/// @param sh_xx             horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
/// @param str_xx            The diagonal term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2]
/// @param BS_coeff_h        A diagnostic array of the backscatter coefficient [L2 T-1 ~> m2 s-1]
///
void hor_visc_backscatter_h(const Box& bxQ,
                            const hv_opts_CS_C& CS,
                            Array4<const Real> const& MEKE_Ku,
                            Array4<const Real> const& BS_struct,
                            const bool use_kh_struct,
                            const int kstart,
                            Array4<const Real> const& visc_limit_h_flag,
                            Array4<const Real> const& sh_xx,
                            Array4<Real>       const& str_xx,
                            Array4<Real>       const& BS_coeff_h);

///
/// @brief Adds the MEKE-based (EY24_EBT_BS) backscatter contribution to the cross term of
///        the stress tensor at q-points. The caller must only invoke this routine when
///        CS%EY24_EBT_BS is true.
///
/// @param bxH               The h-point iteration box of this tile, [isc:iec, jsc:jec, 1:kmax]
/// @param CS                Options used by the horizontal viscosity kernels
/// @param MEKE_Ku           The MEKE-derived lateral viscosity coefficient [L2 T-1 ~> m2 s-1];
///                          may be unallocated; pass a default-constructed Array4 when unallocated
/// @param BS_struct         Vertical structure function used in backscatter [nondim];
///                          may be unallocated; pass a default-constructed Array4 when unallocated
/// @param use_kh_struct     If true, shape the backscatter coefficient with VarMix%BS_struct
/// @param kstart            The first absolute k-layer of the current k-block
/// @param visc_limit_q_flag determines whether backscatter is shut off [nondim]
/// @param sh_xy             horizontal shearing strain (du/dy + dv/dx) including metric terms [T-1 ~> s-1]
/// @param str_xy            The cross term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2]
/// @param BS_coeff_q        A diagnostic array of the backscatter coefficient [L2 T-1 ~> m2 s-1]
///
void hor_visc_backscatter_q(const Box& bxH,
                            const hv_opts_CS_C& CS,
                            Array4<const Real> const& MEKE_Ku,
                            Array4<const Real> const& BS_struct,
                            const bool use_kh_struct,
                            const int kstart,
                            Array4<const Real> const& visc_limit_q_flag,
                            Array4<const Real> const& sh_xy,
                            Array4<Real>       const& str_xy,
                            Array4<Real>       const& BS_coeff_q);

///
/// @brief Calculates the barotropic tension and shearing strain fields and the GME efficiency
///        fields from the barotropic velocities, for use within horizontal_viscosity when GME
///        (CS%use_GME) is active.
///
/// @param bxH0        The h-point box of the whole column, [isc:iec, jsc:jec, 1:nz]
/// @param bxQ0        The B-grid-index box of the whole column, [IscB:IecB, JscB:JecB, 1:nz]
/// @param CS          Options used by the horizontal viscosity kernels
/// @param IdyCu       1/dyCu [L-1 ~> m-1]
/// @param IdxCv       1/dxCv [L-1 ~> m-1]
/// @param IdyCv       1/dyCv [L-1 ~> m-1]
/// @param IdxCu       1/dxCu [L-1 ~> m-1]
/// @param mask2dBu    0 for boundary points and 1 for ocean points on the q grid [nondim]
/// @param mask2dCu    0 for boundary points and 1 for ocean points on the u grid [nondim]
/// @param mask2dCv    0 for boundary points and 1 for ocean points on the v grid [nondim]
/// @param mask2dT     0 for land points and 1 for ocean points on the h-grid [nondim]
/// @param dx_dyT      Pre-calculated dx/dy at h points [nondim]; may be unallocated; pass a default-constructed Array4
///                    when unallocated
/// @param dy_dxT      Pre-calculated dy/dx at h points [nondim]; may be unallocated; pass a default-constructed Array4
///                    when unallocated
/// @param dx_dyBu     Pre-calculated dx/dy at q points [nondim]; may be unallocated; pass a default-constructed Array4
///                    when unallocated
/// @param dy_dxBu     Pre-calculated dy/dx at q points [nondim]; may be unallocated; pass a default-constructed Array4
///                    when unallocated
/// @param h           Layer thicknesses [H ~> m or kg m-2]
/// @param ubtav       zonal barotropic velocity averaged over a baroclinic time-step [L T-1 ~> m s-1]
/// @param vbtav       meridional barotropic velocity averaged over a baroclinic time-step [L T-1 ~> m s-1]
/// @param dudx_bt     x-component in the barotropic horizontal tension [T-1 ~> s-1]
/// @param dvdy_bt     y-component in the barotropic horizontal tension [T-1 ~> s-1]
/// @param dvdx_bt     x-component in the barotropic shearing strain [T-1 ~> s-1]
/// @param dudy_bt     y-component in the barotropic shearing strain [T-1 ~> s-1]
/// @param sh_xx_bt    Barotropic horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
/// @param sh_xy_bt    Barotropic horizontal shearing strain (du/dy + dv/dx) including metric terms [T-1 ~> s-1]
/// @param GME_effic_h The filtered efficiency of the GME terms at h points [nondim]
/// @param GME_effic_q The filtered efficiency of the GME terms at q points [nondim]
///
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
                     Array4<Real>       const& GME_effic_q);

///
/// @brief Calculates the magnitude of the vorticity and divergence gradients, and the
///        Laplacian of vorticity, that are used within horizontal_viscosity by the Leith,
///        modified Leith, Leith+E, and QG Leith viscosity schemes. The caller must only
///        invoke this routine when CS%Leith_Kh, CS%Leith_Ah, or CS%use_Leithy is true.
///
/// @param bxH                The h-point iteration box of this tile, [isc:iec, jsc:jec, 1:kmax]
/// @param bxKh               The thickness-point viscosity box of this tile, [is_Kh:ie_Kh, js_Kh:je_Kh, 1:kmax]
/// @param bxKh_i             The thickness-point viscosity range in i only, [is_Kh:ie_Kh, jsc:jec, 1:kmax]
/// @param bxKh_j             The thickness-point viscosity range in j only, [isc:iec, js_Kh:je_Kh, 1:kmax]
/// @param CS                 Options used by the horizontal viscosity kernels
/// @param dyBu               dyBu is delta y at q points [L ~> m]
/// @param IdxBu              1/dxBu [L-1 ~> m-1]
/// @param IdyCu              1/dyCu [L-1 ~> m-1]
/// @param dxBu               dxBu is delta x at q points [L ~> m]
/// @param IdyBu              1/dyBu [L-1 ~> m-1]
/// @param IdxCv              1/dxCv [L-1 ~> m-1]
/// @param IdyCv              1/dyCv [L-1 ~> m-1]
/// @param IdxCu              1/dxCu [L-1 ~> m-1]
/// @param dF_dx              Derivative d/dx f (Coriolis parameter) at h-points [T-1 L-1 ~> s-1 m-1]
/// @param dF_dy              Derivative d/dy f (Coriolis parameter) at h-points [T-1 L-1 ~> s-1 m-1]
/// @param dudx               x-derivative of the horizontal tension [T-1 ~> s-1]
/// @param dvdy               y-derivative of the horizontal tension [T-1 ~> s-1]
/// @param vort_xy            Vertical vorticity (dv/dx - du/dy) including metric terms [T-1 ~> s-1]
/// @param vort_xy_smooth     Vertical vorticity including metric terms, smoothed [T-1 ~> s-1]
/// @param grad_vort_mag_h_2d Magnitude of 2d vorticity gradient at h-points [L-1 T-1 ~> m-1 s-1]
/// @param grad_div_mag_h     Magnitude of divergence gradient at h-points [L-1 T-1 ~> m-1 s-1]
/// @param grad_vort_mag_q_2d Magnitude of 2d vorticity gradient at q-points [L-1 T-1 ~> m-1 s-1]
/// @param grad_div_mag_q     Magnitude of divergence gradient at q-points [L-1 T-1 ~> m-1 s-1]
/// @param Del2vort_q         Laplacian of vorticity at q-points [L-2 T-1 ~> m-2 s-1]
/// @param vort_xy_dx         x-derivative of vertical vorticity (d/dx(dv/dx - du/dy)) [L-1 T-1 ~> m-1 s-1]
/// @param vort_xy_dx_smooth  x-derivative of smoothed vertical vorticity [L-1 T-1 ~> m-1 s-1]
/// @param div_xx_dy          y-derivative of horizontal divergence (d/dy(du/dx + dv/dy)) [L-1 T-1 ~> m-1 s-1]
/// @param vort_xy_dy         y-derivative of vertical vorticity (d/dy(dv/dx - du/dy)) [L-1 T-1 ~> m-1 s-1]
/// @param vort_xy_dy_smooth  y-derivative of smoothed vertical vorticity [L-1 T-1 ~> m-1 s-1]
/// @param div_xx_dx          x-derivative of horizontal divergence (d/dx(du/dx + dv/dy)) [L-1 T-1 ~> m-1 s-1]
///
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
                           Array4<Real>       const& div_xx_dx);

///
/// @brief Calculates the magnitudes of the vorticity gradients used within horizontal_viscosity by the
///        Leith viscosity schemes, after any QG Leith contribution has been added.
///
/// @param bxH                  The h-point iteration box of this tile, [isc:iec, jsc:jec, 1:kmax]
/// @param bxKh                 The thickness-point viscosity box of this tile, [is_Kh:ie_Kh, js_Kh:je_Kh, 1:kmax]
/// @param CS                   Options used by the horizontal viscosity kernels
/// @param grad_vort_mag_h      Magnitude of vorticity gradient at h-points [L-1 T-1 ~> m-1 s-1]
/// @param grad_vort_mag_q      Magnitude of vorticity gradient at q-points [L-1 T-1 ~> m-1 s-1]
/// @param vert_vort_mag_smooth Magnitude of gradient of smoothed vertical vorticity (h or q) [L-1 T-1 ~> m-1 s-1]
/// @param vort_xy_dx           x-derivative of vertical vorticity (d/dx(dv/dx - du/dy)) [L-1 T-1 ~> m-1 s-1]
/// @param vort_xy_dx_smooth    x-derivative of smoothed vertical vorticity [L-1 T-1 ~> m-1 s-1]
/// @param vort_xy_dy           y-derivative of vertical vorticity (d/dy(dv/dx - du/dy)) [L-1 T-1 ~> m-1 s-1]
/// @param vort_xy_dy_smooth    y-derivative of smoothed vertical vorticity [L-1 T-1 ~> m-1 s-1]
///
void hor_visc_leith_grad_b(const Box& bxH,
                           const Box& bxKh,
                           const hv_opts_CS_C& CS,
                           Array4<Real>       const& grad_vort_mag_h,
                           Array4<Real>       const& grad_vort_mag_q,
                           Array4<Real>       const& vert_vort_mag_smooth,
                           Array4<const Real> const& vort_xy_dx,
                           Array4<const Real> const& vort_xy_dx_smooth,
                           Array4<const Real> const& vort_xy_dy,
                           Array4<const Real> const& vort_xy_dy_smooth);

///
/// @brief Computes the Leith+E antisymmetric-viscosity ratio m_leithy, before any smoothing. The
///        caller must only invoke this routine when CS%use_Leithy is true.
///
/// @param bxKh             The thickness-point viscosity box of this tile, [is_Kh:ie_Kh, js_Kh:je_Kh, 1:kmax]
/// @param CS               Options used by the horizontal viscosity kernels
/// @param mask2dBu         0 for boundary points and 1 for ocean points on the q grid [nondim]
/// @param Ah_bg_xx         The background biharmonic viscosity at h points [L4 T-1 ~> m4 s-1]; may be unallocated; pass
///                         a default-constructed Array4 when unallocated
/// @param m_const_leithy   Pre-calculated .5*sqrt(c_K)*max{dx,dy} [L ~> m]; may be unallocated; pass a
///                         default-constructed Array4 when unallocated
/// @param m_leithy_max     Pre-calculated 4./max(dx,dy)^2 at h points [L-2 ~> m-2]; may be unallocated; pass a
///                         default-constructed Array4 when unallocated
/// @param Biharm6_const_xx Biharmonic metric-dependent constants [L6 ~> m6]; may be unallocated; pass a
///                         default-constructed Array4 when unallocated
/// @param inv_PI6          The inverse of pi to the sixth power [nondim]
/// @param Del2vort_q       Laplacian of vorticity at q-points [L-2 T-1 ~> m-2 s-1]
/// @param vert_vort_mag    Magnitude of the vertical vorticity gradient (h or q) [L-1 T-1 ~> m-1 s-1]
/// @param vort_xy_smooth   Vertical vorticity including metric terms, smoothed [T-1 ~> s-1]
/// @param m_leithy         Kh=m_leithy*Ah in Leith+E parameterization [L-2 ~> m-2]
///
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
                       Array4<Real>       const& m_leithy);

///
/// @brief Computes the Leith+E antisymmetric-viscosity ratio m_leithy, and updates the
///        biharmonic viscosity Ah (and Ah_h) to include the Leith+E contribution, optionally
///        smoothing m_leithy and Ah in the process. The caller must only invoke this
///        routine when CS%use_Leithy is true.
///
/// @param bxKh                 The thickness-point viscosity box of this tile, [is_Kh:ie_Kh, js_Kh:je_Kh, 1:kmax]
/// @param CS                   Options used by the horizontal viscosity kernels
/// @param Ah_bg_xx             The background biharmonic viscosity at h points [L4 T-1 ~> m4 s-1]; may be unallocated;
///                             pass a default-constructed Array4 when unallocated
/// @param Biharm6_const_xx     Biharmonic metric-dependent constants [L6 ~> m6]; may be unallocated; pass a
///                             default-constructed Array4 when unallocated
/// @param inv_PI6              The inverse of pi to the sixth power [nondim]
/// @param Del2vort_q           Laplacian of vorticity at q-points [L-2 T-1 ~> m-2 s-1]
/// @param vert_vort_mag_smooth Magnitude of gradient of smoothed vertical vorticity (h or q) [L-1 T-1 ~> m-1 s-1]
/// @param Ah                   biharmonic viscosity (h or q) [L4 T-1 ~> m4 s-1]
/// @param m_leithy             Kh=m_leithy*Ah in Leith+E parameterization [L-2 ~> m-2]
/// @param Ah_sq                The square of the biharmonic viscosity [L8 T-2 ~> m8 s-2]
///
void hor_visc_leithy_ah(const Box& bxKh,
                        const hv_opts_CS_C& CS,
                        Array4<const Real> const& Ah_bg_xx,
                        Array4<const Real> const& Biharm6_const_xx,
                        const Real inv_PI6,
                        Array4<const Real> const& Del2vort_q,
                        Array4<const Real> const& vert_vort_mag_smooth,
                        Array4<Real>       const& Ah,
                        Array4<const Real> const& m_leithy,
                        Array4<Real>       const& Ah_sq);

///
/// @brief Sets the Leith+E biharmonic viscosity Ah and Ah_h, from the smoothed square of Ah if
///        CS%smooth_Ah is true. The caller must only invoke this routine when CS%use_Leithy is true.
///
/// @param bxKh      The thickness-point viscosity box of this tile, [is_Kh:ie_Kh, js_Kh:je_Kh, 1:kmax]
/// @param CS        Options used by the horizontal viscosity kernels
/// @param Ah_bg_xx  The background biharmonic viscosity at h points [L4 T-1 ~> m4 s-1]; may be unallocated; pass a
///                  default-constructed Array4 when unallocated
/// @param kstart    The first absolute k-layer of the current k-block
/// @param Ah        biharmonic viscosity (h or q) [L4 T-1 ~> m4 s-1]
/// @param Ah_h      biharmonic viscosity at thickness points [L4 T-1 ~> m4 s-1]
/// @param Ah_sq     The square of the biharmonic viscosity [L8 T-2 ~> m8 s-2]
///
void hor_visc_leithy_ah_h(const Box& bxKh,
                          const hv_opts_CS_C& CS,
                          Array4<const Real> const& Ah_bg_xx,
                          const int kstart,
                          Array4<Real>       const& Ah,
                          Array4<Real>       const& Ah_h,
                          Array4<const Real> const& Ah_sq);

///
/// @brief Computes the horizontal tension, shearing strain and vorticity, their smoothed Leith+E
///        counterparts, the thicknesses at velocity points and, for a biharmonic viscosity, the
///        Laplacian of velocity.
///
/// @param bxH              The h-point iteration box of this tile, [isc:iec, jsc:jec, 1:kmax]
/// @param bxQ              The B-grid-index iteration box of this tile, [IscB:IecB, JscB:JecB, 1:kmax]
/// @param bxU              The u-point iteration box of this tile, [IscB:IecB, jsc:jec, 1:kmax]
/// @param bxV              The v-point iteration box of this tile, [isc:iec, JscB:JecB, 1:kmax]
/// @param bxKh             The thickness-point viscosity box of this tile, [is_Kh:ie_Kh, js_Kh:je_Kh, 1:kmax]
/// @param bxVort           The vorticity box of this tile, [is_vort:ie_vort, js_vort:je_vort, 1:kmax]
/// @param CS               Options used by the horizontal viscosity kernels
/// @param IdyCu            1/dyCu [L-1 ~> m-1]
/// @param IdxCv            1/dxCv [L-1 ~> m-1]
/// @param IdyCv            1/dyCv [L-1 ~> m-1]
/// @param IdxCu            1/dxCu [L-1 ~> m-1]
/// @param mask2dT          0 for land points and 1 for ocean points on the h-grid [nondim]
/// @param dxBu             dxBu is delta x at q points [L ~> m]
/// @param mask2dBu         0 for boundary points and 1 for ocean points on the q grid [nondim]
/// @param IareaBu          IareaBu = 1/areaBu [L-2 ~> m-2]
/// @param dyCv             dyCv is delta y at v points [L ~> m]
/// @param dxCu             dxCu is delta x at u points [L ~> m]
/// @param dx2h             Pre-calculated dx^2 at h points [L2 ~> m2]; may be unallocated; pass a default-constructed
///                         Array4 when unallocated
/// @param dy2h             Pre-calculated dy^2 at h points [L2 ~> m2]; may be unallocated; pass a default-constructed
///                         Array4 when unallocated
/// @param dx_dyT           Pre-calculated dx/dy at h points [nondim]; may be unallocated; pass a default-constructed
///                         Array4 when unallocated
/// @param dy_dxT           Pre-calculated dy/dx at h points [nondim]; may be unallocated; pass a default-constructed
///                         Array4 when unallocated
/// @param dx2q             Pre-calculated dx^2 at q points [L2 ~> m2]; may be unallocated; pass a default-constructed
///                         Array4 when unallocated
/// @param dy2q             Pre-calculated dy^2 at q points [L2 ~> m2]; may be unallocated; pass a default-constructed
///                         Array4 when unallocated
/// @param dx_dyBu          Pre-calculated dx/dy at q points [nondim]; may be unallocated; pass a default-constructed
///                         Array4 when unallocated
/// @param dy_dxBu          Pre-calculated dy/dx at q points [nondim]; may be unallocated; pass a default-constructed
///                         Array4 when unallocated
/// @param Idx2dyCu         1/(dx^2 dy) at u points [L-3 ~> m-3]; may be unallocated; pass a default-constructed Array4
///                         when unallocated
/// @param Idxdy2u          1/(dx dy^2) at u points [L-3 ~> m-3]; may be unallocated; pass a default-constructed Array4
///                         when unallocated
/// @param Idx2dyCv         1/(dx^2 dy) at v points [L-3 ~> m-3]; may be unallocated; pass a default-constructed Array4
///                         when unallocated
/// @param Idxdy2v          1/(dx dy^2) at v points [L-3 ~> m-3]; may be unallocated; pass a default-constructed Array4
///                         when unallocated
/// @param OBC              Open boundary condition structure
/// @param kstart           The first absolute k-layer of the current k-block
/// @param use_vort_xy      True if vort_xy must be computed
/// @param apply_OBC        If true, open boundary conditions are applied on this PE
/// @param apply_OBC_strain If true, the specified open boundary strain configuration is applied
/// @param use_cont_huv     If true, use the continuity solver thicknesses at velocity points
/// @param is_vort          Loop ranges for vorticity terms
/// @param ie_vort          Loop ranges for vorticity terms
/// @param js_vort          Loop ranges for vorticity terms
/// @param je_vort          Loop ranges for vorticity terms
/// @param u                The zonal velocity [L T-1 ~> m s-1]
/// @param v                The meridional velocity [L T-1 ~> m s-1]
/// @param h                Layer thicknesses [H ~> m or kg m-2]
/// @param hu_cont          Layer thickness at u-points [H ~> m or kg m-2]; may be unallocated; pass a
///                         default-constructed Array4 when unallocated
/// @param hv_cont          Layer thickness at v-points [H ~> m or kg m-2]; may be unallocated; pass a
///                         default-constructed Array4 when unallocated
/// @param u_smooth         Zonal velocity, smoothed with a spatial low-pass filter [L T-1 ~> m s-1]
/// @param v_smooth         Meridional velocity, smoothed with a spatial low-pass filter [L T-1 ~> m s-1]
/// @param NoSt             A diagnostic array of normal stress [T-1 ~> s-1]
/// @param ShSt             A diagnostic array of shear stress [T-1 ~> s-1]
/// @param Del2u            The u-component of the Laplacian of velocity [L-1 T-1 ~> m-1 s-1]
/// @param h_u              Thickness interpolated to u points [H ~> m or kg m-2]
/// @param Del2v            The v-component of the Laplacian of velocity [L-1 T-1 ~> m-1 s-1]
/// @param h_v              Thickness interpolated to v points [H ~> m or kg m-2]
/// @param sh_xx            horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
/// @param sh_xx_smooth     horizontal tension from smoothed velocity including metric terms [T-1 ~> s-1]
/// @param dudx             components in the horizontal tension [T-1 ~> s-1]
/// @param dvdy             components in the horizontal tension [T-1 ~> s-1]
/// @param dvdx             components in the shearing strain [T-1 ~> s-1]
/// @param dudy             components in the shearing strain [T-1 ~> s-1]
/// @param dvdx_smooth      components in the shearing strain from smoothed velocity [T-1 ~> s-1]
/// @param dudy_smooth      components in the shearing strain from smoothed velocity [T-1 ~> s-1]
/// @param sh_xy            horizontal shearing strain (du/dy + dv/dx) including metric terms [T-1 ~> s-1]
/// @param sh_xy_smooth     horizontal shearing strain from smoothed velocity including metric terms [T-1 ~> s-1]
/// @param vort_xy          Vertical vorticity (dv/dx - du/dy) including metric terms [T-1 ~> s-1]
/// @param vort_xy_smooth   Vertical vorticity including metric terms, smoothed [T-1 ~> s-1]
///
void hor_visc_strain(const Box& bxH,
                     const Box& bxQ,
                     const Box& bxU,
                     const Box& bxV,
                     const Box& bxKh,
                     const Box& bxVort,
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
                     const int is_vort,
                     const int ie_vort,
                     const int js_vort,
                     const int je_vort,
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
                     Array4<Real>       const& vort_xy_smooth);

///
/// @brief Determines the Laplacian and biharmonic viscosities at h points, the shear magnitude and
///        the related stability bounds, and the Laplacian component of str_xx.
///
/// @param bxH                The h-point iteration box of this tile, [isc:iec, jsc:jec, 1:kmax]
/// @param bxQ                The B-grid-index iteration box of this tile, [IscB:IecB, JscB:JecB, 1:kmax]
/// @param bxKh               The thickness-point viscosity box of this tile, [is_Kh:ie_Kh, js_Kh:je_Kh, 1:kmax]
/// @param CS                 Options used by the horizontal viscosity kernels
/// @param MEKE_Ku            The MEKE-derived lateral viscosity coefficient [L2 T-1 ~> m2 s-1]; may be unallocated;
///                           pass a default-constructed Array4 when unallocated
/// @param Res_fn_h           Non-dimensional function of the ratio the first baroclinic deformation radius to the grid
///                           spacing at h points [nondim]; may be unallocated; pass a default-constructed Array4 when
///                           unallocated
/// @param BS_struct          Vertical structure function used in backscatter [nondim]; may be unallocated; pass a
///                           default-constructed Array4 when unallocated
/// @param Kh_bg_xx           The background Laplacian viscosity at h points [L2 T-1 ~> m2 s-1]; may be unallocated;
///                           pass a default-constructed Array4 when unallocated
/// @param Ah_bg_xx           The background biharmonic viscosity at h points [L4 T-1 ~> m4 s-1]; may be unallocated;
///                           pass a default-constructed Array4 when unallocated
/// @param Kh_Max_xx          The maximum permitted Laplacian viscosity [L2 T-1 ~> m2 s-1]; may be unallocated; pass a
///                           default-constructed Array4 when unallocated
/// @param n1n2_h             Factor n1*n2 in the anisotropic direction tensor at h-points [nondim]; may be unallocated;
///                           pass a default-constructed Array4 when unallocated
/// @param n1n1_m_n2n2_h      Factor n1**2-n2**2 in the anisotropic direction tensor at h-points [nondim]; may be
///                           unallocated; pass a default-constructed Array4 when unallocated
/// @param grid_sp_h2         Harmonic mean of the squares of the grid [L2 ~> m2]; may be unallocated; pass a default-
///                           constructed Array4 when unallocated
/// @param Laplac2_const_xx   Laplacian metric-dependent constants [L2 ~> m2]; may be unallocated; pass a default-
///                           constructed Array4 when unallocated
/// @param Biharm6_const_xx   Biharmonic metric-dependent constants [L6 ~> m6]; may be unallocated; pass a default-
///                           constructed Array4 when unallocated
/// @param Laplac3_const_xx   Laplacian metric-dependent constants [L3 ~> m3]; may be unallocated; pass a default-
///                           constructed Array4 when unallocated
/// @param Biharm_const_xx    Biharmonic metric-dependent constants [L4 ~> m4]; may be unallocated; pass a default-
///                           constructed Array4 when unallocated
/// @param Biharm_const2_xx   Biharmonic metric-dependent constants [T L4 ~> s m4]; may be unallocated; pass a default-
///                           constructed Array4 when unallocated
/// @param kstart             The first absolute k-layer of the current k-block
/// @param h_neglect          thickness so small it can be lost in roundoff and so neglected [H ~> m or kg m-2]
/// @param inv_PI3            Powers of the inverse of pi [nondim]
/// @param inv_PI6            Powers of the inverse of pi [nondim]
/// @param use_kh_struct      If true, the MEKE viscosity contributions are shaped with VarMix%BS_struct
/// @param use_MEKE_Ku        If true, add the Laplacian viscosity derived from MEKE
/// @param rescale_Kh         If true, the Laplacian viscosity is scaled by the resolution function
/// @param u                  The zonal velocity [L T-1 ~> m s-1]
/// @param v                  The meridional velocity [L T-1 ~> m s-1]
/// @param h                  Layer thicknesses [H ~> m or kg m-2]
/// @param Kh_h               Laplacian viscosity at thickness points [L2 T-1 ~> m2 s-1]
/// @param grid_Re_Kh         Grid Reynolds number for Laplacian horizontal viscosity at h points [nondim]
/// @param div_xx_h           horizontal divergence [T-1 ~> s-1]
/// @param sh_xx_h            horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
/// @param h_u                Thickness interpolated to u points [H ~> m or kg m-2]
/// @param h_v                Thickness interpolated to v points [H ~> m or kg m-2]
/// @param sh_xx              horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
/// @param str_xx             str_xx is the diagonal term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2], but at
///                           some points in the code it is not yet layer integrated, so is in [L2 T-2 ~> m2 s-2]
/// @param grad_vort_mag_h    Magnitude of vorticity gradient at h-points [L-1 T-1 ~> m-1 s-1]
/// @param grad_vort_mag_h_2d Magnitude of 2d vorticity gradient at h-points [L-1 T-1 ~> m-1 s-1]
/// @param grad_div_mag_h     Magnitude of divergence gradient at h-points [L-1 T-1 ~> m-1 s-1]
/// @param dudx               components in the horizontal tension [T-1 ~> s-1]
/// @param dvdy               components in the horizontal tension [T-1 ~> s-1]
/// @param sh_xy              horizontal shearing strain (du/dy + dv/dx) including metric terms [T-1 ~> s-1]
/// @param Del2vort_q         Laplacian of vorticity at q-points [L-2 T-1 ~> m-2 s-1]
/// @param Ah                 biharmonic viscosity (h or q) [L4 T-1 ~> m4 s-1]
/// @param Kh                 Laplacian viscosity (h or q) [L2 T-1 ~> m2 s-1]
/// @param Shear_mag          magnitude of the shear (h or q) [T-1 ~> s-1]
/// @param vert_vort_mag      magnitude of the vertical vorticity gradient (h or q) [L-1 T-1 ~> m-1 s-1]
/// @param hrat_min           h_min divided by the thickness at the stress point (h or q) [nondim]
/// @param visc_bound_rem     fraction of overall viscous bounds that remain to be applied (h or q) [nondim]
///
void hor_visc_h_coeffs_a(const Box& bxH,
                         const Box& bxQ,
                         const Box& bxKh,
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
                         Array4<Real>       const& visc_bound_rem);

///
/// @brief Converts the stress tensor components to layer-integrated stresses, takes their divergence to
///        get the horizontal viscous accelerations, and diagnoses the work done by lateral friction and its
///        contributions to the MEKE and SKEBS source terms.
///
/// @param bxH                The h-point iteration box of this tile, [isc:iec, jsc:jec, 1:kmax]
/// @param bxQ                The B-grid-index iteration box of this tile, [IscB:IecB, JscB:JecB, 1:kmax]
/// @param bxU                The u-point iteration box of this tile, [IscB:IecB, jsc:jec, 1:kmax]
/// @param bxV                The v-point iteration box of this tile, [isc:iec, JscB:JecB, 1:kmax]
/// @param CS                 Options used by the horizontal viscosity kernels
/// @param MEKE_mom_src       MEKE source from lateral friction in the momentum equations [R Z L2 T-3 ~> W m-2]; may be
///                           unallocated; pass a default-constructed Array4 when unallocated
/// @param MEKE_mom_src_bh    MEKE source from the biharmonic part of the lateral friction in the momentum equations
///                           [R Z L2 T-3 ~> W m-2]; may be unallocated; pass a default-constructed Array4 when
///                           unallocated
/// @param MEKE_GME_snk       MEKE sink from GME backscatter in the momentum equations [R Z L2 T-3 ~> W m-2]; may be
///                           unallocated; pass a default-constructed Array4 when unallocated
/// @param backscatter_Ro_pow Power in Rossby number function for backscatter [nondim]
/// @param backscatter_Ro_c   Coefficient in Rossby number function for backscatter [nondim]
/// @param skeb_diss          Dissipation rate used to set amplitude of SKEBS [L2 T-3 ~> m2 s-3]; may be unallocated;
///                           pass a default-constructed Array4 when unallocated
/// @param skeb_frict_coef    If skeb_use_frict is true, then skeb_gm_coef * GM_work is added to the dissipation rate
///                           used to set the amplitude of SKEBS [nondim]
/// @param mask2dBu           0 for boundary points and 1 for ocean points on the q grid [nondim]
/// @param IdxCu              1/dxCu [L-1 ~> m-1]
/// @param IdyCu              1/dyCu [L-1 ~> m-1]
/// @param IareaCu            The masked inverse areas of u-grid cells [L-2 ~> m-2]
/// @param IdyCv              1/dyCv [L-1 ~> m-1]
/// @param IdxCv              1/dxCv [L-1 ~> m-1]
/// @param IareaCv            The masked inverse areas of v-grid cells [L-2 ~> m-2]
/// @param IdxT               1/dxT [L-1 ~> m-1]
/// @param IdyT               IdyT is 1/dyT [L-1 ~> m-1]
/// @param IdyBu              1/dyBu [L-1 ~> m-1]
/// @param IdxBu              1/dxBu [L-1 ~> m-1]
/// @param IareaT             1/areaT [L-2 ~> m-2]
/// @param dxCu               dxCu is delta x at u points [L ~> m]
/// @param dyCv               dyCv is delta y at v points [L ~> m]
/// @param CoriolisBu         The Coriolis parameter at corner points [T-1 ~> s-1]
/// @param reduction_xx       The amount by which stresses through h points are reduced due to partial barriers
///                           [nondim]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param reduction_xy       The amount by which stresses through q points are reduced due to partial barriers
///                           [nondim]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param dx2h               Pre-calculated dx^2 at h points [L2 ~> m2]; may be unallocated; pass a default-constructed
///                           Array4 when unallocated
/// @param dy2h               Pre-calculated dy^2 at h points [L2 ~> m2]; may be unallocated; pass a default-constructed
///                           Array4 when unallocated
/// @param dx2q               Pre-calculated dx^2 at q points [L2 ~> m2]; may be unallocated; pass a default-constructed
///                           Array4 when unallocated
/// @param dy2q               Pre-calculated dy^2 at q points [L2 ~> m2]; may be unallocated; pass a default-constructed
///                           Array4 when unallocated
/// @param OBC                Open boundary condition structure
/// @param kstart             The first absolute k-layer of the current k-block
/// @param h_neglect          thickness so small it can be lost in roundoff and so neglected [H ~> m or kg m-2]
/// @param apply_OBC          If true, open boundary conditions are applied on this PE
/// @param find_FrictWork     If true, the work done by lateral friction is calculated
/// @param find_FrictWork_bh  If true, the work done by biharmonic lateral friction is calculated
/// @param backscat_subround  The ratio of f over Shear_mag that is so small that the backscatter calculation gives the
///                           same value as if f were 0 [nondim].
/// @param skeb_use_frict     If true, accumulate the frictional dissipation for stochastic backscatter
/// @param u                  The zonal velocity [L T-1 ~> m s-1]
/// @param v                  The meridional velocity [L T-1 ~> m s-1]
/// @param h                  Layer thicknesses [H ~> m or kg m-2]
/// @param uh                 The zonal volume transport [H L2 T-1 ~> m3 s-1]
/// @param vh                 The meridional volume transport [H L2 T-1 ~> m3 s-1]
/// @param diffu              Zonal acceleration due to convergence of along-coordinate stress tensor [L T-2 ~> m s-2]
/// @param diffv              Meridional acceleration due to convergence of along-coordinate stress tensor
///                           [L T-2 ~> m s-2]
/// @param str_xx_GME         smoothed diagonal term in the stress tensor from GME [L2 T-2 ~> m2 s-2]
/// @param str_xy_GME         smoothed cross term in the stress tensor from GME [L2 T-2 ~> m2 s-2]
/// @param FrictWork          work done by MKE dissipation mechanisms [R Z L2 T-3 ~> W m-2]
/// @param FrictWork_bh       work done by the biharmonic MKE dissipation mechanisms [R Z L2 T-3 ~> W m-2]
/// @param FrictWork_GME      work done by GME [R Z L2 T-3 ~> W m-2]
/// @param visc_limit_h_flag  determines whether backscatter is shut off [nondim]
/// @param h_u                Thickness interpolated to u points [H ~> m or kg m-2]
/// @param h_v                Thickness interpolated to v points [H ~> m or kg m-2]
/// @param sh_xx              horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
/// @param str_xx             str_xx is the diagonal term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2], but at
///                           some points in the code it is not yet layer integrated, so is in [L2 T-2 ~> m2 s-2]
/// @param bhstr_xx           A copy of str_xx that only contains the biharmonic contribution
///                           [H L2 T-2 ~> m3 s-2 or kg s-2]
/// @param sh_xy              horizontal shearing strain (du/dy + dv/dx) including metric terms [T-1 ~> s-1]
/// @param str_xy             str_xy is the cross term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2], but at some
///                           points in the code it is not yet layer integrated, so is in [L2 T-2 ~> m2 s-2]
/// @param bhstr_xy           A copy of str_xy that only contains the biharmonic contribution
///                           [H L2 T-2 ~> m3 s-2 or kg s-2]
/// @param hq                 harmonic mean of the harmonic means of the u- & v point thicknesses [H ~> m or kg m-2]
///                           This form guarantees that hq/hu < 4.
/// @param H_to_RZ            A constant that translates thickness units to the units of mass per unit area
///                           [R Z H-1 ~> kg m-3 or 1]
/// @param s_to_T             A constant that translates time intervals in seconds to the units of time [T s-1 ~> 1]
///
void hor_visc_stress_div(const Box& bxH,
                         const Box& bxQ,
                         const Box& bxU,
                         const Box& bxV,
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
                         const Real s_to_T);

///
/// @brief Completes the biharmonic viscosity at h points and its contribution to str_xx, adds the
///        MEKE backscatter, and computes the q-point shear magnitude, thicknesses and stability bounds.
///
/// @param bxH               The h-point iteration box of this tile, [isc:iec, jsc:jec, 1:kmax]
/// @param bxQ               The B-grid-index iteration box of this tile, [IscB:IecB, JscB:JecB, 1:kmax]
/// @param bxKh              The thickness-point viscosity box of this tile, [is_Kh:ie_Kh, js_Kh:je_Kh, 1:kmax]
/// @param CS                Options used by the horizontal viscosity kernels
/// @param MEKE_Au           The MEKE-derived lateral biharmonic viscosity coefficient [L4 T-1 ~> m4 s-1]; may be
///                          unallocated; pass a default-constructed Array4 when unallocated
/// @param MEKE_Ku           The MEKE-derived lateral viscosity coefficient [L2 T-1 ~> m2 s-1]; may be unallocated; pass
///                          a default-constructed Array4 when unallocated
/// @param BS_struct         Vertical structure function used in backscatter [nondim]; may be unallocated; pass a
///                          default-constructed Array4 when unallocated
/// @param IdyCu             1/dyCu [L-1 ~> m-1]
/// @param IdxCv             1/dxCv [L-1 ~> m-1]
/// @param IdyCv             1/dyCv [L-1 ~> m-1]
/// @param IdxCu             1/dxCu [L-1 ~> m-1]
/// @param mask2dBu          0 for boundary points and 1 for ocean points on the q grid [nondim]
/// @param mask2dCu          0 for boundary points and 1 for ocean points on the u grid [nondim]
/// @param mask2dCv          0 for boundary points and 1 for ocean points on the v grid [nondim]
/// @param reduction_xx      The amount by which stresses through h points are reduced due to partial barriers [nondim];
///                          may be unallocated; pass a default-constructed Array4 when unallocated
/// @param Ah_Max_xx         The maximum permitted biharmonic viscosity [L4 T-1 ~> m4 s-1]; may be unallocated; pass a
///                          default-constructed Array4 when unallocated
/// @param Ah_Max_xx_KS      The maximum permitted biharmonic viscosity for the kill switch [L4 T-1 ~> m4 s-1]; may be
///                          unallocated; pass a default-constructed Array4 when unallocated
/// @param grid_sp_h3        Harmonic mean of the squares of the grid^(3/2) [L3 ~> m3]; may be unallocated; pass a
///                          default-constructed Array4 when unallocated
/// @param dx_dyT            Pre-calculated dx/dy at h points [nondim]; may be unallocated; pass a default-constructed
///                          Array4 when unallocated
/// @param dy_dxT            Pre-calculated dy/dx at h points [nondim]; may be unallocated; pass a default-constructed
///                          Array4 when unallocated
/// @param dx_dyBu           Pre-calculated dx/dy at q points [nondim]; may be unallocated; pass a default-constructed
///                          Array4 when unallocated
/// @param dy_dxBu           Pre-calculated dy/dx at q points [nondim]; may be unallocated; pass a default-constructed
///                          Array4 when unallocated
/// @param Re_Ah_const_xx    Biharmonic metric-dependent constants [L3 ~> m3]; may be unallocated; pass a default-
///                          constructed Array4 when unallocated
/// @param OBC               Open boundary condition structure
/// @param kstart            The first absolute k-layer of the current k-block
/// @param h_neglect         thickness so small it can be lost in roundoff and so neglected [H ~> m or kg m-2]
/// @param h_neglect3        h_neglect^3 [H3 ~> m3 or kg3 m-6]
/// @param apply_OBC         If true, open boundary conditions are applied on this PE
/// @param use_kh_struct     If true, the MEKE viscosity contributions are shaped with VarMix%BS_struct
/// @param use_MEKE_Au       If true, add the biharmonic viscosity derived from MEKE
/// @param u                 The zonal velocity [L T-1 ~> m s-1]
/// @param v                 The meridional velocity [L T-1 ~> m s-1]
/// @param h                 Layer thicknesses [H ~> m or kg m-2]
/// @param Ah_h              biharmonic viscosity at thickness points [L4 T-1 ~> m4 s-1]
/// @param Kh_h              Laplacian viscosity at thickness points [L2 T-1 ~> m2 s-1]
/// @param grid_Re_Ah        Grid Reynolds number for Biharmonic horizontal viscosity at h points [nondim]
/// @param visc_limit_h      Used to stabilize the EY24_EBT_BS backscatter [nondim]
/// @param visc_limit_h_flag determines whether backscatter is shut off [nondim]
/// @param visc_limit_h_frac determines how close backscatter is to shutting off [nondim]
/// @param BS_coeff_h        A diagnostic array of the backscatter coefficient [L2 T-1 ~> m2 s-1]
/// @param Del2u             The u-component of the Laplacian of velocity [L-1 T-1 ~> m-1 s-1]
/// @param h_u               Thickness interpolated to u points [H ~> m or kg m-2]
/// @param Del2v             The v-component of the Laplacian of velocity [L-1 T-1 ~> m-1 s-1]
/// @param h_v               Thickness interpolated to v points [H ~> m or kg m-2]
/// @param sh_xx             horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
/// @param sh_xx_smooth      horizontal tension from smoothed velocity including metric terms [T-1 ~> s-1]
/// @param str_xx            str_xx is the diagonal term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2], but at
///                          some points in the code it is not yet layer integrated, so is in [L2 T-2 ~> m2 s-2]
/// @param bhstr_xx          A copy of str_xx that only contains the biharmonic contribution
///                          [H L2 T-2 ~> m3 s-2 or kg s-2]
/// @param m_leithy          Kh=m_leithy*Ah in Leith+E parameterization [L-2 ~> m-2]
/// @param dDel2vdx          Components in the biharmonic equivalent of the shearing strain [L-2 T-1 ~> m-2 s-1]
/// @param dDel2udy          Components in the biharmonic equivalent of the shearing strain [L-2 T-1 ~> m-2 s-1]
/// @param sh_xy             horizontal shearing strain (du/dy + dv/dx) including metric terms [T-1 ~> s-1]
/// @param hq                harmonic mean of the harmonic means of the u- & v point thicknesses [H ~> m or kg m-2] This
///                          form guarantees that hq/hu < 4.
/// @param Ah                biharmonic viscosity (h or q) [L4 T-1 ~> m4 s-1]
/// @param Kh                Laplacian viscosity (h or q) [L2 T-1 ~> m2 s-1]
/// @param Shear_mag         magnitude of the shear (h or q) [T-1 ~> s-1]
/// @param hrat_min          h_min divided by the thickness at the stress point (h or q) [nondim]
/// @param visc_bound_rem    fraction of overall viscous bounds that remain to be applied (h or q) [nondim]
///
void hor_visc_h_coeffs_b(const Box& bxH,
                         const Box& bxQ,
                         const Box& bxKh,
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
                         Array4<const Real> const& visc_bound_rem);

///
/// @brief Determines the Laplacian and biharmonic viscosities at q points and the harmonic and
///        biharmonic parts of str_xy, adds the MEKE backscatter, and computes the GME stresses.
///
/// @param bxH                The h-point iteration box of this tile, [isc:iec, jsc:jec, 1:kmax]
/// @param bxQ                The B-grid-index iteration box of this tile, [IscB:IecB, JscB:JecB, 1:kmax]
/// @param CS                 Options used by the horizontal viscosity kernels
/// @param MEKE_Ku            The MEKE-derived lateral viscosity coefficient [L2 T-1 ~> m2 s-1]; may be unallocated;
///                           pass a default-constructed Array4 when unallocated
/// @param MEKE_Au            The MEKE-derived lateral biharmonic viscosity coefficient [L4 T-1 ~> m4 s-1]; may be
///                           unallocated; pass a default-constructed Array4 when unallocated
/// @param Res_fn_q           Non-dimensional function of the ratio the first baroclinic deformation radius to the grid
///                           spacing at q points [nondim]; may be unallocated; pass a default-constructed Array4 when
///                           unallocated
/// @param BS_struct          Vertical structure function used in backscatter [nondim]; may be unallocated; pass a
///                           default-constructed Array4 when unallocated
/// @param mask2dBu           0 for boundary points and 1 for ocean points on the q grid [nondim]
/// @param Kh_bg_xy           The background Laplacian viscosity at q points [L2 T-1 ~> m2 s-1]; may be unallocated;
///                           pass a default-constructed Array4 when unallocated
/// @param Ah_bg_xy           The background biharmonic viscosity at q points [L4 T-1 ~> m4 s-1]; may be unallocated;
///                           pass a default-constructed Array4 when unallocated
/// @param reduction_xy       The amount by which stresses through q points are reduced due to partial barriers
///                           [nondim]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param Kh_Max_xy          The maximum permitted Laplacian viscosity [L2 T-1 ~> m2 s-1]; may be unallocated; pass a
///                           default-constructed Array4 when unallocated
/// @param Ah_Max_xy          The maximum permitted biharmonic viscosity [L4 T-1 ~> m4 s-1]; may be unallocated; pass a
///                           default-constructed Array4 when unallocated
/// @param Ah_Max_xy_KS       The maximum permitted biharmonic viscosity for the kill switch [L4 T-1 ~> m4 s-1]; may be
///                           unallocated; pass a default-constructed Array4 when unallocated
/// @param n1n2_q             Factor n1*n2 in the anisotropic direction tensor at q-points [nondim]; may be unallocated;
///                           pass a default-constructed Array4 when unallocated
/// @param n1n1_m_n2n2_q      Factor n1**2-n2**2 in the anisotropic direction tensor at q-points [nondim]; may be
///                           unallocated; pass a default-constructed Array4 when unallocated
/// @param Laplac2_const_xy   Laplacian metric-dependent constants [L2 ~> m2]; may be unallocated; pass a default-
///                           constructed Array4 when unallocated
/// @param Biharm6_const_xy   Biharmonic metric-dependent constants [L6 ~> m6]; may be unallocated; pass a default-
///                           constructed Array4 when unallocated
/// @param Laplac3_const_xy   Laplacian metric-dependent constants [L3 ~> m3]; may be unallocated; pass a default-
///                           constructed Array4 when unallocated
/// @param Biharm_const_xy    Biharmonic metric-dependent constants [L4 ~> m4]; may be unallocated; pass a default-
///                           constructed Array4 when unallocated
/// @param Biharm_const2_xy   Biharmonic metric-dependent constants [T L4 ~> s m4]; may be unallocated; pass a default-
///                           constructed Array4 when unallocated
/// @param Re_Ah_const_xy     Biharmonic metric-dependent constants [L3 ~> m3]; may be unallocated; pass a default-
///                           constructed Array4 when unallocated
/// @param kstart             The first absolute k-layer of the current k-block
/// @param inv_PI3            Powers of the inverse of pi [nondim]
/// @param inv_PI6            Powers of the inverse of pi [nondim]
/// @param find_FrictWork     If true, the work done by lateral friction is calculated
/// @param use_kh_struct      If true, the MEKE viscosity contributions are shaped with VarMix%BS_struct
/// @param use_MEKE_Ku        If true, add the Laplacian viscosity derived from MEKE
/// @param use_MEKE_Au        If true, add the biharmonic viscosity derived from MEKE
/// @param rescale_Kh         If true, the Laplacian viscosity is scaled by the resolution function
/// @param u                  The zonal velocity [L T-1 ~> m s-1]
/// @param v                  The meridional velocity [L T-1 ~> m s-1]
/// @param sh_xx_bt           barotropic horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
/// @param sh_xy_bt           barotropic horizontal shearing strain (du/dy + dv/dx) inc. metric terms [T-1 ~> s-1]
/// @param GME_effic_h        The filtered efficiency of the GME terms at h points [nondim]
/// @param GME_effic_q        The filtered efficiency of the GME terms at q points [nondim]
/// @param KH_u_GME           Isopycnal height diffusivities in u-columns [L2 T-1 ~> m2 s-1]
/// @param KH_v_GME           Isopycnal height diffusivities in v-columns [L2 T-1 ~> m2 s-1]
/// @param GME_coeff_h        GME coefficient at h-points [L2 T-1 ~> m2 s-1]
/// @param GME_coeff_q        GME coeff. at q-points [L2 T-1 ~> m2 s-1]
/// @param str_xx_GME         smoothed diagonal term in the stress tensor from GME [L2 T-2 ~> m2 s-2]
/// @param str_xy_GME         smoothed cross term in the stress tensor from GME [L2 T-2 ~> m2 s-2]
/// @param Ah_h               biharmonic viscosity at thickness points [L4 T-1 ~> m4 s-1]
/// @param Ah_q               biharmonic viscosity at corner points [L4 T-1 ~> m4 s-1]
/// @param Kh_h               Laplacian viscosity at thickness points [L2 T-1 ~> m2 s-1]
/// @param Kh_q               Laplacian viscosity at corner points [L2 T-1 ~> m2 s-1]
/// @param vort_xy_q          vertical vorticity at corner points [T-1 ~> s-1]
/// @param sh_xy_q            horizontal shearing strain at corner points [T-1 ~> s-1]
/// @param visc_limit_q       used to stabilize the EY24_EBT_BS backscatter [nondim]
/// @param visc_limit_q_flag  determines whether backscatter is shut off [nondim]
/// @param visc_limit_q_frac  determines how close backscatter is to shutting off [nondim]
/// @param BS_coeff_q         A diagnostic array of the backscatter coefficient [L2 T-1 ~> m2 s-1]
/// @param sh_xx              horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
/// @param dDel2vdx           Components in the biharmonic equivalent of the shearing strain [L-2 T-1 ~> m-2 s-1]
/// @param dDel2udy           Components in the biharmonic equivalent of the shearing strain [L-2 T-1 ~> m-2 s-1]
/// @param sh_xy              horizontal shearing strain (du/dy + dv/dx) including metric terms [T-1 ~> s-1]
/// @param sh_xy_smooth       horizontal shearing strain from smoothed velocity including metric terms [T-1 ~> s-1]
/// @param str_xy             str_xy is the cross term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2], but at some
///                           points in the code it is not yet layer integrated, so is in [L2 T-2 ~> m2 s-2]
/// @param bhstr_xy           A copy of str_xy that only contains the biharmonic contribution
///                           [H L2 T-2 ~> m3 s-2 or kg s-2]
/// @param vort_xy            Vertical vorticity (dv/dx - du/dy) including metric terms [T-1 ~> s-1]
/// @param grad_vort_mag_q    Magnitude of vorticity gradient at q-points [L-1 T-1 ~> m-1 s-1]
/// @param grad_vort_mag_q_2d Magnitude of 2d vorticity gradient at q-points [L-1 T-1 ~> m-1 s-1]
/// @param Del2vort_q         Laplacian of vorticity at q-points [L-2 T-1 ~> m-2 s-1]
/// @param grad_div_mag_q     Magnitude of divergence gradient at q-points [L-1 T-1 ~> m-1 s-1]
/// @param hq                 harmonic mean of the harmonic means of the u- & v point thicknesses [H ~> m or kg m-2]
///                           This form guarantees that hq/hu < 4.
/// @param Ah                 biharmonic viscosity (h or q) [L4 T-1 ~> m4 s-1]
/// @param Kh                 Laplacian viscosity (h or q) [L2 T-1 ~> m2 s-1]
/// @param Shear_mag          magnitude of the shear (h or q) [T-1 ~> s-1]
/// @param vert_vort_mag      magnitude of the vertical vorticity gradient (h or q) [L-1 T-1 ~> m-1 s-1]
/// @param hrat_min           h_min divided by the thickness at the stress point (h or q) [nondim]
/// @param visc_bound_rem     fraction of overall viscous bounds that remain to be applied (h or q) [nondim]
///
void hor_visc_q_coeffs(const Box& bxH,
                       const Box& bxQ,
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
                       Array4<Real>       const& visc_bound_rem);

/// @brief Device views of the core ocean grid fields used by the horizontal viscosity kernels.
struct hv_grid_dev {
    Array4<const Real> mask2dT;    ///< 0 for land points and 1 for ocean points on the h-grid [nondim].
    Array4<const Real> IareaT;     ///< 1/areaT [L-2 ~> m-2].
    Array4<const Real> IdxT;       ///< 1/dxT [L-1 ~> m-1].
    Array4<const Real> IdyT;       ///< IdyT is 1/dyT [L-1 ~> m-1].
    Array4<const Real> mask2dCu;   ///< 0 for boundary points and 1 for ocean points on the u grid [nondim].
    Array4<const Real> IdxCu;      ///< 1/dxCu [L-1 ~> m-1].
    Array4<const Real> dxCu;       ///< dxCu is delta x at u points [L ~> m].
    Array4<const Real> mask2dCv;   ///< 0 for boundary points and 1 for ocean points on the v grid [nondim].
    Array4<const Real> IdyCv;      ///< 1/dyCv [L-1 ~> m-1].
    Array4<const Real> dyCv;       ///< dyCv is delta y at v points [L ~> m].
    Array4<const Real> CoriolisBu; ///< The Coriolis parameter at corner points [T-1 ~> s-1].
    Array4<const Real> IdyCu;      ///< 1/dyCu [L-1 ~> m-1].
    Array4<const Real> IdxCv;      ///< 1/dxCv [L-1 ~> m-1].
    Array4<const Real> mask2dBu;   ///< 0 for boundary points and 1 for ocean points on the q grid [nondim].
    Array4<const Real> dxBu;       ///< dxBu is delta x at q points [L ~> m].
    Array4<const Real> dyBu;       ///< dyBu is delta y at q points [L ~> m].
    Array4<const Real> IareaBu;    ///< IareaBu = 1/areaBu [L-2 ~> m-2].
    Array4<const Real> df_dx;      ///< Derivative d/dx f (Coriolis parameter) at h-points [T-1 L-1 ~> s-1 m-1].
    Array4<const Real> df_dy;      ///< Derivative d/dy f (Coriolis parameter) at h-points [T-1 L-1 ~> s-1 m-1].
    Array4<const Real> IareaCu;    ///< The masked inverse areas of u-grid cells [L-2 ~> m-2].
    Array4<const Real> IareaCv;    ///< The masked inverse areas of v-grid cells [L-2 ~> m-2].
    Array4<const Real> IdxBu;      ///< 1/dxBu [L-1 ~> m-1].
    Array4<const Real> IdyBu;      ///< 1/dyBu [L-1 ~> m-1].
};

/// @brief Device views of the horizontal viscosity metric and coefficient arrays; each may be unallocated
///        (a default-constructed Array4).
struct hv_metrics_dev {
    Array4<const Real> Kh_bg_xx;         ///< The background Laplacian viscosity at h points [L2 T-1 ~> m2 s-1]. The
                                         ///< actual viscosity may be the larger of this viscosity and the Smagorinsky
                                         ///< and Leith viscosities.
    Array4<const Real> Ah_bg_xx;         ///< The background biharmonic viscosity at h points [L4 T-1 ~> m4 s-1]. The
                                         ///< actual viscosity may be the larger of this viscosity and the Smagorinsky
                                         ///< and Leith viscosities.
    Array4<const Real> reduction_xx;     ///< The amount by which stresses through h points are reduced due to partial
                                         ///< barriers [nondim].
    Array4<const Real> Kh_Max_xx;        ///< The maximum permitted Laplacian viscosity [L2 T-1 ~> m2 s-1].
    Array4<const Real> Ah_Max_xx;        ///< The maximum permitted biharmonic viscosity [L4 T-1 ~> m4 s-1].
    Array4<const Real> Ah_Max_xx_KS;     ///< The maximum permitted biharmonic viscosity for the kill switch
                                         ///< [L4 T-1 ~> m4 s-1].
    Array4<const Real> n1n2_h;           ///< Factor n1*n2 in the anisotropic direction tensor at h-points [nondim]
    Array4<const Real> n1n1_m_n2n2_h;    ///< Factor n1**2-n2**2 in the anisotropic direction tensor at h-points
                                         ///< [nondim]
    Array4<const Real> grid_sp_h2;       ///< Harmonic mean of the squares of the grid [L2 ~> m2]
    Array4<const Real> grid_sp_h3;       ///< Harmonic mean of the squares of the grid^(3/2) [L3 ~> m3]
    Array4<const Real> Kh_bg_xy;         ///< The background Laplacian viscosity at q points [L2 T-1 ~> m2 s-1]. The
                                         ///< actual viscosity may be the larger of this viscosity and the Smagorinsky
                                         ///< and Leith viscosities.
    Array4<const Real> Ah_bg_xy;         ///< The background biharmonic viscosity at q points [L4 T-1 ~> m4 s-1]. The
                                         ///< actual viscosity may be the larger of this viscosity and the Smagorinsky
                                         ///< and Leith viscosities.
    Array4<const Real> reduction_xy;     ///< The amount by which stresses through q points are reduced due to partial
                                         ///< barriers [nondim].
    Array4<const Real> Kh_Max_xy;        ///< The maximum permitted Laplacian viscosity [L2 T-1 ~> m2 s-1].
    Array4<const Real> Ah_Max_xy;        ///< The maximum permitted biharmonic viscosity [L4 T-1 ~> m4 s-1].
    Array4<const Real> Ah_Max_xy_KS;     ///< The maximum permitted biharmonic viscosity for the  kill switch
                                         ///< [L4 T-1 ~> m4 s-1].
    Array4<const Real> n1n2_q;           ///< Factor n1*n2 in the anisotropic direction tensor at q-points [nondim]
    Array4<const Real> n1n1_m_n2n2_q;    ///< Factor n1**2-n2**2 in the anisotropic direction tensor at q-points
                                         ///< [nondim]
    Array4<const Real> dx2h;             ///< Pre-calculated dx^2 at h points [L2 ~> m2]
    Array4<const Real> dy2h;             ///< Pre-calculated dy^2 at h points [L2 ~> m2]
    Array4<const Real> dx_dyT;           ///< Pre-calculated dx/dy at h points [nondim]
    Array4<const Real> dy_dxT;           ///< Pre-calculated dy/dx at h points [nondim]
    Array4<const Real> m_const_leithy;   ///< Pre-calculated .5*sqrt(c_K)*max{dx,dy} [L ~> m]
    Array4<const Real> m_leithy_max;     ///< Pre-calculated 4./max(dx,dy)^2 at h points [L-2 ~> m-2]
    Array4<const Real> dx2q;             ///< Pre-calculated dx^2 at q points [L2 ~> m2]
    Array4<const Real> dy2q;             ///< Pre-calculated dy^2 at q points [L2 ~> m2]
    Array4<const Real> dx_dyBu;          ///< Pre-calculated dx/dy at q points [nondim]
    Array4<const Real> dy_dxBu;          ///< Pre-calculated dy/dx at q points [nondim]
    Array4<const Real> Idx2dyCu;         ///< 1/(dx^2 dy) at u points [L-3 ~> m-3]
    Array4<const Real> Idxdy2u;          ///< 1/(dx dy^2) at u points [L-3 ~> m-3]
    Array4<const Real> Idx2dyCv;         ///< 1/(dx^2 dy) at v points [L-3 ~> m-3]
    Array4<const Real> Idxdy2v;          ///< 1/(dx dy^2) at v points [L-3 ~> m-3]
    Array4<const Real> Laplac2_const_xx; ///< Laplacian metric-dependent constants [L2 ~> m2]
    Array4<const Real> Biharm6_const_xx; ///< Biharmonic metric-dependent constants [L6 ~> m6]
    Array4<const Real> Laplac3_const_xx; ///< Laplacian metric-dependent constants [L3 ~> m3]
    Array4<const Real> Biharm_const_xx;  ///< Biharmonic metric-dependent constants [L4 ~> m4]
    Array4<const Real> Biharm_const2_xx; ///< Biharmonic metric-dependent constants [T L4 ~> s m4]
    Array4<const Real> Re_Ah_const_xx;   ///< Biharmonic metric-dependent constants [L3 ~> m3]
    Array4<const Real> Laplac2_const_xy; ///< Laplacian metric-dependent constants [L2 ~> m2]
    Array4<const Real> Biharm6_const_xy; ///< Biharmonic metric-dependent constants [L6 ~> m6]
    Array4<const Real> Laplac3_const_xy; ///< Laplacian metric-dependent constants [L3 ~> m3]
    Array4<const Real> Biharm_const_xy;  ///< Biharmonic metric-dependent constants [L4 ~> m4]
    Array4<const Real> Biharm_const2_xy; ///< Biharmonic metric-dependent constants [T L4 ~> s m4]
    Array4<const Real> Re_Ah_const_xy;   ///< Biharmonic metric-dependent constants [L3 ~> m3]
};

/// @brief Device views of the MEKE, VarMix and STOCH fields used by the horizontal viscosity kernels.
struct hv_mixing_dev {
    Array4<const Real> MEKE_Ku;   ///< The MEKE-derived lateral viscosity coefficient [L2 T-1 ~> m2 s-1]; may be
                                  ///< unallocated
    Array4<const Real> MEKE_Au;   ///< The MEKE-derived lateral biharmonic viscosity coefficient [L4 T-1 ~> m4 s-1]; may
                                  ///< be unallocated
    Array4<Real> MEKE_mom_src;    ///< MEKE source from lateral friction in the momentum equations
                                  ///< [R Z L2 T-3 ~> W m-2]; may be unallocated
    Array4<Real> MEKE_mom_src_bh; ///< MEKE source from the biharmonic part of the lateral friction in the momentum
                                  ///< equations [R Z L2 T-3 ~> W m-2]; may be unallocated
    Array4<Real> MEKE_GME_snk;    ///< MEKE sink from GME backscatter in the momentum equations [R Z L2 T-3 ~> W m-2];
                                  ///< may be unallocated
    Real backscatter_Ro_pow;      ///< Power in Rossby number function for backscatter [nondim].
    Real backscatter_Ro_c;        ///< Coefficient in Rossby number function for backscatter [nondim].
    Array4<const Real> Res_fn_h;  ///< Non-dimensional function of the ratio the first baroclinic deformation radius to
                                  ///< the grid spacing at h points [nondim]; may be unallocated
    Array4<const Real> Res_fn_q;  ///< Non-dimensional function of the ratio the first baroclinic deformation radius to
                                  ///< the grid spacing at q points [nondim]; may be unallocated
    Array4<const Real> BS_struct; ///< Vertical structure function used in backscatter [nondim]; may be unallocated
    bool use_variable_mixing;     ///< If true, use the variable mixing.
    bool Resoln_scaled_Kh;        ///< If true, scale away the Laplacian viscosity when the deformation radius is well
                                  ///< resolved.
    Array4<Real> skeb_diss;       ///< Dissipation rate used to set amplitude of SKEBS [L2 T-3 ~> m2 s-3]; may be
                                  ///< unallocated
    bool skeb_use_frict;          ///< If true, adds viscous dissipation rate to the amplitude of SKEBS
    Real skeb_frict_coef;         ///< If skeb_use_frict is true, then skeb_gm_coef * GM_work is added to the
                                  ///< dissipation rate used to set the amplitude of SKEBS [nondim]
};

///
/// @brief Computes the convergence of the horizontal viscous stresses (diffu, diffv) and the
///        associated diagnostics, block by block in the vertical.
///
/// @param bxH0              The h-point box of the whole column, [isc:iec, jsc:jec, 1:nz]
/// @param bxQ0              The B-grid-index box of the whole column, [IscB:IecB, JscB:JecB, 1:nz]
/// @param CS                Options used by the horizontal viscosity kernels
/// @param G                 Device views of the core ocean grid fields
/// @param M                 Device views of the horizontal viscosity metric and coefficient arrays
/// @param X                 Device views of the MEKE, VarMix and STOCH fields
/// @param OBC               Open boundary condition structure
/// @param nkblock           The effective k-block size [nondim]
/// @param H_subroundoff     A thickness that is so small that it can be added to a thickness of Angstrom or larger
///                          without changing it at the bit level [H ~> m or kg m-2]
/// @param H_to_RZ           A constant that translates thickness units to the units of mass per unit area
///                          [R Z H-1 ~> kg m-3 or 1]
/// @param s_to_T            A constant that translates time intervals in seconds to the units of time [T s-1 ~> 1]
/// @param u                 The zonal velocity [L T-1 ~> m s-1]
/// @param v                 The meridional velocity [L T-1 ~> m s-1]
/// @param h                 Layer thicknesses [H ~> m or kg m-2]
/// @param uh                The zonal volume transport [H L2 T-1 ~> m3 s-1]
/// @param vh                The meridional volume transport [H L2 T-1 ~> m3 s-1]
/// @param diffu             Zonal acceleration due to convergence of along-coordinate stress tensor [L T-2 ~> m s-2]
/// @param diffv             Meridional acceleration due to convergence of along-coordinate stress tensor
///                          [L T-2 ~> m s-2]
/// @param hu_cont           Layer thickness at u-points [H ~> m or kg m-2]; may be unallocated; pass a default-
///                          constructed Array4 when unallocated
/// @param hv_cont           Layer thickness at v-points [H ~> m or kg m-2]; may be unallocated; pass a default-
///                          constructed Array4 when unallocated
/// @param u_smooth          Zonal velocity, smoothed with a spatial low-pass filter [L T-1 ~> m s-1]
/// @param v_smooth          Meridional velocity, smoothed with a spatial low-pass filter [L T-1 ~> m s-1]
/// @param sh_xx_bt          barotropic horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
/// @param sh_xy_bt          barotropic horizontal shearing strain (du/dy + dv/dx) inc. metric terms [T-1 ~> s-1]
/// @param GME_effic_h       The filtered efficiency of the GME terms at h points [nondim]
/// @param GME_effic_q       The filtered efficiency of the GME terms at q points [nondim]
/// @param KH_u_GME          Isopycnal height diffusivities in u-columns [L2 T-1 ~> m2 s-1]
/// @param KH_v_GME          Isopycnal height diffusivities in v-columns [L2 T-1 ~> m2 s-1]
/// @param GME_coeff_h       GME coefficient at h-points [L2 T-1 ~> m2 s-1]
/// @param GME_coeff_q       GME coeff. at q-points [L2 T-1 ~> m2 s-1]
/// @param str_xx_GME        smoothed diagonal term in the stress tensor from GME [L2 T-2 ~> m2 s-2]
/// @param str_xy_GME        smoothed cross term in the stress tensor from GME [L2 T-2 ~> m2 s-2]
/// @param NoSt              A diagnostic array of normal stress [T-1 ~> s-1]
/// @param ShSt              A diagnostic array of shear stress [T-1 ~> s-1]
/// @param FrictWork         work done by MKE dissipation mechanisms [R Z L2 T-3 ~> W m-2]
/// @param FrictWork_bh      work done by the biharmonic MKE dissipation mechanisms [R Z L2 T-3 ~> W m-2]
/// @param FrictWork_GME     work done by GME [R Z L2 T-3 ~> W m-2]
/// @param Ah_h              biharmonic viscosity at thickness points [L4 T-1 ~> m4 s-1]
/// @param Ah_q              biharmonic viscosity at corner points [L4 T-1 ~> m4 s-1]
/// @param Kh_h              Laplacian viscosity at thickness points [L2 T-1 ~> m2 s-1]
/// @param Kh_q              Laplacian viscosity at corner points [L2 T-1 ~> m2 s-1]
/// @param grid_Re_Kh        Grid Reynolds number for Laplacian horizontal viscosity at h points [nondim]
/// @param grid_Re_Ah        Grid Reynolds number for Biharmonic horizontal viscosity at h points [nondim]
/// @param div_xx_h          horizontal divergence [T-1 ~> s-1]
/// @param sh_xx_h           horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
/// @param vort_xy_q         vertical vorticity at corner points [T-1 ~> s-1]
/// @param sh_xy_q           horizontal shearing strain at corner points [T-1 ~> s-1]
/// @param visc_limit_h      Used to stabilize the EY24_EBT_BS backscatter [nondim]
/// @param visc_limit_q      used to stabilize the EY24_EBT_BS backscatter [nondim]
/// @param visc_limit_h_flag determines whether backscatter is shut off [nondim]
/// @param visc_limit_q_flag determines whether backscatter is shut off [nondim]
/// @param visc_limit_h_frac determines how close backscatter is to shutting off [nondim]
/// @param visc_limit_q_frac determines how close backscatter is to shutting off [nondim]
/// @param BS_coeff_h        A diagnostic array of the backscatter coefficient [L2 T-1 ~> m2 s-1]
/// @param BS_coeff_q        A diagnostic array of the backscatter coefficient [L2 T-1 ~> m2 s-1]
///
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
                             Array4<Real>       const& BS_coeff_q);

} // namespace MOM
