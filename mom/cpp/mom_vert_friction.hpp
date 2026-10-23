// SKILLS: 0.3.3
#pragma once
///
/// @file mom_vert_friction.hpp
/// @brief Box-level AMReX kernel declarations for the MOM6 vertical friction (viscosity) terms.
///
#include <AMReX_Box.H>
#include <AMReX_Array4.H>

#include "mom_vert_friction_kernel.hpp"
#include "turbotmp_bridge_c_types.h"

/// @brief Options used in the vertical viscosity solver and the velocity truncations.
/// Field-for-field mirror of the Fortran `bind(C)` type `vertvisc_opts_CS_C` -- order and
/// types must not change.
struct vertvisc_opts_CS_C {
    double Hmix_stress;        ///< The mixed layer thickness over which the wind
                               ///< stress is applied with direct_stress [H ~> m or kg m-2].
    bool   direct_stress;      ///< If true, the wind stress is distributed over the topmost
                               ///< Hmix_stress of fluid, and an added mixed layer viscosity or a
                               ///< physically based boundary layer turbulence parameterization is
                               ///< not needed for stability.
    bool   StokesMixing;       ///< If true, do Stokes drift mixing via the Lagrangian current
                               ///< (Eulerian plus Stokes drift).  False by default and set
                               ///< via STOKES_MIXING_COMBINED.
    double vel_underflow;      ///< Velocity components smaller than vel_underflow
                               ///< are set to 0 [L T-1 ~> m s-1].
    double CFL_trunc;          ///< Velocity components will be truncated when they
                               ///< are large enough that the corresponding CFL number
                               ///< exceeds this value [nondim].
    double CFL_report;         ///< The value of the CFL number that will cause the
                               ///< accelerations to be reported [nondim].  CFL_report
                               ///< will often equal CFL_trunc.
    bool   do_report_u;        ///< If true, the columns with large u-velocities are reported.
    bool   do_report_v;        ///< If true, the columns with large v-velocities are reported.
    bool   do_GLwork;          ///< If true, the KE source from GL90 vertical viscosity is needed for
                               ///< diagnostics.
    bool   do_du_dt_visc_gl90; ///< If true, the zonal acceleration from GL90 vertical
                               ///< viscosity is needed for diagnostics.
    bool   do_dv_dt_visc_gl90; ///< If true, the meridional acceleration from GL90 vertical
                               ///< viscosity is needed for diagnostics.
};

/// @brief Options used in the GL90 vertical viscosity parameterization.
/// Field-for-field mirror of the Fortran `bind(C)` type `GL90_CS_C` -- order and types must not
/// change.
struct GL90_CS_C {
    bool   use_GL90_N2;                ///< If true, use GL90 vertical viscosity coefficient that is
                                       ///< depth-independent; this corresponds to a kappa_GM that
                                       ///< scales as N^2 with depth.
    double kappa_gl90;                 ///< The scalar diffusivity used in the GL90 vertical
                                       ///< viscosity scheme [L2 H Z-1 T-1 ~> m2 s-1 or Pa s]
    bool   read_kappa_gl90;            ///< If true, read a file containing the spatially varying
                                       ///< kappa_gl90
    double alpha_gl90;                 ///< Coefficient used to compute a depth-independent GL90
                                       ///< vertical viscosity via Kv_gl90 = alpha_gl90 * f^2. Note
                                       ///< that the implied Kv_gl90 corresponds to a kappa_gl90
                                       ///< that scales as N^2 with depth. [H Z T ~> m2 s or kg s m-1]
    struct RealArray_C kappa_gl90_2d;  ///< 2D kappa_gl90 at h-points
                                       ///< [L2 H Z-1 T-1 ~> m2 s-1 or Pa s]
};

/// @brief Options used in the vertical viscosity coupling coefficient calculations.
/// Field-for-field mirror of the Fortran `bind(C)` type `coef_opts_CS_C` -- order and types must
/// not change.
struct coef_opts_CS_C {
    coupling_CS_C coupling_CS; ///< Options used in the calculation of the coupling coefficients
                               ///< across interfaces.
    double Hbbl;               ///< The static bottom boundary layer thickness [Z ~> m].
    double Hbbl_gl90;          ///< The static bottom boundary layer thickness used for GL90 [Z ~> m].
    bool   use_GL90_in_SSW;    ///< If true, use the GL90 parameterization in stacked shallow water
                               ///< mode (SSW). The calculation of the GL90 viscosity coefficient
                               ///< uses the fact that in SSW we simply have 1/N^2 = h/g^prime,
                               ///< where g^prime is the reduced gravity. This identity does not
                               ///< generalize to non-SSW setups.
    GL90_CS_C GL90_CS;         ///< Options used in the GL90 vertical viscosity parameterization.
    bool   harmonic_visc;      ///< If true, the harmonic mean thicknesses are used
                               ///< to calculate the viscous coupling between layers
                               ///< except near the bottom.  Otherwise the arithmetic
                               ///< mean thickness is used except near the bottom.
    double harm_BL_val;        ///< A scale to determine when water is in the boundary
                               ///< layers based solely on harmonic mean thicknesses
                               ///< for the purpose of determining the extent to which
                               ///< the thicknesses used in the viscosities are upwinded [nondim].
};

namespace MOM {
using amrex::Box;
using amrex::Array4;
using amrex::Real;

///
/// @brief Velocity components which exceed a threshold for physically reasonable values are truncated,
///        and the running sum of the number of trunctionas within the non-symmetric memory computational
///        domain is incremented.  Optionally, any column with excessive velocities may be sent
///        to a diagnostic reporting subroutine.
///
/// @param bxU0         The u-point box of the whole column, [IscB:IecB, jsc:jec, 1:nz]
/// @param bxV0         The v-point box of the whole column, [isc:iec, JscB:JecB, 1:nz]
/// @param u            Zonal velocity [L T-1 ~> m s-1]
/// @param v            Meridional velocity [L T-1 ~> m s-1]
/// @param dt           Time increment [T ~> s]
/// @param areaT        The area of an h-cell [L2 ~> m2]
/// @param IareaT       1/areaT [L-2 ~> m-2]
/// @param dx_Cv        The unblocked lengths of the v-faces of the h-cell [L ~> m]
/// @param dy_Cu        The unblocked lengths of the u-faces of the h-cell [L ~> m]
/// @param ke           The number of layers/levels in the vertical
/// @param Angstrom_H   A one-Angstrom thickness in the model thickness units [H ~> m or kg m-2]
/// @param m_s_to_L_T   Convert lateral velocities from m s-1 to L T-1 [L s T-1 m-1 ~> 1]
/// @param CS           Vertical viscosity control structure
/// @param h_u          The effective layer thickness at u-points [H ~> m or kg m-2]
/// @param h_v          The effective layer thickness at v-points [H ~> m or kg m-2]
/// @param u_old        The previous u-velocity [L T-1 ~> m s-1]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param v_old        The previous v-velocity [L T-1 ~> m s-1]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param vel_report_u The u-velocity to report [L T-1 ~> m s-1]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param vel_report_v The v-velocity to report [L T-1 ~> m s-1]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param dowrite_u    True for the u-point columns to report; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param dowrite_v    True for the v-point columns to report; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param ntrunc_sum   The number of velocity truncations, added to the value on entry [nondim]
///
void vertvisc_limit_vel(const Box& bxU0,
                        const Box& bxV0,
                        Array4<Real>       const& u,
                        Array4<Real>       const& v,
                        const Real dt,
                        Array4<const Real> const& areaT,
                        Array4<const Real> const& IareaT,
                        Array4<const Real> const& dx_Cv,
                        Array4<const Real> const& dy_Cu,
                        const int ke,
                        const Real Angstrom_H,
                        const Real m_s_to_L_T,
                        const vertvisc_opts_CS_C& CS,
                        Array4<const Real> const& h_u,
                        Array4<const Real> const& h_v,
                        Array4<Real>       const& u_old,
                        Array4<Real>       const& v_old,
                        Array4<Real>       const& vel_report_u,
                        Array4<Real>       const& vel_report_v,
                        Array4<int>        const& dowrite_u,
                        Array4<int>        const& dowrite_v,
                        int& ntrunc_sum);

///
/// @brief The implementation of vertvisc_remnant (the root of its call tree).
///
/// @param bxU0       The u-point box of the whole column, [IscB:IecB, jsc:jec, 1:nz]
/// @param bxV0       The v-point box of the whole column, [isc:iec, JscB:JecB, 1:nz]
/// @param Ray_u      The Rayleigh drag velocity to be applied to each layer at u-points
///                   [H T-1 ~> m s-1 or Pa s m-1]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param Ray_v      The Rayleigh drag velocity to be applied to each layer at v-points
///                   [H T-1 ~> m s-1 or Pa s m-1]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param visc_rem_u Fraction of a time-step's worth of a barotropic acceleration that a layer
///                   experiences after viscosity is applied in the zonal direction [nondim]
/// @param visc_rem_v Fraction of a time-step's worth of a barotropic acceleration that a layer
///                   experiences after viscosity is applied in the meridional direction [nondim]
/// @param dt         Time increment [T ~> s]
/// @param mask2dCu   0 for boundary points and 1 for ocean points on the u grid [nondim]
/// @param mask2dCv   0 for boundary points and 1 for ocean points on the v grid [nondim]
/// @param ke         The number of layers/levels in the vertical
/// @param a_u        The u-drag coefficient across an interface [H T-1 ~> m s-1 or Pa s m-1]
/// @param a_v        The v-drag coefficient across an interface [H T-1 ~> m s-1 or Pa s m-1]
/// @param h_u        The effective layer thickness at u-points [H ~> m or kg m-2]
/// @param h_v        The effective layer thickness at v-points [H ~> m or kg m-2]
///
void vertvisc_remnant_tr(const Box& bxU0,
                         const Box& bxV0,
                         Array4<const Real> const& Ray_u,
                         Array4<const Real> const& Ray_v,
                         Array4<Real>       const& visc_rem_u,
                         Array4<Real>       const& visc_rem_v,
                         const Real dt,
                         Array4<const Real> const& mask2dCu,
                         Array4<const Real> const& mask2dCv,
                         const int ke,
                         Array4<const Real> const& a_u,
                         Array4<const Real> const& a_v,
                         Array4<const Real> const& h_u,
                         Array4<const Real> const& h_v);

///
/// @brief Adds the nonlocal stress increments to ui^n and vi^n, and returns the Eulerian velocities and
///        the increments for the diagnostics.
///
/// @param bxU0     The u-point box of the whole column, [IscB:IecB, jsc:jec, 1:nz]
/// @param bxV0     The v-point box of the whole column, [isc:iec, JscB:JecB, 1:nz]
/// @param ui       Zonal velocity after vertvisc [L T-1 ~> m s-1]
/// @param vi       Meridional velocity after vertvisc [L T-1 ~> m s-1]
/// @param hbl_h    boundary layer depth [H ~> m]
/// @param taux     zonal wind stress [R L Z T-2 ~> Pa]
/// @param tauy     meridional wind stress [R L Z T-2 ~> Pa]
/// @param Us_x     3d zonal Stokes drift profile [L T-1 ~> m s-1]
/// @param Us_y     3d meridional Stokes drift profile [L T-1 ~> m s-1]
/// @param dt       Time increment [T ~> s]
/// @param Cemp_NL  empirical coefficient of non-local momentum mixing [nondim]
/// @param Irho0    Inverse of the mean density rescaled to [Z L-1 R-1 ~> m3 kg-1]
/// @param mask2dCu 0 for boundary points and 1 for ocean points on the u grid [nondim]
/// @param mask2dCv 0 for boundary points and 1 for ocean points on the v grid [nondim]
/// @param mask2dT  0 for land points and 1 for ocean points on the h-grid [nondim]
/// @param h_u      The effective layer thickness at u-points [H ~> m or kg m-2]
/// @param h_v      The effective layer thickness at v-points [H ~> m or kg m-2]
/// @param uE_u     zonal Eulerian u-pts [L T-1 ~> m s-1]
/// @param vE_v     merid Eulerian v-pts [L T-1 ~> m s-1]
/// @param uInc_u   zonal Eulerian u-pts [L T-1 ~> m s-1]
/// @param vInc_v   merid Eulerian v-pts [L T-1 ~> m s-1]
///
void vertfpmix_tr(const Box& bxU0,
                  const Box& bxV0,
                  Array4<Real>       const& ui,
                  Array4<Real>       const& vi,
                  Array4<const Real> const& hbl_h,
                  Array4<const Real> const& taux,
                  Array4<const Real> const& tauy,
                  Array4<const Real> const& Us_x,
                  Array4<const Real> const& Us_y,
                  const Real dt,
                  const Real Cemp_NL,
                  const Real Irho0,
                  Array4<const Real> const& mask2dCu,
                  Array4<const Real> const& mask2dCv,
                  Array4<const Real> const& mask2dT,
                  Array4<const Real> const& h_u,
                  Array4<const Real> const& h_v,
                  Array4<Real>       const& uE_u,
                  Array4<Real>       const& vE_v,
                  Array4<Real>       const& uInc_u,
                  Array4<Real>       const& vInc_v);

///
/// @brief The implementation of vertvisc_coef (the root of its call tree).
///
/// @param bxU0              The u-point box of the whole column, [IscB:IecB, jsc:jec, 1:nz]
/// @param bxV0              The v-point box of the whole column, [isc:iec, JscB:JecB, 1:nz]
/// @param u                 Zonal velocity [L T-1 ~> m s-1]
/// @param v                 Meridional velocity [L T-1 ~> m s-1]
/// @param h                 Layer thickness [H ~> m or kg m-2]
/// @param dz                Vertical distance across layers [Z ~> m]
/// @param frac_shelf_u      Fractional ice shelf coverage of u-cells, nondimensional from 0 to 1
///                          [nondim]. This is only associated if ice shelves are enabled, and is
///                          exactly 0 away from shelves or on land; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param frac_shelf_v      Fractional ice shelf coverage of v-cells, nondimensional from 0 to 1
///                          [nondim]. This is only associated if ice shelves are enabled, and is
///                          exactly 0 away from shelves or on land; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param Ustar_2d          The wind friction velocity, calculated using the Boussinesq reference
///                          density or the time-evolving surface density in non-Boussinesq mode
///                          [Z T-1 ~> m s-1]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param Kv_bbl_u          The bottom boundary layer viscosity at the u-points
///                          [H Z T-1 ~> m2 s-1 or Pa s]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param Kv_bbl_v          The bottom boundary layer viscosity at the v-points
///                          [H Z T-1 ~> m2 s-1 or Pa s]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param bbl_thick_u       The bottom boundary layer thickness at the u-points [Z ~> m];
///                          may be unallocated; pass a default-constructed Array4 when unallocated
/// @param bbl_thick_v       The bottom boundary layer thickness at the v-points [Z ~> m];
///                          may be unallocated; pass a default-constructed Array4 when unallocated
/// @param tbl_thick_shelf_u Thickness of the viscous top boundary layer under ice shelves at
///                          u-points [Z ~> m]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param tbl_thick_shelf_v Thickness of the viscous top boundary layer under ice shelves at
///                          v-points [Z ~> m]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param Kv_tbl_shelf_u    Viscosity in the viscous top boundary layer under ice shelves at
///                          u-points [H Z T-1 ~> m2 s-1 or Pa s]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param Kv_tbl_shelf_v    Viscosity in the viscous top boundary layer under ice shelves at
///                          v-points [H Z T-1 ~> m2 s-1 or Pa s]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param nkml_visc_u       The number of layers in the viscous surface mixed layer at u-points
///                          [nondim]. This is not an integer because there may be fractional
///                          layers, and it is stored in terms of layers, not depth, to facilitate
///                          the movement of the viscous boundary layer with the flow;
///                          may be unallocated; pass a default-constructed Array4 when unallocated
/// @param nkml_visc_v       The number of layers in the viscous surface mixed layer at v-points
///                          [nondim]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param Kv_shear          The shear-driven turbulent vertical viscosity at the interfaces between
///                          layers in tracer columns [H Z T-1 ~> m2 s-1 or Pa s]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param Kv_shear_Bu       The shear-driven turbulent vertical viscosity at the interfaces between
///                          layers in corner columns [H Z T-1 ~> m2 s-1 or Pa s]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param SpV_avg           The layer averaged in situ specific volume [R-1 ~> m3 kg-1];
///                          may be unallocated; pass a default-constructed Array4 when unallocated
/// @param dt                Time increment [T ~> s]
/// @param mask2dCu          0 for boundary points and 1 for ocean points on the u grid [nondim]
/// @param mask2dCv          0 for boundary points and 1 for ocean points on the v grid [nondim]
/// @param bathyT            Ocean bottom depth, referenced to Z_ref at tracer points [Z ~> m]
/// @param CoriolisBu        The Coriolis parameter at corner points [T-1 ~> s-1]
/// @param ke                The number of layers/levels in the vertical
/// @param Boussinesq        If true, make the Boussinesq approximation.
/// @param nkml              The number of layers at the top that should be treated as parts of a
///                          homogeneous region.
/// @param dZ_subroundoff    A thickness in height units that is so small that it can be added to a
///                          vertical distance of Angstrom_Z or 1e-17 m without changing it at the
///                          bit level [Z ~> m].  This is the height equivalent of H_subroundoff.
/// @param H_subroundoff     A thickness that is so small that it can be added to a thickness of
///                          Angstrom or larger without changing it at the bit level
///                          [H ~> m or kg m-2]. If Angstrom is 0 or exceedingly small, this is
///                          negligible compared to 1e-17 m.
/// @param H_to_m            A constant that translates distances in the units of thickness to m
///                          [m H-1 ~> 1 or m3 kg-1].
/// @param m_to_H            A constant that translates distances in m to the units of thickness
///                          [H m-1 ~> 1 or kg m-3].
/// @param RZ_to_H           A constant that translates mass per unit area units to thickness
///                          units [H R-1 Z-1 ~> m3 kg-2 or 1].
/// @param Z_to_H            A constant that translates depth units to thickness units depth
///                          [H Z-1 ~> 1 or kg m-3].
/// @param T_to_s            A constant that translates the units of time to seconds [s T-1 ~> 1]
/// @param g_prime           The reduced gravity at each interface [L2 Z-1 T-2 ~> m s-2];
///                          may be unallocated; pass a default-constructed Array4 when unallocated
/// @param CS                Vertical viscosity control structure
/// @param kappa_gl90_2d     2D kappa_gl90 at h-points [L2 H Z-1 T-1 ~> m2 s-1 or Pa s];
///                          may be unallocated; pass a default-constructed Array4 when unallocated
/// @param a_u               The u-drag coefficient across an interface [H T-1 ~> m s-1 or Pa s m-1]
/// @param a_u_gl90          The u-drag coefficient associated with GL90 across an interface
///                          [H T-1 ~> m s-1 or Pa s m-1]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param h_u               The effective layer thickness at u-points [H ~> m or kg m-2].
/// @param a_v               The v-drag coefficient across an interface [H T-1 ~> m s-1 or Pa s m-1]
/// @param a_v_gl90          The v-drag coefficient associated with GL90 across an interface
///                          [H T-1 ~> m s-1 or Pa s m-1]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param h_v               The effective layer thickness at v-points [H ~> m or kg m-2].
/// @param a1_shelf_u        The u-momentum coupling coefficient under ice shelves
///                          [H T-1 ~> m s-1 or Pa s m-1]. Retained to determine stress under
///                          shelves; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param a1_shelf_v        The v-momentum coupling coefficient under ice shelves
///                          [H T-1 ~> m s-1 or Pa s m-1]. Retained to determine stress under
///                          shelves; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param OBC               Open boundary condition structure
/// @param kdgl90_struct     Vertical structure function used in GL90 diffusivity [nondim];
///                          may be unallocated; pass a default-constructed Array4 when unallocated
/// @param hML_u             Diagnostic of the mixed layer depth at u points [Z ~> m];
///                          may be unallocated; pass a default-constructed Array4 when unallocated
/// @param hML_v             Diagnostic of the mixed layer depth at v points [Z ~> m];
///                          may be unallocated; pass a default-constructed Array4 when unallocated
/// @param Kv_u              Total vertical viscosity at u-points in thickness-based units
///                          [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1];
///                          may be unallocated; pass a default-constructed Array4 when unallocated
/// @param Kv_v              Total vertical viscosity at v-points in thickness-based units
///                          [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1];
///                          may be unallocated; pass a default-constructed Array4 when unallocated
/// @param Kv_gl90_u         GL90 vertical viscosity at u-points in thickness-based units
///                          [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1];
///                          may be unallocated; pass a default-constructed Array4 when unallocated
/// @param Kv_gl90_v         GL90 vertical viscosity at v-points in thickness-based units
///                          [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1];
///                          may be unallocated; pass a default-constructed Array4 when unallocated
///
void vertvisc_coef_tr(const Box& bxU0,
                      const Box& bxV0,
                      Array4<const Real> const& u,
                      Array4<const Real> const& v,
                      Array4<const Real> const& h,
                      Array4<const Real> const& dz,
                      Array4<const Real> const& frac_shelf_u,
                      Array4<const Real> const& frac_shelf_v,
                      Array4<const Real> const& Ustar_2d,
                      Array4<const Real> const& Kv_bbl_u,
                      Array4<const Real> const& Kv_bbl_v,
                      Array4<const Real> const& bbl_thick_u,
                      Array4<const Real> const& bbl_thick_v,
                      Array4<const Real> const& tbl_thick_shelf_u,
                      Array4<const Real> const& tbl_thick_shelf_v,
                      Array4<const Real> const& Kv_tbl_shelf_u,
                      Array4<const Real> const& Kv_tbl_shelf_v,
                      Array4<const Real> const& nkml_visc_u,
                      Array4<const Real> const& nkml_visc_v,
                      Array4<const Real> const& Kv_shear,
                      Array4<const Real> const& Kv_shear_Bu,
                      Array4<const Real> const& SpV_avg,
                      const Real dt,
                      Array4<const Real> const& mask2dCu,
                      Array4<const Real> const& mask2dCv,
                      Array4<const Real> const& bathyT,
                      Array4<const Real> const& CoriolisBu,
                      const int ke,
                      const bool Boussinesq,
                      const int nkml,
                      const Real dZ_subroundoff,
                      const Real H_subroundoff,
                      const Real H_to_m,
                      const Real m_to_H,
                      const Real RZ_to_H,
                      const Real Z_to_H,
                      const Real T_to_s,
                      Array4<const Real> const& g_prime,
                      const coef_opts_CS_C& CS,
                      Array4<const Real> const& kappa_gl90_2d,
                      Array4<Real>       const& a_u,
                      Array4<Real>       const& a_u_gl90,
                      Array4<Real>       const& h_u,
                      Array4<Real>       const& a_v,
                      Array4<Real>       const& a_v_gl90,
                      Array4<Real>       const& h_v,
                      Array4<Real>       const& a1_shelf_u,
                      Array4<Real>       const& a1_shelf_v,
                      OceanOBC* OBC,
                      Array4<const Real> const& kdgl90_struct,
                      Array4<Real>       const& hML_u,
                      Array4<Real>       const& hML_v,
                      Array4<Real>       const& Kv_u,
                      Array4<Real>       const& Kv_v,
                      Array4<Real>       const& Kv_gl90_u,
                      Array4<Real>       const& Kv_gl90_v);

///
/// @brief The implementation of vertvisc (the root of its call tree).
///
/// @param bxU0            The u-point box of the whole column, [IscB:IecB, jsc:jec, 1:nz]
/// @param bxV0            The v-point box of the whole column, [isc:iec, JscB:JecB, 1:nz]
/// @param u               Zonal velocity [L T-1 ~> m s-1]
/// @param v               Meridional velocity [L T-1 ~> m s-1]
/// @param h               Layer thickness [H ~> m or kg m-2]
/// @param taux            zonal wind stress [R L Z T-2 ~> Pa]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param tauy            meridional wind stress [R L Z T-2 ~> Pa]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param Ray_u           The Rayleigh drag velocity to be applied to each layer
///                        at u-points [H T-1 ~> m s-1 or Pa s m-1]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param Ray_v           The Rayleigh drag velocity to be applied to each layer
///                        at v-points [H T-1 ~> m s-1 or Pa s m-1]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param taux_shelf      The zonal stresses on the ocean under shelves
///                        [R Z L T-2 ~> Pa]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param tauy_shelf      The meridional stresses on the ocean under shelves
///                        [R Z L T-2 ~> Pa]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param dt              Time increment [T ~> s]
/// @param OBC             Open boundary condition structure
/// @param du_dt_visc      Zonal acceleration due to vertical viscosity
///                        [L T-2 ~> m s-2]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param dv_dt_visc      Meridional acceleration due to vertical viscosity
///                        [L T-2 ~> m s-2]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param du_dt_visc_gl90 Zonal acceleration due to GL90 vertical viscosity
///                        (is included in du_dt_visc) [L T-2 ~> m s-2]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param dv_dt_visc_gl90 Meridional acceleration due to GL90 vertical viscosity
///                        (is included in dv_dt_visc) [L T-2 ~> m s-2]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param du_dt_str       Zonal acceleration due to the surface stress (included
///                        in du_dt_visc) [L T-2 ~> m s-2]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param dv_dt_str       Meridional acceleration due to the surface stress
///                        (included in dv_dt_visc) [L T-2 ~> m s-2]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param mask2dCu        0 for boundary points and 1 for ocean points on the u grid [nondim].
/// @param mask2dCv        0 for boundary points and 1 for ocean points on the v grid [nondim].
/// @param areaCu          The areas of the u-grid cells [L2 ~> m2].
/// @param areaCv          The areas of the v-grid cells [L2 ~> m2].
/// @param areaT           The area of an h-cell [L2 ~> m2].
/// @param IareaT          1/areaT [L-2 ~> m-2].
/// @param dx_Cv           The unblocked lengths of the v-faces of the h-cell [L ~> m].
/// @param dy_Cu           The unblocked lengths of the u-faces of the h-cell [L ~> m].
/// @param ke              The number of layers/levels in the vertical
/// @param H_subroundoff   A thickness that is so small that it can be added to a thickness of
///                        Angstrom or larger without changing it at the bit level
///                        [H ~> m or kg m-2]. If Angstrom is 0 or exceedingly small, this is
///                        negligible compared to 1e-17 m.
/// @param Angstrom_H      A one-Angstrom thickness in the model thickness units
///                        [H ~> m or kg m-2].
/// @param H_to_RZ         A constant that translates thickness units to the units of
///                        mass per unit area [R Z H-1 ~> kg m-3 or 1].
/// @param m_s_to_L_T      Convert lateral velocities from m s-1 to L T-1
///                        [L s T-1 m-1 ~> 1]
/// @param CS              Vertical viscosity control structure
/// @param a_u             The u-drag coefficient across an interface
///                        [H T-1 ~> m s-1 or Pa s m-1]
/// @param a_u_gl90        The u-drag coefficient associated with GL90 across an interface
///                        [H T-1 ~> m s-1 or Pa s m-1]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param h_u             The effective layer thickness at u-points [H ~> m or kg m-2].
/// @param a_v             The v-drag coefficient across an interface
///                        [H T-1 ~> m s-1 or Pa s m-1]
/// @param a_v_gl90        The v-drag coefficient associated with GL90 across an interface
///                        [H T-1 ~> m s-1 or Pa s m-1]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param h_v             The effective layer thickness at v-points [H ~> m or kg m-2].
/// @param a1_shelf_u      The u-momentum coupling coefficient under ice shelves
///                        [H T-1 ~> m s-1 or Pa s m-1]. Retained to determine stress
///                        under shelves; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param a1_shelf_v      The v-momentum coupling coefficient under ice shelves
///                        [H T-1 ~> m s-1 or Pa s m-1]. Retained to determine stress
///                        under shelves; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param taux_bot        Zonal bottom stress from ocean to rock [R L Z T-2 ~> Pa];
///                        may be unallocated; pass a default-constructed Array4 when unallocated
/// @param tauy_bot        Meridional bottom stress from ocean to rock [R L Z T-2 ~> Pa];
///                        may be unallocated; pass a default-constructed Array4 when unallocated
/// @param KE_u            The area integral of a KE term in a layer at u-points
///                        [H L4 T-3 ~> m5 s-3 or kg m2 s-3]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param KE_v            The area integral of a KE term in a layer at v-points
///                        [H L4 T-3 ~> m5 s-3 or kg m2 s-3]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param u_old           The previous u-velocity [L T-1 ~> m s-1]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param v_old           The previous v-velocity [L T-1 ~> m s-1]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param vel_report_u    The u-velocity to report [L T-1 ~> m s-1]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param vel_report_v    The v-velocity to report [L T-1 ~> m s-1]; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param dowrite_u       True for the u-point columns to report; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param dowrite_v       True for the v-point columns to report; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param ntrunc_sum      The number of velocity truncations, added to the value on entry [nondim]
/// @param fpmix           fpmix along Eulerian shear
/// @param Us_x            3d zonal Stokes drift profile [L T-1 ~> m s-1]
///                        Horizontal -> U points
///                        Vertical -> Mid-points; may be unallocated; pass a default-constructed Array4 when unallocated
/// @param Us_y            3d meridional Stokes drift profile [L T-1 ~> m s-1]
///                        Horizontal -> V points
///                        Vertical -> Mid-points; may be unallocated; pass a default-constructed Array4 when unallocated
///
void vertvisc_tr(const Box& bxU0,
                 const Box& bxV0,
                 Array4<Real>       const& u,
                 Array4<Real>       const& v,
                 Array4<const Real> const& h,
                 Array4<const Real> const& taux,
                 Array4<const Real> const& tauy,
                 Array4<const Real> const& Ray_u,
                 Array4<const Real> const& Ray_v,
                 Array4<Real>       const& taux_shelf,
                 Array4<Real>       const& tauy_shelf,
                 const Real dt,
                 OceanOBC* OBC,
                 Array4<Real>       const& du_dt_visc,
                 Array4<Real>       const& dv_dt_visc,
                 Array4<Real>       const& du_dt_visc_gl90,
                 Array4<Real>       const& dv_dt_visc_gl90,
                 Array4<Real>       const& du_dt_str,
                 Array4<Real>       const& dv_dt_str,
                 Array4<const Real> const& mask2dCu,
                 Array4<const Real> const& mask2dCv,
                 Array4<const Real> const& areaCu,
                 Array4<const Real> const& areaCv,
                 Array4<const Real> const& areaT,
                 Array4<const Real> const& IareaT,
                 Array4<const Real> const& dx_Cv,
                 Array4<const Real> const& dy_Cu,
                 const int ke,
                 const Real H_subroundoff,
                 const Real Angstrom_H,
                 const Real H_to_RZ,
                 const Real m_s_to_L_T,
                 const vertvisc_opts_CS_C& CS,
                 Array4<const Real> const& a_u,
                 Array4<const Real> const& a_u_gl90,
                 Array4<const Real> const& h_u,
                 Array4<const Real> const& a_v,
                 Array4<const Real> const& a_v_gl90,
                 Array4<const Real> const& h_v,
                 Array4<const Real> const& a1_shelf_u,
                 Array4<const Real> const& a1_shelf_v,
                 Array4<Real>       const& taux_bot,
                 Array4<Real>       const& tauy_bot,
                 Array4<Real>       const& KE_u,
                 Array4<Real>       const& KE_v,
                 Array4<Real>       const& u_old,
                 Array4<Real>       const& v_old,
                 Array4<Real>       const& vel_report_u,
                 Array4<Real>       const& vel_report_v,
                 Array4<int>        const& dowrite_u,
                 Array4<int>        const& dowrite_v,
                 int& ntrunc_sum,
                 const bool fpmix,
                 Array4<const Real> const& Us_x,
                 Array4<const Real> const& Us_y);

} // namespace MOM
