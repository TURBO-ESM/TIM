// SKILLS: 0.3.2
#pragma once
///
/// @file mom_coriolisadv.hpp
/// @brief Box-level AMReX kernel declarations for the MOM6 Coriolis and momentum advection terms.
///
#include <AMReX_Box.H>
#include <AMReX_Array4.H>

#include "mom_coriolisadv_kernel.hpp"

struct OceanOBC;          // Undefined at the moment
struct WaveParametersCS;  // Undefined at the moment

/// @brief Options selecting the discretization of the Coriolis terms. Field-for-field
/// mirror of the Fortran `bind(C)` type `Coriolis_scheme_CS_C` -- order and types must
/// not change.
struct Coriolis_scheme_CS_C {
    int  Coriolis_Scheme;  ///< Selects the discretization for the Coriolis terms.
    bool Coriolis_En_Dis;  ///< If true, use the energy-dissipating biased Coriolis scheme.
};

/// @brief Options used in the Coriolis and momentum advection calculations.
/// Field-for-field mirror of the Fortran `bind(C)` type `CorAdv_opts_CS_C` -- order and
/// types must not change.
struct CorAdv_opts_CS_C {
    Coriolis_scheme_CS_C Coriolis_scheme_CS; ///< Options selecting the discretization of the
                                             ///< Coriolis terms.
    int    KE_Scheme;            ///< KE_SCHEME selects the discretization for
                                 ///< the kinetic energy. Valid values are:
                                 ///<  KE_ARAKAWA, KE_SIMPLE_GUDONOV, KE_GUDONOV
    bool   KE_use_limiter;       ///< If true, use the Koren limiter for KE_UP3 scheme
    int    PV_Adv_Scheme;        ///< PV_ADV_SCHEME selects the discretization for PV advection
                                 ///< Valid values are:
                                 ///< - PV_ADV_CENTERED - centered (aka Sadourny, 75)
                                 ///< - PV_ADV_UPWIND1  - upwind, first order
    int    nkblock;              ///< The k block size used in Coriolis/advection calculations [nondim].
    double F_eff_max_blend;      ///< The factor by which the maximum effective Coriolis
                                 ///< acceleration from any point can be increased when
                                 ///< blending different discretizations with the
                                 ///< ARAKAWA_LAMB_BLEND Coriolis scheme [nondim].
                                 ///< This must be greater than 2.0, and is 4.0 by default.
    double wt_lin_blend;         ///< A weighting value beyond which the blending between
                                 ///< Sadourny and Arakawa & Hsu goes linearly to 0 [nondim].
                                 ///< This must be between 1 and 1e-15, often 1/8.
    bool   no_slip;              ///< If true, no slip boundary conditions are used.
                                 ///< Otherwise free slip boundary conditions are assumed.
                                 ///< The implementation of the free slip boundary
                                 ///< conditions on a C-grid is much cleaner than the
                                 ///< no slip boundary conditions. The use of free slip
                                 ///< b.c.s is strongly encouraged. The no slip b.c.s
                                 ///< are not implemented with the biharmonic viscosity.
    bool   bound_Coriolis;       ///< If true, the Coriolis terms at u points are
                                 ///< bounded by the four estimates of (f+rv)v from the
                                 ///< four neighboring v points, and similarly at v
                                 ///< points.  This option would have no effect on the
                                 ///< SADOURNY75_ENERGY scheme if it were possible to
                                 ///< use centered difference thickness fluxes.
    bool   weno_velocity_smooth; ///< If true, use velocity to compute the smoothness
                                 ///< indicator for WENO
    bool   do_CAS_diag;          ///< If true, the Stokes contribution to CAu or CAv is needed
                                 ///< for diagnostics.
    bool   do_RV_diag;           ///< If true, the relative vorticity is needed for diagnostics.
    bool   do_PV_diag;           ///< If true, the potential vorticity is needed for diagnostics.
};

namespace MOM {
using amrex::Box;
using amrex::Array4;

/// @name Enumeration values for Coriolis_Scheme
/// @{
constexpr int SADOURNY75_ENERGY   = 1; ///< Sadourny (1975) energy-conserving scheme ("SADOURNY75_ENERGY")
constexpr int ARAKAWA_HSU90       = 2; ///< Arakawa and Hsu (1990) energy-enstrophy scheme ("ARAKAWA_HSU90")
constexpr int ROBUST_ENSTRO       = 3; ///< Enstrophy-conserving scheme robust to vanishing layers ("ROBUST_ENSTRO")
constexpr int SADOURNY75_ENSTRO   = 4; ///< Sadourny (1975) enstrophy-conserving scheme ("SADOURNY75_ENSTRO")
constexpr int ARAKAWA_LAMB81      = 5; ///< Arakawa and Lamb (1981) energy-enstrophy scheme ("ARAKAWA_LAMB81")
constexpr int AL_BLEND            = 6; ///< Blend of Arakawa-Lamb (1981) and Arakawa-Hsu (1990) ("ARAKAWA_LAMB_BLEND")
constexpr int wenovi7th_PV_ENSTRO = 7; ///< 7th-order WENO PV-enstrophy scheme ("WENOVI7TH_PV_ENSTRO")
constexpr int wenovi5th_PV_ENSTRO = 8; ///< 5th-order WENO PV-enstrophy scheme ("WENOVI5TH_PV_ENSTRO")
constexpr int wenovi3rd_PV_ENSTRO = 9; ///< 3rd-order WENO PV-enstrophy scheme ("WENOVI3RD_PV_ENSTRO")
/// @}

/// @name Enumeration values for PV_Adv_Scheme
/// @{
constexpr int PV_ADV_CENTERED = 21; ///< Centered PV advection ("PV_ADV_CENTERED")
constexpr int PV_ADV_UPWIND1  = 22; ///< First-order upwind PV advection ("PV_ADV_UPWIND1")
/// @}

/// @name Enumeration values for KE_Scheme
/// @{
constexpr int KE_ARAKAWA        = 10; ///< Arakawa and Lamb (1977) kinetic energy discretization ("KE_ARAKAWA")
constexpr int KE_SIMPLE_GUDONOV = 11; ///< Simple Godunov kinetic energy discretization ("KE_SIMPLE_GUDONOV")
constexpr int KE_GUDONOV        = 12; ///< Godunov kinetic energy discretization ("KE_GUDONOV")
constexpr int KE_UP3            = 13; ///< Third-order upwind kinetic energy discretization ("KE_UP3")
/// @}

///
/// @brief Calculates the acceleration due to the gradient of kinetic energy over one iteration tile.
///
/// @param bxQ            The B-grid-index iteration box of this tile, [IscB:IecB, JscB:JecB, ksc:kec]
/// @param bxU            The u-point iteration box of this tile, [IscB:IecB, jsc:jec, ksc:kec]
/// @param bxV            The v-point iteration box of this tile, [isc:iec, JscB:JecB, ksc:kec]
/// @param u              Zonal velocity [L T-1 ~> m s-1]
/// @param v              Meridional velocity [L T-1 ~> m s-1]
/// @param KE             Kinetic energy per unit mass [L2 T-2 ~> m2 s-2]
/// @param KEx            Zonal acceleration due to kinetic energy gradient [L T-2 ~> m s-2]
/// @param KEy            Meridional acceleration due to kinetic energy gradient [L T-2 ~> m s-2]
/// @param areaCu         The areas of the u-grid cells [L2 ~> m2].
/// @param areaCv         The areas of the v-grid cells [L2 ~> m2].
/// @param IareaT         1/areaT [L-2 ~> m-2].
/// @param mask2dCu       0 for boundary points and 1 for ocean points on the u grid [nondim].
/// @param mask2dCv       0 for boundary points and 1 for ocean points on the v grid [nondim].
/// @param IdxCu_OBCmask  1/dxCu or 0 at boundary or OBC points [L-1 ~> m-1].
/// @param IdyCv_OBCmask  1/dxCv or 0 at boundary or OBC points [L-1 ~> m-1].
/// @param KE_Scheme      Selects the discretization for the kinetic energy: KE_ARAKAWA,
///                       KE_SIMPLE_GUDONOV, KE_GUDONOV or KE_UP3
/// @param KE_use_limiter If true, use the Koren limiter for the KE_UP3 scheme
///
void gradke(const Box& bxQ,
            const Box& bxU,
            const Box& bxV,
            Array4<const Real> const& u,
            Array4<const Real> const& v,
            Array4<Real>       const& KE,
            Array4<Real>       const& KEx,
            Array4<Real>       const& KEy,
            Array4<const Real> const& areaCu,
            Array4<const Real> const& areaCv,
            Array4<const Real> const& IareaT,
            Array4<const Real> const& mask2dCu,
            Array4<const Real> const& mask2dCv,
            Array4<const Real> const& IdxCu_OBCmask,
            Array4<const Real> const& IdyCv_OBCmask,
            const int KE_Scheme,
            const bool KE_use_limiter);

///
/// @brief Calculates the potential vorticity and the related thickness and transport fields
///        used by the Coriolis and momentum advection terms over one iteration tile.
///
/// @param bxH               The h-point iteration box of this tile, [isc:iec, jsc:jec, ksc:kec]
/// @param bxQ               The B-grid-index iteration box of this tile, [IscB:IecB, JscB:JecB, ksc:kec]
/// @param bxU               The u-point iteration box of this tile, [IscB:IecB, jsc:jec, ksc:kec]
/// @param bxV               The v-point iteration box of this tile, [isc:iec, JscB:JecB, ksc:kec]
/// @param bxQs              The scheme-dependent vorticity-point box, [Is_q:Ie_q, Js_q:Je_q, ksc:kec]
/// @param u                 Zonal velocity [L T-1 ~> m s-1]
/// @param v                 Meridional velocity [L T-1 ~> m s-1]
/// @param h                 Layer thickness [H ~> m or kg m-2]
/// @param OBC               Open boundary control structure
/// @param rv_x_v_associated True if the rv_x_v acceleration diagnostic is in use
/// @param rv_x_u_associated True if the rv_x_u acceleration diagnostic is in use
/// @param por_face_areaU    fractional open area of U-faces [nondim]
/// @param por_face_areaV    fractional open area of V-faces [nondim]
/// @param Waves             An optional pointer to Stokes drift CS
/// @param dyCv              dyCv is delta y at v points [L ~> m].
/// @param dxCu              dxCu is delta x at u points [L ~> m].
/// @param dy_Cu             The unblocked lengths of the u-faces of the h-cell [L ~> m].
/// @param dx_Cv             The unblocked lengths of the v-faces of the h-cell [L ~> m].
/// @param dyBu              dyBu is delta y at q points [L ~> m].
/// @param dxBu              dxBu is delta x at q points [L ~> m].
/// @param mask2dBu          0 for boundary points and 1 for ocean points on the q grid [nondim].
/// @param IareaBu           IareaBu = 1/areaBu [L-2 ~> m-2].
/// @param CoriolisBu        The Coriolis parameter at corner points [T-1 ~> s-1].
/// @param Area_h            The ocean area at h points [L2 ~> m2]
/// @param Area_q            The sum of the ocean areas at the 4 adjacent thickness points [L2 ~> m2]
/// @param Stokes_VF         If true, include the Stokes drift
/// @param use_weno          True if using one of the WENO schemes
/// @param vol_neglect       A negligible volume [H L2 ~> m3 or kg]
/// @param area_neglect      A negligible area [L2 ~> m2]
/// @param q                 Layer potential vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
/// @param Ih_q              Inverse of thickness interpolated to q points [H-1 ~> m-1 or m2 kg-1]
/// @param abs_vort          Absolute vorticity at q-points [T-1 ~> s-1]
/// @param h_q               Thickness interpolated to q points [H ~> m or kg m-2]
/// @param qS                Layer Stokes vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
/// @param q2                Relative vorticity over thickness [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
/// @param uh_center         Transport based on arithmetic mean h at u-points
///                          [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param vh_center         Transport based on arithmetic mean h at v-points
///                          [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param RV                Diagnostic relative vorticity [T-1 ~> s-1]
/// @param PV                Diagnostic potential vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
/// @param CS                Options used in the Coriolis and momentum advection calculations
///
void coradv_setup(const Box& bxH,
                  const Box& bxQ,
                  const Box& bxU,
                  const Box& bxV,
                  const Box& bxQs,
                  Array4<const Real> const& u,
                  Array4<const Real> const& v,
                  Array4<const Real> const& h,
                  OceanOBC* OBC,
                  const bool rv_x_v_associated,
                  const bool rv_x_u_associated,
                  Array4<const Real> const& por_face_areaU,
                  Array4<const Real> const& por_face_areaV,
                  WaveParametersCS* Waves,
                  Array4<const Real> const& dyCv,
                  Array4<const Real> const& dxCu,
                  Array4<const Real> const& dy_Cu,
                  Array4<const Real> const& dx_Cv,
                  Array4<const Real> const& dyBu,
                  Array4<const Real> const& dxBu,
                  Array4<const Real> const& mask2dBu,
                  Array4<const Real> const& IareaBu,
                  Array4<const Real> const& CoriolisBu,
                  Array4<const Real> const& Area_h,
                  Array4<const Real> const& Area_q,
                  const bool Stokes_VF,
                  const bool use_weno,
                  const Real vol_neglect,
                  const Real area_neglect,
                  Array4<Real>       const& q,
                  Array4<Real>       const& Ih_q,
                  Array4<Real>       const& abs_vort,
                  Array4<Real>       const& h_q,
                  Array4<Real>       const& qS,
                  Array4<Real>       const& q2,
                  Array4<Real>       const& uh_center,
                  Array4<Real>       const& vh_center,
                  Array4<Real>       const& RV,
                  Array4<Real>       const& PV,
                  const CorAdv_opts_CS_C& CS);

///
/// @brief Adds the Coriolis bounds and the kinetic energy gradient to the Coriolis and momentum
///        advection accelerations, and calculates the Stokes and acceleration diagnostics, over
///        one iteration tile.
///
/// @param bxU       The u-point iteration box of this tile, [IscB:IecB, jsc:jec, ksc:kec]
/// @param bxV       The v-point iteration box of this tile, [isc:iec, JscB:JecB, ksc:kec]
/// @param u         Zonal velocity [L T-1 ~> m s-1]
/// @param v         Meridional velocity [L T-1 ~> m s-1]
/// @param uh        Zonal transport u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param vh        Meridional transport v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param abs_vort  Absolute vorticity at q-points [T-1 ~> s-1]
/// @param qS        Layer Stokes vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
/// @param q2        Relative vorticity over thickness [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
/// @param KEx       Zonal gradient of kinetic energy per unit mass [L T-2 ~> m s-2]
/// @param KEy       Meridional gradient of kinetic energy per unit mass [L T-2 ~> m s-2]
/// @param Stokes_VF If true, include the Stokes drift
/// @param IdxCu     1/dxCu [L-1 ~> m-1].
/// @param IdyCv     1/dyCv [L-1 ~> m-1].
/// @param CAu       Zonal acceleration due to Coriolis and momentum advection [L T-2 ~> m s-2].
/// @param CAv       Meridional acceleration due to Coriolis and momentum advection [L T-2 ~> m s-2].
/// @param CAuS      Stokes contribution to CAu [L T-2 ~> m s-2]
/// @param CAvS      Stokes contribution to CAv [L T-2 ~> m s-2]
/// @param gradKEu   gradKEu = - d/dx(u2) [L T-2 ~> m s-2]; may be absent (.p == nullptr)
/// @param gradKEv   gradKEv = - d/dy(u2) [L T-2 ~> m s-2]; may be absent (.p == nullptr)
/// @param rv_x_v    rv_x_v = rv * v at u [L T-2 ~> m s-2]; may be absent (.p == nullptr)
/// @param rv_x_u    rv_x_u = rv * u at v [L T-2 ~> m s-2]; may be absent (.p == nullptr)
/// @param CS        Options used in the Coriolis and momentum advection calculations
///
void coradv_common_terms(const Box& bxU,
                         const Box& bxV,
                         Array4<const Real> const& u,
                         Array4<const Real> const& v,
                         Array4<const Real> const& uh,
                         Array4<const Real> const& vh,
                         Array4<const Real> const& abs_vort,
                         Array4<const Real> const& qS,
                         Array4<const Real> const& q2,
                         Array4<const Real> const& KEx,
                         Array4<const Real> const& KEy,
                         const bool Stokes_VF,
                         Array4<const Real> const& IdxCu,
                         Array4<const Real> const& IdyCv,
                         Array4<Real>       const& CAu,
                         Array4<Real>       const& CAv,
                         Array4<Real>       const& CAuS,
                         Array4<Real>       const& CAvS,
                         Array4<Real>       const& gradKEu,
                         Array4<Real>       const& gradKEv,
                         Array4<Real>       const& rv_x_v,
                         Array4<Real>       const& rv_x_u,
                         const CorAdv_opts_CS_C& CS);

///
/// @brief Calculates the Coriolis and momentum advection accelerations with the Arakawa & Hsu
///        (1990), Arakawa & Lamb (1981) or ARAKAWA_LAMB_BLEND schemes over one iteration tile.
///
/// @param bxH          The h-point iteration box of this tile, [isc:iec, jsc:jec, ksc:kec]
/// @param bxQ          The B-grid-index iteration box of this tile, [IscB:IecB, JscB:JecB, ksc:kec]
/// @param bxU          The u-point iteration box of this tile, [IscB:IecB, jsc:jec, ksc:kec]
/// @param bxV          The v-point iteration box of this tile, [isc:iec, JscB:JecB, ksc:kec]
/// @param uh           Zonal transport u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param vh           Meridional transport v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param q            Layer potential vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
/// @param Ih_q         Inverse of thickness interpolated to q points [H-1 ~> m-1 or m2 kg-1]
/// @param Fe_m2        ARAKAWA_LAMB_BLEND blending parameter [nondim]
/// @param rat_lin      ARAKAWA_LAMB_BLEND blending parameter [nondim]
/// @param IdxCu        1/dxCu [L-1 ~> m-1].
/// @param IdyCv        1/dyCv [L-1 ~> m-1].
/// @param CAu          Zonal acceleration due to Coriolis and momentum advection [L T-2 ~> m s-2].
/// @param CAv          Meridional acceleration due to Coriolis and momentum advection [L T-2 ~> m s-2].
/// @param wt_lin_blend A weighting value beyond which the blending between Sadourny and
///                     Arakawa & Hsu goes linearly to 0 [nondim].
/// @param CS           Options selecting the Coriolis discretization
///
void coradv_arakawa(const Box& bxH,
                    const Box& bxQ,
                    const Box& bxU,
                    const Box& bxV,
                    Array4<const Real> const& uh,
                    Array4<const Real> const& vh,
                    Array4<const Real> const& q,
                    Array4<const Real> const& Ih_q,
                    const Real Fe_m2,
                    const Real rat_lin,
                    Array4<const Real> const& IdxCu,
                    Array4<const Real> const& IdyCv,
                    Array4<Real>       const& CAu,
                    Array4<Real>       const& CAv,
                    const Real wt_lin_blend,
                    const Coriolis_scheme_CS_C& CS);

///
/// @brief Calculates the Coriolis and momentum advection accelerations with an enstrophy
///        conserving scheme robust to vanishing layers over one iteration tile.
///
/// @param bxU           The u-point iteration box of this tile, [IscB:IecB, jsc:jec, ksc:kec]
/// @param bxV           The v-point iteration box of this tile, [isc:iec, JscB:JecB, ksc:kec]
/// @param u             Zonal velocity [L T-1 ~> m s-1]
/// @param v             Meridional velocity [L T-1 ~> m s-1]
/// @param h             Layer thickness [H ~> m or kg m-2]
/// @param uh            Zonal transport u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param vh            Meridional transport v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param abs_vort      Absolute vorticity at q-points [T-1 ~> s-1]
/// @param eps_vel       A tiny, positive velocity [L T-1 ~> m s-1]
/// @param h_tiny        A very small thickness [H ~> m or kg m-2]
/// @param IdxCv         1/dxCv [L-1 ~> m-1].
/// @param IdxCu         1/dxCu [L-1 ~> m-1].
/// @param IdyCu         1/dyCu [L-1 ~> m-1].
/// @param IdyCv         1/dyCv [L-1 ~> m-1].
/// @param CAu           Zonal acceleration due to Coriolis and momentum advection [L T-2 ~> m s-2].
/// @param CAv           Meridional acceleration due to Coriolis and momentum advection [L T-2 ~> m s-2].
/// @param PV_Adv_Scheme Selects the discretization for PV advection: PV_ADV_CENTERED or
///                      PV_ADV_UPWIND1
///
void coradv_robust_enstro(const Box& bxU,
                          const Box& bxV,
                          Array4<const Real> const& u,
                          Array4<const Real> const& v,
                          Array4<const Real> const& h,
                          Array4<const Real> const& uh,
                          Array4<const Real> const& vh,
                          Array4<const Real> const& abs_vort,
                          const Real eps_vel,
                          const Real h_tiny,
                          Array4<const Real> const& IdxCv,
                          Array4<const Real> const& IdxCu,
                          Array4<const Real> const& IdyCu,
                          Array4<const Real> const& IdyCv,
                          Array4<Real>       const& CAu,
                          Array4<Real>       const& CAv,
                          const int PV_Adv_Scheme);

///
/// @brief Calculates the Coriolis and momentum advection accelerations with the Sadourny (1975)
///        energy or enstrophy conserving schemes over one iteration tile.
///
/// @param bxH       The h-point iteration box of this tile, [isc:iec, jsc:jec, ksc:kec]
/// @param bxU       The u-point iteration box of this tile, [IscB:IecB, jsc:jec, ksc:kec]
/// @param bxV       The v-point iteration box of this tile, [isc:iec, JscB:JecB, ksc:kec]
/// @param u         Zonal velocity [L T-1 ~> m s-1]
/// @param v         Meridional velocity [L T-1 ~> m s-1]
/// @param uh        Zonal transport u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param vh        Meridional transport v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param q         Layer potential vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
/// @param uh_center Transport based on arithmetic mean h at u-points
///                  [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param vh_center Transport based on arithmetic mean h at v-points
///                  [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param dy_Cu     The unblocked lengths of the u-faces of the h-cell [L ~> m].
/// @param dx_Cv     The unblocked lengths of the v-faces of the h-cell [L ~> m].
/// @param IdxCu     1/dxCu [L-1 ~> m-1].
/// @param IdyCv     1/dyCv [L-1 ~> m-1].
/// @param CAu       Zonal acceleration due to Coriolis and momentum advection [L T-2 ~> m s-2].
/// @param CAv       Meridional acceleration due to Coriolis and momentum advection [L T-2 ~> m s-2].
/// @param CS        Options selecting the Coriolis discretization
///
void coradv_sadourny(const Box& bxH,
                     const Box& bxU,
                     const Box& bxV,
                     Array4<const Real> const& u,
                     Array4<const Real> const& v,
                     Array4<const Real> const& uh,
                     Array4<const Real> const& vh,
                     Array4<const Real> const& q,
                     Array4<const Real> const& uh_center,
                     Array4<const Real> const& vh_center,
                     Array4<const Real> const& dy_Cu,
                     Array4<const Real> const& dx_Cv,
                     Array4<const Real> const& IdxCu,
                     Array4<const Real> const& IdyCv,
                     Array4<Real>       const& CAu,
                     Array4<Real>       const& CAv,
                     const Coriolis_scheme_CS_C& CS);

///
/// @brief Calculates the Coriolis and momentum advection accelerations with the 3rd, 5th or
///        7th-order WENO PV reconstruction schemes over one iteration tile.
///
/// @param bxU                  The u-point iteration box of this tile, [IscB:IecB, jsc:jec, ksc:kec]
/// @param bxV                  The v-point iteration box of this tile, [isc:iec, JscB:JecB, ksc:kec]
/// @param u                    Zonal velocity [L T-1 ~> m s-1]
/// @param v                    Meridional velocity [L T-1 ~> m s-1]
/// @param uh                   Zonal transport u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param vh                   Meridional transport v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param q                    Layer potential vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
/// @param abs_vort             Absolute vorticity at q-points [T-1 ~> s-1]
/// @param h_q                  Thickness interpolated to q points [H ~> m or kg m-2]
/// @param IdxCu                1/dxCu [L-1 ~> m-1].
/// @param mask2dCu             0 for boundary points and 1 for ocean points on the u grid [nondim].
/// @param IdyCv                1/dyCv [L-1 ~> m-1].
/// @param mask2dCv             0 for boundary points and 1 for ocean points on the v grid [nondim].
/// @param H_subroundoff        A thickness that is so small that it can be added to a thickness of
///                             Angstrom or larger without changing it at the bit level
///                             [H ~> m or kg m-2].
/// @param CAu                  Zonal acceleration due to Coriolis and momentum advection [L T-2 ~> m s-2].
/// @param CAv                  Meridional acceleration due to Coriolis and momentum advection
///                             [L T-2 ~> m s-2].
/// @param weno_velocity_smooth If true, use velocity to compute the smoothness indicator for WENO
/// @param CS                   Options selecting the Coriolis discretization
///
void coradv_weno(const Box& bxU,
                 const Box& bxV,
                 Array4<const Real> const& u,
                 Array4<const Real> const& v,
                 Array4<const Real> const& uh,
                 Array4<const Real> const& vh,
                 Array4<const Real> const& q,
                 Array4<const Real> const& abs_vort,
                 Array4<const Real> const& h_q,
                 Array4<const Real> const& IdxCu,
                 Array4<const Real> const& mask2dCu,
                 Array4<const Real> const& IdyCv,
                 Array4<const Real> const& mask2dCv,
                 const Real H_subroundoff,
                 Array4<Real>       const& CAu,
                 Array4<Real>       const& CAv,
                 const bool weno_velocity_smooth,
                 const Coriolis_scheme_CS_C& CS);

///
/// @brief Calculates the Coriolis and momentum advection contributions to the acceleration over
///        the whole column, in iteration tiles of nkblock layers.
///
/// @param bxH0            The h-point box of the whole column, [isc:iec, jsc:jec, 1:nz]
/// @param bxQ0            The B-grid-index box of the whole column, [IscB:IecB, JscB:JecB, 1:nz]
/// @param u               Zonal velocity [L T-1 ~> m s-1]
/// @param v               Meridional velocity [L T-1 ~> m s-1]
/// @param h               Layer thickness [H ~> m or kg m-2]
/// @param uh              Zonal transport u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param vh              Meridional transport v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1]
/// @param CAu             Zonal acceleration due to Coriolis and momentum advection [L T-2 ~> m s-2].
/// @param CAv             Meridional acceleration due to Coriolis and momentum advection
///                        [L T-2 ~> m s-2].
/// @param RV              Diagnostic relative vorticity [T-1 ~> s-1]
/// @param PV              Diagnostic potential vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
/// @param CAuS            Stokes contribution to CAu [L T-2 ~> m s-2]
/// @param CAvS            Stokes contribution to CAv [L T-2 ~> m s-2]
/// @param OBC             Open boundary control structure
/// @param gradKEu         gradKEu = - d/dx(u2) [L T-2 ~> m s-2]; may be absent (.p == nullptr)
/// @param gradKEv         gradKEv = - d/dy(u2) [L T-2 ~> m s-2]; may be absent (.p == nullptr)
/// @param rv_x_v          rv_x_v = rv * v at u [L T-2 ~> m s-2]; may be absent (.p == nullptr)
/// @param rv_x_u          rv_x_u = rv * u at v [L T-2 ~> m s-2]; may be absent (.p == nullptr)
/// @param mask2dT         0 for land points and 1 for ocean points on the h-grid [nondim].
/// @param areaT           The area of an h-cell [L2 ~> m2].
/// @param areaCu          The areas of the u-grid cells [L2 ~> m2].
/// @param areaCv          The areas of the v-grid cells [L2 ~> m2].
/// @param IareaT          1/areaT [L-2 ~> m-2].
/// @param mask2dCu        0 for boundary points and 1 for ocean points on the u grid [nondim].
/// @param mask2dCv        0 for boundary points and 1 for ocean points on the v grid [nondim].
/// @param dyCv            dyCv is delta y at v points [L ~> m].
/// @param dxCu            dxCu is delta x at u points [L ~> m].
/// @param dy_Cu           The unblocked lengths of the u-faces of the h-cell [L ~> m].
/// @param dx_Cv           The unblocked lengths of the v-faces of the h-cell [L ~> m].
/// @param dyBu            dyBu is delta y at q points [L ~> m].
/// @param dxBu            dxBu is delta x at q points [L ~> m].
/// @param mask2dBu        0 for boundary points and 1 for ocean points on the q grid [nondim].
/// @param IareaBu         IareaBu = 1/areaBu [L-2 ~> m-2].
/// @param CoriolisBu      The Coriolis parameter at corner points [T-1 ~> s-1].
/// @param IdxCu           1/dxCu [L-1 ~> m-1].
/// @param IdyCv           1/dyCv [L-1 ~> m-1].
/// @param IdxCv           1/dxCv [L-1 ~> m-1].
/// @param IdyCu           1/dyCu [L-1 ~> m-1].
/// @param IdxCu_OBCmask   1/dxCu or 0 at boundary or OBC points [L-1 ~> m-1].
/// @param IdyCv_OBCmask   1/dxCv or 0 at boundary or OBC points [L-1 ~> m-1].
/// @param nz              The number of layers/levels in the vertical
/// @param H_subroundoff   A thickness that is so small that it can be added to a thickness of
///                        Angstrom or larger without changing it at the bit level
///                        [H ~> m or kg m-2].
/// @param Angstrom_H      A one-Angstrom thickness in the model thickness units [H ~> m or kg m-2].
/// @param m_to_L          A constant that translates lengths in meters to the units of
///                        horizontal lengths [L m-1 ~> 1]
/// @param m_s_to_L_T      Convert lateral velocities from m s-1 to L T-1 [L s T-1 m-1 ~> 1]
/// @param CS              Options used in the Coriolis and momentum advection calculations
/// @param por_face_areaU  fractional open area of U-faces [nondim]
/// @param por_face_areaV  fractional open area of V-faces [nondim]
/// @param Waves           An optional pointer to Stokes drift CS
/// @param Stokes_VF       If true, include the Stokes drift
///
void coradcalc_tr(const Box& bxH0,
                  const Box& bxQ0,
                  Array4<const Real> const& u,
                  Array4<const Real> const& v,
                  Array4<const Real> const& h,
                  Array4<const Real> const& uh,
                  Array4<const Real> const& vh,
                  Array4<Real>       const& CAu,
                  Array4<Real>       const& CAv,
                  Array4<Real>       const& RV,
                  Array4<Real>       const& PV,
                  Array4<Real>       const& CAuS,
                  Array4<Real>       const& CAvS,
                  OceanOBC* OBC,
                  Array4<Real>       const& gradKEu,
                  Array4<Real>       const& gradKEv,
                  Array4<Real>       const& rv_x_v,
                  Array4<Real>       const& rv_x_u,
                  Array4<const Real> const& mask2dT,
                  Array4<const Real> const& areaT,
                  Array4<const Real> const& areaCu,
                  Array4<const Real> const& areaCv,
                  Array4<const Real> const& IareaT,
                  Array4<const Real> const& mask2dCu,
                  Array4<const Real> const& mask2dCv,
                  Array4<const Real> const& dyCv,
                  Array4<const Real> const& dxCu,
                  Array4<const Real> const& dy_Cu,
                  Array4<const Real> const& dx_Cv,
                  Array4<const Real> const& dyBu,
                  Array4<const Real> const& dxBu,
                  Array4<const Real> const& mask2dBu,
                  Array4<const Real> const& IareaBu,
                  Array4<const Real> const& CoriolisBu,
                  Array4<const Real> const& IdxCu,
                  Array4<const Real> const& IdyCv,
                  Array4<const Real> const& IdxCv,
                  Array4<const Real> const& IdyCu,
                  Array4<const Real> const& IdxCu_OBCmask,
                  Array4<const Real> const& IdyCv_OBCmask,
                  const int nz,
                  const Real H_subroundoff,
                  const Real Angstrom_H,
                  const Real m_to_L,
                  const Real m_s_to_L_T,
                  const CorAdv_opts_CS_C& CS,
                  Array4<const Real> const& por_face_areaU,
                  Array4<const Real> const& por_face_areaV,
                  WaveParametersCS* Waves,
                  const bool Stokes_VF);

}
