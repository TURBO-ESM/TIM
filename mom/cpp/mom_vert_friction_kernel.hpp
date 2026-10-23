// SKILLS: 0.3.3
#pragma once
///
/// @file mom_vert_friction_kernel.hpp
/// @brief Per-point device primitives for the MOM6 vertical friction (viscosity) kernels.
///
#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Array4.H>
#include <AMReX_Math.H>

#include <cmath>

struct OceanOBC;  // Undefined at the moment

/// @brief Options used in the calculation of the coupling coefficients across interfaces.
/// Field-for-field mirror of the Fortran `bind(C)` type `coupling_CS_C` -- order and types must
/// not change.
struct coupling_CS_C {
    double Hmix;               ///< The mixed layer thickness [Z ~> m].
    double Kvml_invZ2;         ///< The extra vertical viscosity scale in [H Z T-1 ~> m2 s-1 or Pa s]
                               ///< in a surface mixed layer with a characteristic thickness given by
                               ///< Hmix, and scaling proportional to (Hmix/z)^2, where z is the
                               ///< distance from the surface; this can get very large with thin
                               ///< layers.
    double Kv;                 ///< The interior vertical viscosity [H Z T-1 ~> m2 s-1 or Pa s].
    double Kv_extra_bbl;       ///< An extra vertical viscosity in the bottom boundary layer of
                               ///< thickness Hbbl when there is not a bottom drag law in use
                               ///< [H Z T-1 ~> m2 s-1 or Pa s].
    double vonKar;             ///< The von Karman constant as used for mixed layer viscosity [nondim]
    bool   bottomdraglaw;      ///< If true, the  bottom stress is calculated with a
                               ///< drag law c_drag*|u|*u. The velocity magnitude
                               ///< may be an assumed value or it may be based on the
                               ///< actual velocity in the bottommost HBBL, depending
                               ///< on whether linear_drag is true.
    bool   dynamic_viscous_ML; ///< If true, use the results from a dynamic
                               ///< calculation, perhaps based on a bulk Richardson
                               ///< number criterion, to determine the mixed layer
                               ///< thickness for viscosity.
    bool   fixed_LOTW_ML;      ///< If true, use a Law-of-the-wall prescription for the mixed
                               ///< layer viscosity within a boundary layer that is the lesser of
                               ///< Hmix and the total depth of the ocean in a column.
    bool   apply_LOTW_floor;   ///< If true, use a Law-of-the-wall prescription to set a lower
                               ///< bound on the viscous coupling between layers within the
                               ///< surface boundary layer, based the distance of interfaces
                               ///< from the surface.  This only acts when there are large
                               ///< changes in the thicknesses of successive layers or when the
                               ///< viscosity is set externally and the wind stress has
                               ///< subsequently increased.
    int    answer_date;        ///< The vintage of the order of arithmetic and expressions in the
                               ///< viscous calculations.  Values below 20190101 recover the answers
                               ///< from the end of 2018, while higher values use expressions that do
                               ///< not use an arbitrary and hard-coded maximum viscous coupling
                               ///< coefficient between layers.  In non-Boussinesq cases, values
                               ///< below 20230601 recover a form of the viscosity within  the mixed
                               ///< layer that breaks up the magnitude of the wind stress with
                               ///< BULKMIXEDLAYER, DYNAMIC_VISCOUS_ML or FIXED_DEPTH_LOTW_ML, but
                               ///< not LOTW_VISCOUS_ML_FLOOR.
};

namespace MOM {
using amrex::Real;
using amrex::Array4;

///
/// @brief Computes the composite shape function and its derivative at sigma. This is a pure copy of
///        cvmix_kpp_composite_Gshape from pkg/CVMix-src, so that it can be called inside do concurrent.
///
/// @param sigma  Normalized boundary layer coordinate [nondim]
/// @param Gat1   The shape function at sigma = 1 [nondim]
/// @param Gsig   The shape function at sigma [nondim]
/// @param dGdsig The derivative of the shape function with sigma [nondim]
///
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void fpmix_gshape_point(Real const sigma, Real const Gat1, Real& Gsig, Real& dGdsig) noexcept
{
    /// Shape parameters [nondim]
    constexpr Real a2Gsig = -2.1637;
    constexpr Real a3Gsig =  0.5831;
    constexpr Real sig_m  =  0.35;
    constexpr Real G_m    =  0.11;

    /// Intermediate shape values [nondim]
    Real bGsig, G_1, sig;

    G_1    =  amrex::max( 0.0  ,  amrex::min( Gat1 , G_m ) );

    if (sigma < sig_m) {
        sig    = amrex::max( sigma , 0.0 );
        Gsig   = sig + sig * sig * (a2Gsig + a3Gsig * sig);
        dGdsig = 1.0 + sig * (2.0 * a2Gsig + 3.0 * a3Gsig * sig);
    } else {
        bGsig  =  (G_m-G_1) / ((1.0-sig_m) * (1.0-sig_m));
        sig    = amrex::min( sigma , 1.0 );
        Gsig   = G_1 + bGsig * (1.0 - sig) * (1.0 - sig);
        dGdsig = bGsig * 2.0 * (sig - 1.0);
    }
}

///
/// @brief Compute coupling coefficient associated with vertical viscosity parameterization as in Greatbatch and Lamb
///        (1990), hereafter referred to as the GL90 vertical viscosity parameterization. This vertical viscosity scheme
///        redistributes momentum in the vertical, and is the equivalent of the Gent & McWilliams (1990) parameterization,
///        but in a TWA (thickness-weighted averaged) set of equations. The vertical viscosity coefficient nu is computed
///        from kappa_GM via thermal wind balance, and the following relation:
///        nu = kappa_GM * f^2 / N^2.
///        In the following subroutine kappa_GM is assumed either (a) constant or (b) horizontally varying. In both cases,
///        (a) and  (b), one can additionally impose an EBT structure in the vertical for kappa_GM.
///        A third possible formulation of nu is depth-independent:
///        nu = f^2 * alpha
///        The latter formulation would be equivalent to a kappa_GM that varies as N^2 with depth.
///        The vertical viscosity del_z ( nu del_z u) is applied to the momentum equation with stress-free boundary
///        conditions at the top and bottom.
///
///        In SSW mode, we have 1/N^2 = h/g'. The coupling coefficient is therefore equal to
///        a_cpl_gl90 = nu / h = kappa_GM * f^2 / g'
///        or
///        a_cpl_gl90 = nu / h = f^2 * alpha / h
///
/// @tparam work_on_u      If true, u-points are being calculated, otherwise they are v-points.
/// @param i               Column i-index
/// @param j               Column j-index
/// @param a_cpl_gl90      Coupling coefficient associated with GL90 across interfaces; is not
///                        included in a_cpl [H T-1 ~> m s-1 or Pa s m-1].
/// @param hvel            Distance between interfaces at velocity points [Z ~> m]
/// @param z_i             Estimate of interface heights above the bottom, normalized by the GL90
///                        bottom boundary layer thickness [nondim]
/// @param ke              The number of layers/levels in the vertical
/// @param dZ_subroundoff  A thickness in height units that is so small that it can be added to a
///                        vertical distance of Angstrom_Z or 1e-17 m without changing it at the
///                        bit level [Z ~> m].  This is the height equivalent of H_subroundoff.
/// @param g_prime         The reduced gravity at each interface [L2 Z-1 T-2 ~> m s-2].
/// @param CoriolisBu      The Coriolis parameter at corner points [T-1 ~> s-1].
/// @param use_GL90_N2     If true, use GL90 vertical viscosity coefficient that is
///                        depth-independent; this corresponds to a kappa_GM that scales as N^2
///                        with depth.
/// @param kappa_gl90      The scalar diffusivity used in the GL90 vertical viscosity scheme
///                        [L2 H Z-1 T-1 ~> m2 s-1 or Pa s]
/// @param read_kappa_gl90 If true, read a file containing the spatially varying kappa_gl90
/// @param alpha_gl90      Coefficient used to compute a depth-independent GL90 vertical
///                        viscosity via Kv_gl90 = alpha_gl90 * f^2. Note that the implied
///                        Kv_gl90 corresponds to a kappa_gl90 that scales as N^2 with depth.
///                        [H Z T ~> m2 s or kg s m-1]
/// @param kappa_gl90_2d   2D kappa_gl90 at h-points [L2 H Z-1 T-1 ~> m2 s-1 or Pa s]
/// @param kdgl90_struct   Vertical structure function used in GL90 diffusivity [nondim];
///                        may be unallocated; pass a default-constructed Array4 when unallocated
///
template <bool work_on_u>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void find_coupling_coef_gl90_point(int i, int j,
                                   Array4<Real>       const& a_cpl_gl90,
                                   Array4<const Real> const& hvel,
                                   Array4<const Real> const& z_i,
                                   int const ke,
                                   Real const dZ_subroundoff,
                                   Array4<const Real> const& g_prime,
                                   Array4<const Real> const& CoriolisBu,
                                   bool const use_GL90_N2,
                                   Real const kappa_gl90,
                                   bool const read_kappa_gl90,
                                   Real const alpha_gl90,
                                   Array4<const Real> const& kappa_gl90_2d,
                                   Array4<const Real> const& kdgl90_struct) noexcept
{
    bool kdgl90_use_vert_struct;  ///< use vertical structure for GL90 coefficient
    Real f2;                      ///< Squared Coriolis parameter at a velocity grid point [T-2 ~> s-2].
    Real h_neglect;               ///< A vertical distance that is so small it is usually lost in
                                  ///< roundoff error and can be neglected [Z ~> m].
    Real botfn;                   ///< A function that is 1 at the bottom and small far from it [nondim]
    Real z2;                      ///< The distance from the bottom, normalized by Hbbl_gl90 [nondim]

    const int nz = ke;
    h_neglect = dZ_subroundoff;
    kdgl90_use_vert_struct = false;

    kdgl90_use_vert_struct = (kdgl90_struct.p != nullptr);

    for (int k = 0; k <= nz; ++k) a_cpl_gl90(i,j,k) = 0.0;

    for (int k = 1; k < nz; ++k) {
        if constexpr (work_on_u) {
            /// compute coupling coefficient at u-points
            Real const fsum = (CoriolisBu(i,j-1,0) + CoriolisBu(i,j,0));
            f2 = 0.25 * (fsum * fsum);
            if (use_GL90_N2) {
                a_cpl_gl90(i,j,k) = 2.0 * f2 * alpha_gl90 / (hvel(i,j,k) + hvel(i,j,k-1) + h_neglect);
            } else {
                if (read_kappa_gl90) {
                    a_cpl_gl90(i,j,k) = f2 * 0.5 * (kappa_gl90_2d(i,j,0) + kappa_gl90_2d(i+1,j,0)) / g_prime(k,0,0);
                } else {
                    a_cpl_gl90(i,j,k) = f2 * kappa_gl90 / g_prime(k,0,0);
                }
                if (kdgl90_use_vert_struct) {
                    a_cpl_gl90(i,j,k) = a_cpl_gl90(i,j,k) * 0.5
                        * (kdgl90_struct(i,j,k-1) + kdgl90_struct(i+1,j,k-1));
                }
            }
            /// botfn determines when a point is within the influence of the GL90 bottom boundary
            /// layer, going from 1 at the bottom to 0 in the interior.
            z2 = z_i(i,j,k);
            botfn = 1.0 / (1.0 + 0.09 * z2 * z2 * z2 * z2 * z2 * z2);

            a_cpl_gl90(i,j,k) = a_cpl_gl90(i,j,k) * (1.0 - botfn);
        } else {
            /// compute viscosities at v-points
            Real const fsum = (CoriolisBu(i-1,j,0) + CoriolisBu(i,j,0));
            f2 = 0.25 * (fsum * fsum);

            if (use_GL90_N2) {
                a_cpl_gl90(i,j,k) = 2.0 * f2 * alpha_gl90 / (hvel(i,j,k) + hvel(i,j,k-1) + h_neglect);
            } else {
                if (read_kappa_gl90) {
                    a_cpl_gl90(i,j,k) = f2 * 0.5 * (kappa_gl90_2d(i,j,0) + kappa_gl90_2d(i,j+1,0)) / g_prime(k,0,0);
                } else {
                    a_cpl_gl90(i,j,k) = f2 * kappa_gl90 / g_prime(k,0,0);
                }
                if (kdgl90_use_vert_struct) {
                    a_cpl_gl90(i,j,k) = a_cpl_gl90(i,j,k) * 0.5
                        * (kdgl90_struct(i,j,k-1) + kdgl90_struct(i,j+1,k-1));
                }
            }
            /// botfn determines when a point is within the influence of the GL90 bottom boundary
            /// layer, going from 1 at the bottom to 0 in the interior.
            z2 = z_i(i,j,k);
            botfn = 1.0 / (1.0 + 0.09 * z2 * z2 * z2 * z2 * z2 * z2);

            a_cpl_gl90(i,j,k) = a_cpl_gl90(i,j,k) * (1.0 - botfn);
        }
    }
}

///
/// @brief Calculate the 'coupling coefficient' (a_cpl) at the interfaces.
///        If BOTTOMDRAGLAW is defined, the minimum of Hbbl and half the adjacent
///        layer thicknesses are used to calculate a_cpl near the bottom.
///
/// @tparam work_on_u       If true, u-points are being calculated, otherwise they are v-points
/// @tparam do_shelf        If true, use a surface boundary condition appropriate for an ice shelf.
/// @param i                Column i-index
/// @param j                Column j-index
/// @param a_cpl            Coupling coefficient across interfaces [H T-1 ~> m s-1 or Pa s m-1]
/// @param hvel             Distance between interfaces at velocity points [Z ~> m]
/// @param h_harm           Harmonic mean of thicknesses around a velocity grid point [Z ~> m]
/// @param bbl_thick        Bottom boundary layer thickness [Z ~> m]
/// @param kv_bbl           Bottom boundary layer viscosity, exclusive of any depth-dependent
///                         contributions from visc%Kv_shear [H Z T-1 ~> m2 s-1 or Pa s]
/// @param z_i              Estimate of interface heights above the bottom, normalized by the
///                         bottom boundary layer thickness [nondim]
/// @param h_ml             Mixed layer depth [Z ~> m]
/// @param dt               Time increment [T ~> s]
/// @param ke               The number of layers/levels in the vertical
/// @param Boussinesq       If true, make the Boussinesq approximation.
/// @param nkml             The number of layers at the top that should be treated as parts of a
///                         homogeneous region.
/// @param dZ_subroundoff   A thickness in height units that is so small that it can be added to a
///                         vertical distance of Angstrom_Z or 1e-17 m without changing it at the
///                         bit level [Z ~> m].  This is the height equivalent of H_subroundoff.
/// @param H_to_m           A constant that translates distances in the units of thickness to m
///                         [m H-1 ~> 1 or m3 kg-1].
/// @param RZ_to_H          A constant that translates mass per unit area units to thickness
///                         units [H R-1 Z-1 ~> m3 kg-2 or 1].
/// @param Z_to_H           A constant that translates depth units to thickness units depth
///                         [H Z-1 ~> 1 or kg m-3].
/// @param CoriolisBu       The Coriolis parameter at corner points [T-1 ~> s-1].
/// @param CS               Vertical viscosity control structure
/// @param Kv_shear         The shear-driven turbulent vertical viscosity at the interfaces between
///                         layers in tracer columns [H Z T-1 ~> m2 s-1 or Pa s];
///                         may be unallocated; pass a default-constructed Array4 when unallocated
/// @param Kv_shear_Bu      The shear-driven turbulent vertical viscosity at the interfaces between
///                         layers in corner columns [H Z T-1 ~> m2 s-1 or Pa s];
///                         may be unallocated; pass a default-constructed Array4 when unallocated
/// @param Kv_tbl_shelf_u   Viscosity in the viscous top boundary layer under ice shelves at
///                         u-points [H Z T-1 ~> m2 s-1 or Pa s]
/// @param Kv_tbl_shelf_v   Viscosity in the viscous top boundary layer under ice shelves at
///                         v-points [H Z T-1 ~> m2 s-1 or Pa s]
/// @param tbl_thick_shelf_u Thickness of the viscous top boundary layer under ice shelves at
///                         u-points [Z ~> m].
/// @param tbl_thick_shelf_v Thickness of the viscous top boundary layer under ice shelves at
///                         v-points [Z ~> m].
/// @param nkml_visc_u      The number of layers in the viscous surface mixed layer at u-points
///                         [nondim]. This is not an integer because there may be fractional
///                         layers, and it is stored in terms of layers, not depth, to facilitate
///                         the movement of the viscous boundary layer with the flow.
/// @param nkml_visc_v      The number of layers in the viscous surface mixed layer at v-points
///                         [nondim].
/// @param Ustar_2d         The wind friction velocity, calculated using the Boussinesq reference
///                         density or the time-evolving surface density in non-Boussinesq mode
///                         [Z T-1 ~> m s-1]
/// @param SpV_avg          The layer averaged in situ specific volume [R-1 ~> m3 kg-1];
///                         may be unallocated; pass a default-constructed Array4 when unallocated
///
template <bool work_on_u, bool do_shelf>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void find_coupling_coef_point(int i, int j,
                              Array4<Real>       const& a_cpl,
                              Array4<const Real> const& hvel,
                              Array4<const Real> const& h_harm,
                              Real const bbl_thick,
                              Real const kv_bbl,
                              Array4<const Real> const& z_i,
                              Real& h_ml,
                              Real const dt,
                              int const ke,
                              bool const Boussinesq,
                              int const nkml,
                              Real const dZ_subroundoff,
                              Real const H_to_m,
                              Real const RZ_to_H,
                              Real const Z_to_H,
                              Array4<const Real> const& CoriolisBu,
                              coupling_CS_C const& CS,
                              Array4<const Real> const& Kv_shear,
                              Array4<const Real> const& Kv_shear_Bu,
                              Array4<const Real> const& Kv_tbl_shelf_u,
                              Array4<const Real> const& Kv_tbl_shelf_v,
                              Array4<const Real> const& tbl_thick_shelf_u,
                              Array4<const Real> const& tbl_thick_shelf_v,
                              Array4<const Real> const& nkml_visc_u,
                              Array4<const Real> const& nkml_visc_v,
                              Array4<const Real> const& Ustar_2d,
                              Array4<const Real> const& SpV_avg) noexcept
{
    Real u_star = 0.0;   ///< ustar at a velocity point [Z T-1 ~> m s-1]
    Real tau_mag = 0.0;  ///< The magnitude of the wind stress at a velocity point including
                         ///< gustiness [H Z T-2 ~> m2 s-2 or Pa]
    Real absf = 0.0;     ///< The average of the neighboring absolute values of f [T-1 ~> s-1].
    Real rho_av1;        ///< The harmonic mean surface layer density at velocity points [R ~> kg m-3]
    Real z_t = 0.0;      ///< The distance from the top, sometimes normalized
                         ///< by Hmix, [Z ~> m] or [nondim].
    Real kv_TBL;         ///< The viscosity in a top boundary layer under ice [H Z T-1 ~> m2 s-1 or Pa s]
    Real tbl_thick;      ///< The thickness of the top boundary layer [Z ~> m]
    Real Kv_add;         ///< A viscosity to add [H Z T-1 ~> m2 s-1 or Pa s]
    Real Kv_tot;         ///< The total viscosity at an interface [H Z T-1 ~> m2 s-1 or Pa s]
    int  nk_in_ml;       ///< The index of the deepest interface in the mixed layer.
    Real h_shear;        ///< The distance over which shears occur [Z ~> m].
    Real dhc;            ///< The distance between the center of adjacent layers [Z ~> m].
    Real visc_ml;        ///< The mixed layer viscosity [H Z T-1 ~> m2 s-1 or Pa s].
    Real I_Hmix = 0.0;   ///< The inverse of the mixed layer thickness [Z-1 ~> m-1].
    Real a_ml;           ///< The layer coupling coefficient across an interface in
                         ///< the mixed layer [H T-1 ~> m s-1 or Pa s m-1].
    Real a_floor;        ///< A lower bound on the layer coupling coefficient across an interface in
                         ///< the mixed layer [H T-1 ~> m s-1 or Pa s m-1].
    Real I_amax;         ///< The inverse of the maximum coupling coefficient [T H-1 ~> s m-1 or s m2 kg-1].
    Real temp1;          ///< A temporary variable [Z2 ~> m2]
    Real ustar2_denom;   ///< A temporary variable in the surface boundary layer turbulence
                         ///< calculations [H Z-1 T-1 ~> s-1 or kg m-3 s-1]
    Real h_neglect;      ///< A vertical distance that is so small it is usually lost
                         ///< in roundoff and can be neglected [Z ~> m].
    Real z2;             ///< A copy of z_i [nondim]
    Real botfn;          ///< A function that is 1 at the bottom and small far from it [nondim]
    Real topfn;          ///< A function that is 1 at the top and small far from it [nondim]
    Real kv_top;         ///< A viscosity associated with the top boundary layer [H Z T-1 ~> m2 s-1 or Pa s]
    bool do_OBCs, can_exit;
    int  max_nk = 0;

    const int nz = ke;

    h_neglect = dZ_subroundoff;

    if (CS.answer_date < 20190101) {
        ///   The maximum coupling coefficient was originally introduced to avoid
        /// truncation error problems in the tridiagonal solver. Effectively, the 1e-10
        /// sets the maximum coupling coefficient increment to 1e10 m per timestep.
        I_amax = (1.0e-10*H_to_m) * dt;
    } else {
        I_amax = 0.0;
    }

    do_OBCs = false;
    /*
    if (OBC != nullptr) {
        if constexpr (work_on_u) {
            do_OBCs = OBC->u_E_OBCs_on_PE || OBC->u_W_OBCs_on_PE;
        } else {
            do_OBCs = OBC->v_N_OBCs_on_PE || OBC->v_S_OBCs_on_PE;
        }
    }
    */

    for (int k = 0; k <= nz; ++k) a_cpl(i,j,k) = 0.0;
    h_ml = 0.0;

    if (CS.Kvml_invZ2 > 0.0 && !do_shelf) {
        I_Hmix = 1.0 / (CS.Hmix + h_neglect);
        z_t = h_neglect * I_Hmix;
    }

    for (int k = 1; k < nz; ++k) {
        Kv_tot = CS.Kv;

        if (CS.Kvml_invZ2 > 0.0 && !do_shelf) {
            /// This is an older (vintage ~1997) way to prevent wind stresses from driving very
            /// large flows in nearly massless near-surface layers when there is not a physically-
            /// based surface boundary layer parameterization.  It does not have a plausible
            /// physical basis, and probably should not be used.
            z_t = z_t + h_harm(i,j,k-1) * I_Hmix;
            Kv_tot = CS.Kv + CS.Kvml_invZ2 / ((z_t * z_t) *
                     (1.0 + 0.09 * z_t * z_t * z_t * z_t * z_t * z_t));
        }

        if (Kv_shear.p != nullptr) {
            /// Add in viscosities that are determined by physical processes that are handled in
            /// other modules, and which do not respond immediately to the changing layer
            /// thicknesses.  These processes may include shear-driven mixing or contributions
            /// from some boundary layer turbulence schemes.  Other viscosity contributions that
            /// respond to the evolving layer thicknesses or the surface wind stresses are added
            /// later.
            if constexpr (work_on_u) {
                Kv_add = 0.5 * (Kv_shear(i,j,k) + Kv_shear(i+1,j,k));
                /*
                if (do_OBCs) {
                    if (OBC->u_E_OBCs_on_PE) {
                        if (OBC->segnum_u(i,j,0) > 0) Kv_add = Kv_shear(i,j,k);
                    }
                    if (OBC->u_W_OBCs_on_PE) {
                        if (OBC->segnum_u(i,j,0) < 0) Kv_add = Kv_shear(i+1,j,k);
                    }
                }
                */
                Kv_tot = Kv_tot + Kv_add;
            } else {
                Kv_add = 0.5 * (Kv_shear(i,j,k) + Kv_shear(i,j+1,k));
                /*
                if (do_OBCs) {
                    if (OBC->v_N_OBCs_on_PE) {
                        if (OBC->segnum_v(i,j,0) > 0) Kv_add = Kv_shear(i,j,k);
                    }
                    if (OBC->v_S_OBCs_on_PE) {
                        if (OBC->segnum_v(i,j,0) < 0) Kv_add = Kv_shear(i,j+1,k);
                    }
                }
                */
                Kv_tot = Kv_tot + Kv_add;
            }
        }

        if (Kv_shear_Bu.p != nullptr) {
            /// This is similar to what was done above, but for contributions coming from the
            /// corner (vorticity) points.  Because OBCs run through the faces and corners there
            /// is no need to further modify these viscosities here to take OBCs into account.
            if constexpr (work_on_u) {
                Kv_tot = Kv_tot + 0.5 * (Kv_shear_Bu(i,j-1,k) + Kv_shear_Bu(i,j,k));
            } else {
                Kv_tot = Kv_tot + 0.5 * (Kv_shear_Bu(i-1,j,k) + Kv_shear_Bu(i,j,k));
            }
        }

        /// Set the viscous coupling coefficients, excluding surface mixed layer contributions
        /// for now, but including viscous bottom drag, working up from the bottom.
        if (CS.bottomdraglaw) {
            ///    botfn determines when a point is within the influence of the bottom
            ///  boundary layer, going from 1 at the bottom to 0 in the interior.
            z2 = z_i(i,j,k);
            botfn = 1.0 / (1.0 + 0.09 * z2 * z2 * z2 * z2 * z2 * z2);

            Kv_tot = Kv_tot + (kv_bbl - CS.Kv) * botfn;
            dhc = 0.5 * (hvel(i,j,k) + hvel(i,j,k-1));
            if (dhc > bbl_thick) {
                h_shear = ((1.0 - botfn) * dhc + botfn * bbl_thick) + h_neglect;
            } else {
                h_shear = dhc + h_neglect;
            }

            /// Calculate the coupling coefficients from the viscosities.
            a_cpl(i,j,k) = Kv_tot / (h_shear + (I_amax * Kv_tot));
        } else if (amrex::Math::abs(CS.Kv_extra_bbl) > 0.0) {
            /// There is a simple enhancement of the near-bottom viscosities, but no
            /// adjustment of the viscous coupling length scales to give a particular
            /// bottom stress.

            ///    botfn determines when a point is within the influence of the bottom
            ///  boundary layer, going from 1 at the bottom to 0 in the interior.
            z2 = z_i(i,j,k);
            botfn = 1.0 / (1.0 + 0.09 * z2 * z2 * z2 * z2 * z2 * z2);

            Kv_tot = Kv_tot + CS.Kv_extra_bbl * botfn;
            h_shear = 0.5 * (hvel(i,j,k) + hvel(i,j,k-1) + h_neglect);

            /// Calculate the coupling coefficients from the viscosities.
            a_cpl(i,j,k) = Kv_tot / (h_shear + I_amax * Kv_tot);
        } else {
            /// Any near-bottom viscous enhancements were already incorporated into
            /// Kv_tot, and there is no adjustment of the viscous coupling length
            /// scales to give a particular bottom stress.

            h_shear = 0.5 * (hvel(i,j,k) + hvel(i,j,k-1) + h_neglect);
            /// Calculate the coupling coefficients from the viscosities.
            a_cpl(i,j,k) = Kv_tot / (h_shear + I_amax * Kv_tot);
        }
    }

    /// Assign the bottom coupling coefficients
    if (CS.bottomdraglaw) {
        dhc = hvel(i,j,nz-1) * 0.5;
        a_cpl(i,j,nz) = kv_bbl / ((amrex::min(dhc, bbl_thick) + h_neglect) + I_amax * kv_bbl);
    } else if (amrex::Math::abs(CS.Kv_extra_bbl) > 0.0) {
        a_cpl(i,j,nz) = (CS.Kv + CS.Kv_extra_bbl)
            / ((0.5 * hvel(i,j,nz-1) + h_neglect) + I_amax * (CS.Kv + CS.Kv_extra_bbl));
    } else {
        a_cpl(i,j,nz) = CS.Kv / ((0.5 * hvel(i,j,nz-1) + h_neglect) + I_amax * CS.Kv);
    }

    /// Add surface intensified viscous coupling, either as a no-slip boundary condition under a
    /// rigid ice-shelf, or due to wind-stress driven surface boundary layer mixing that has not
    /// already been added via visc%Kv_shear.
    if constexpr (do_shelf) {
        /// Set the coefficients to include the no-slip surface stress.
        if constexpr (work_on_u) {
            kv_TBL = Kv_tbl_shelf_u(i,j,0);
            tbl_thick = tbl_thick_shelf_u(i,j,0) + h_neglect;
        } else {
            kv_TBL = Kv_tbl_shelf_v(i,j,0);
            tbl_thick = tbl_thick_shelf_v(i,j,0) + h_neglect;
        }
        z_t = 0.0;

        /// If a_cpl(1) were not already 0, it would be added here.
        if (0.5 * hvel(i,j,0) > tbl_thick) {
            a_cpl(i,j,0) = kv_TBL / (tbl_thick + I_amax * kv_TBL);
        } else {
            a_cpl(i,j,0) = kv_TBL / ((0.5 * hvel(i,j,0) + h_neglect) + I_amax * kv_TBL);
        }

        for (int k = 1; k < nz; ++k) {
            z_t = z_t + hvel(i,j,k-1) / tbl_thick;
            Real const z_t3 = (z_t * z_t) * z_t;
            topfn = 1.0 / (1.0 + 0.09 * (z_t3 * z_t3));

            dhc = 0.5 * (hvel(i,j,k) + hvel(i,j,k-1));
            if (dhc > tbl_thick) {
                h_shear = ((1.0 - topfn) * dhc + topfn * tbl_thick) + h_neglect;
            } else {
                h_shear = dhc + h_neglect;
            }

            kv_top = topfn * kv_TBL;
            a_cpl(i,j,k) = a_cpl(i,j,k) + kv_top / (h_shear + I_amax * kv_top);
        }
    } else if (CS.dynamic_viscous_ML || (nkml > 0) || CS.fixed_LOTW_ML || CS.apply_LOTW_floor) {

        /// Find the friction velocity and the absolute value of the Coriolis parameter at this
        /// point.  Zero out the friction velocity on land points.
        u_star = 0.0;
        tau_mag = 0.0;

        if (SpV_avg.p != nullptr) {
            rho_av1 = 0.0;

            if constexpr (work_on_u) {
                u_star = 0.5 * (Ustar_2d(i,j,0) + Ustar_2d(i+1,j,0));
                rho_av1 = 2.0 / (SpV_avg(i,j,0) + SpV_avg(i+1,j,0));
                absf = 0.5 * (amrex::Math::abs(CoriolisBu(i,j-1,0)) + amrex::Math::abs(CoriolisBu(i,j,0)));
                /*
                if (do_OBCs) {
                    if (OBC->u_E_OBCs_on_PE) {
                        if (OBC->segnum_u(i,j,0) > 0) {
                            u_star = Ustar_2d(i,j,0);
                            rho_av1 = 1.0 / SpV_avg(i,j,0);
                        }
                    }
                    if (OBC->u_W_OBCs_on_PE) {
                        if (OBC->segnum_u(i,j,0) < 0) {
                            u_star = Ustar_2d(i+1,j,0);
                            rho_av1 = 1.0 / SpV_avg(i+1,j,0);
                        }
                    }
                }
                */
            } else {
                u_star = 0.5 * (Ustar_2d(i,j,0) + Ustar_2d(i,j+1,0));
                rho_av1 = 2.0 / (SpV_avg(i,j,0) + SpV_avg(i,j+1,0));
                absf = 0.5 * (amrex::Math::abs(CoriolisBu(i-1,j,0)) + amrex::Math::abs(CoriolisBu(i,j,0)));
                /*
                if (do_OBCs) {
                    if (OBC->v_N_OBCs_on_PE) {
                        if (OBC->segnum_v(i,j,0) > 0) {
                            u_star = Ustar_2d(i,j,0);
                            rho_av1 = 1.0 / SpV_avg(i,j,0);
                        }
                    }
                    if (OBC->v_S_OBCs_on_PE) {
                        if (OBC->segnum_v(i,j,0) < 0) {
                            u_star = Ustar_2d(i,j+1,0);
                            rho_av1 = 1.0 / SpV_avg(i,j+1,0);
                        }
                    }
                }
                */
            }

            tau_mag = RZ_to_H * rho_av1 * (u_star * u_star);
        } else {  /// (.not.allocated(tv%SpV_avg))
            if constexpr (work_on_u) {
                u_star = 0.5 * (Ustar_2d(i,j,0) + Ustar_2d(i+1,j,0));
                absf = 0.5 * (amrex::Math::abs(CoriolisBu(i,j-1,0)) + amrex::Math::abs(CoriolisBu(i,j,0)));
                /*
                if (do_OBCs) {
                    if (OBC->u_E_OBCs_on_PE) {
                        if (OBC->segnum_u(i,j,0) > 0) u_star = Ustar_2d(i,j,0);
                    }
                    if (OBC->u_W_OBCs_on_PE) {
                        if (OBC->segnum_u(i,j,0) < 0) u_star = Ustar_2d(i+1,j,0);
                    }
                }
                */
            } else {
                u_star = 0.5 * (Ustar_2d(i,j,0) + Ustar_2d(i,j+1,0));
                absf = 0.5 * (amrex::Math::abs(CoriolisBu(i-1,j,0)) + amrex::Math::abs(CoriolisBu(i,j,0)));
                /*
                if (do_OBCs) {
                    if (OBC->v_N_OBCs_on_PE) {
                        if (OBC->segnum_v(i,j,0) > 0) u_star = Ustar_2d(i,j,0);
                    }
                    if (OBC->v_S_OBCs_on_PE) {
                        if (OBC->segnum_v(i,j,0) < 0) u_star = Ustar_2d(i,j+1,0);
                    }
                }
                */
            }

            tau_mag = Z_to_H * (u_star * u_star);
        }

        /// Determine the thickness of the surface ocean boundary layer and its extent in index
        /// space.
        nk_in_ml = 0;

        if (CS.dynamic_viscous_ML) {
            /// The fractional number of layers that are within the viscous boundary layer were
            /// previously stored in visc%nkml_visc_[uv].
            h_ml = h_neglect;
            max_nk = 0;

            if constexpr (work_on_u) {
                nk_in_ml = static_cast<int>(std::ceil(nkml_visc_u(i,j,0)));
                max_nk = amrex::max(max_nk, nk_in_ml);

                for (int k = 0; k < max_nk; ++k) {
                    if ((k+1) <= nkml_visc_u(i,j,0)) {  /// This layer is all in the ML.
                        h_ml = h_ml + hvel(i,j,k);
                    } else if ((k+1) < nkml_visc_u(i,j,0) + 1.0) {  /// Part of this layer is in the ML.
                        h_ml = h_ml + ((nkml_visc_u(i,j,0) + 1.0) - (k+1)) * hvel(i,j,k);
                    }
                }
            } else {
                nk_in_ml = static_cast<int>(std::ceil(nkml_visc_v(i,j,0)));
                max_nk = amrex::max(max_nk, nk_in_ml);

                for (int k = 0; k < max_nk; ++k) {
                    if ((k+1) <= nkml_visc_v(i,j,0)) {  /// This layer is all in the ML.
                        h_ml = h_ml + hvel(i,j,k);
                    } else if ((k+1) < nkml_visc_v(i,j,0) + 1.0) {  /// Part of this layer is in the ML.
                        h_ml = h_ml + ((nkml_visc_v(i,j,0) + 1.0) - (k+1)) * hvel(i,j,k);
                    }
                }
            }
        } else if (nkml > 0) {
            /// This is a simple application of a refined-bulk mixed layer with GV%nkml sublayers.
            max_nk = nkml;
            nk_in_ml = nkml;

            h_ml = h_neglect;

            for (int k = 0; k < nkml; ++k) {
                h_ml = h_ml + hvel(i,j,k);
            }
        } else if (CS.fixed_LOTW_ML || CS.apply_LOTW_floor) {
            /// Determine which interfaces are within CS%Hmix of the surface, and set the viscous
            /// boundary layer thickness to the the smaller of CS%Hmix and the depth of the ocean.
            h_ml = 0.0;
            for (int k = 0; k < nz; ++k) {
                can_exit = true;
                if (h_ml < CS.Hmix) {
                    nk_in_ml = k+1;

                    if (h_ml + hvel(i,j,k) < CS.Hmix) {
                        h_ml = h_ml + hvel(i,j,k);
                        can_exit = false;  /// Part of the next deeper layer is also in the mixed layer.
                    } else {
                        h_ml = CS.Hmix;
                    }
                }

                if (can_exit) break;  /// All remaining layers in this row are below the mixed layer depth.
            }

            /// This made more sense in the slab/layer solvers...
            max_nk = 0;
            max_nk = amrex::max(max_nk, nk_in_ml);
        }

        /// Avoid working on land or on columns where the viscous coupling could not be increased.
        if (u_star <= 0.0) nk_in_ml = 0;

        /// Set the viscous coupling at the interfaces as the larger of what was previously
        /// set and the contributions from the surface boundary layer.
        z_t = 0.0;
        if (CS.apply_LOTW_floor &&
            (CS.dynamic_viscous_ML || nkml > 0 || CS.fixed_LOTW_ML)) {
            for (int k = 1; k < max_nk; ++k) {
                if ((k+1) <= nk_in_ml) {
                    z_t = z_t + hvel(i,j,k-1);

                    ///   The viscosity in visc_ml is set to go to 0 at the mixed layer top and
                    /// bottom (in a log-layer) and be further limited by rotation to give the
                    /// natural Ekman length.
                    temp1 = (z_t * h_ml - z_t * z_t);
                    if (Boussinesq) {
                        ustar2_denom = (CS.vonKar * Z_to_H * (u_star * u_star))
                            / (absf * temp1 + (h_ml + h_neglect) * u_star);
                    } else {
                        ustar2_denom = (CS.vonKar * tau_mag)
                            / (absf * temp1 + (h_ml + h_neglect) * u_star);
                    }

                    visc_ml = temp1 * ustar2_denom;
                    /// Set the viscous coupling based on the model's vertical resolution.  The
                    /// omission of the I_amax factor here is consistent with answer dates above
                    /// 20190101.
                    a_ml = visc_ml / (0.25 * (hvel(i,j,k) + hvel(i,j,k-1) + h_neglect));

                    /// As a floor on the viscous coupling, assume that the length scale in the
                    /// denominator can not be larger than the distance from the surface,
                    /// consistent with a logarithmic velocity profile.  This is consistent with
                    /// visc_ml, but cancels out common factors of z_t.
                    a_floor = (h_ml - z_t) * ustar2_denom;

                    /// Choose the largest estimate of a_cpl.
                    a_cpl(i,j,k) = amrex::max(amrex::max(a_cpl(i,j,k), a_ml), a_floor);
                    /// An option could be added to change this to: a_cpl(i,K) = max(a_cpl(i,K) + a_ml, a_floor)
                }
            }
        } else if (CS.apply_LOTW_floor) {
            for (int k = 1; k < max_nk; ++k) {
                if ((k+1) <= nk_in_ml) {
                    z_t = z_t + hvel(i,j,k-1);

                    temp1 = (z_t * h_ml - z_t * z_t);
                    if (Boussinesq) {
                        ustar2_denom = (CS.vonKar * Z_to_H * (u_star * u_star))
                            / (absf * temp1 + (h_ml + h_neglect) * u_star);
                    } else {
                        ustar2_denom = (CS.vonKar * tau_mag)
                            / (absf * temp1 + (h_ml + h_neglect) * u_star);
                    }

                    /// As a floor on the viscous coupling, assume that the length scale in the
                    /// denominator can not be larger than the distance from the surface,
                    /// consistent with a logarithmic velocity profile.
                    a_cpl(i,j,k) = amrex::max(a_cpl(i,j,k), (h_ml - z_t) * ustar2_denom);
                }
            }
        } else {
            for (int k = 1; k < max_nk; ++k) {
                if ((k+1) <= nk_in_ml) {
                    z_t = z_t + hvel(i,j,k-1);

                    temp1 = (z_t * h_ml - z_t * z_t);
                    ///   This viscosity is set to go to 0 at the mixed layer top and bottom (in a
                    /// log-layer) and be further limited by rotation to give the natural Ekman
                    /// length.  The following expressions are mathematically equivalent.
                    if (Boussinesq || (CS.answer_date < 20230601)) {
                        visc_ml = u_star * CS.vonKar * (Z_to_H * temp1 * u_star)
                            / (absf * temp1 + (h_ml + h_neglect) * u_star);
                    } else {
                        visc_ml = CS.vonKar * (temp1 * tau_mag)
                            / (absf * temp1 + (h_ml + h_neglect) * u_star);
                    }
                    a_ml = visc_ml / (0.25 * (hvel(i,j,k) + hvel(i,j,k-1) + h_neglect) + 0.5 * I_amax * visc_ml);

                    /// Choose the largest estimate of a_cpl, but these could be changed to be
                    /// additive.
                    a_cpl(i,j,k) = amrex::max(a_cpl(i,j,k), a_ml);
                    /// An option could be added to change this to: a_cpl(i,K) = a_cpl(i,K) + a_ml
                }
            }
        }
    }
    (void)do_OBCs;
}

} // namespace MOM
