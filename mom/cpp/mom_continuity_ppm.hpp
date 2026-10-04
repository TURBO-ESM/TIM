// mom_continuity_ppm.hpp
// SKILLS: 0.3.1
#pragma once
/// @file mom_continuity_ppm.hpp
/// @brief Box-level AMReX kernel declarations for MOM6 PPM continuity
///        (piecewise parabolic reconstruction and edge-thickness routines).

#include "mom_continuity_ppm_kernel.hpp"

struct OceanOBC;    // Undefined at the moment

/// @brief Options controlling the transport adjustment and barotropic-consistency
/// iteration used by the continuity solver. Field-for-field mirror of the Fortran
/// `bind(C)` type `transport_adjust_CS_C` -- order and types must not change.
struct transport_adjust_CS_C {
    double tol_eta;            ///< Tolerance for free-surface height discrepancies.
    double tol_vel;            ///< Tolerance for barotropic velocity discrepancies.
    double CFL_limit_adjust;   ///< Maximum CFL of the adjusted velocities.
    bool   aggress_adjust;     ///< If true, allow a larger relative CFL change.
    bool   vol_CFL;            ///< If true, use the ratio of open face lengths to
                                ///< tracer cell areas when estimating CFL numbers.
    bool   better_iter;        ///< If true, use a velocity-based iteration criterion.
    bool   use_visc_rem_max;   ///< If true, use limiting bounds for viscous columns.
    bool   marginal_faces;     ///< If true, use marginal face areas as barotropic weights.
};

/// @brief Options controlling the edge-value reconstruction scheme used by the
/// continuity solver. Field-for-field mirror of the Fortran `bind(C)` type
/// `reconstruction_CS_C` -- order and types must not change.
struct reconstruction_CS_C {
    bool upwind_1st;  ///< If true, use a first-order upwind scheme.
    bool monotonic;   ///< If true, use the Colella & Woodward monotonic limiter.
    bool simple_2nd;  ///< If true, use a simple second order interpolation.
};

/// @brief AMReX ports of MOM6 numerical kernels.
namespace MOM {
using amrex::Box;
using amrex::Array4;
///  @brief Piecewise parabolic limiter
void ppm_limit_pos(const Box &,
                   Array4<const Real> const&,
                   Array4<Real> const&,
                   Array4<Real> const&,
                   const Real);

/// @brief Piecewise parabolic limiter of Colella and Woodward, 1984
void ppm_limit_cw84(const Box&,
                    Array4<const Real> const&,
                    Array4<Real> const&,
                    Array4<Real> const&);

/// @brief Piecewise reconstruction in the y dimension
void PPM_reconstruction_y(
    const Box&,
    Array4<const Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<const Real> const&,
    Real,
    bool,
    bool,
    OceanOBC*);

/// @brief Piecewise reconstruction in the x dimension
void PPM_reconstruction_x(
    const Box&,
    Array4<const Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<const Real> const&,
    Real,
    bool,
    bool,
    OceanOBC*);

/// @brief Zonal edge thickness — upwind copy or x-direction PPM reconstruction
void zonal_edge_thickness(
    const Box&,
    Array4<const Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<const Real> const&,
    Real,
    bool,
    bool,
    bool,
    OceanOBC*);

/// @brief Meridional edge thickness — upwind copy or y-direction PPM reconstruction
void meridional_edge_thickness(
    const Box&,
    Array4<const Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<const Real> const&,
    Real,
    bool,
    bool,
    bool,
    OceanOBC*);

/// @brief Zonal volume/thickness flux — PPM-reconstructed edge thickness
/// advected by the zonal velocity, scaled by viscosity remnant and
/// open-face area
void zonal_flux_thickness(
    const Box&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<Real> const&,
    Real,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    bool,
    bool,
    OceanOBC*,
    Array4<const Real> const&,
    Array4<const Real> const&);

/// @brief Meridional volume/thickness flux — PPM-reconstructed edge thickness
/// advected by the meridional velocity, scaled by viscosity remnant and
/// open-face area
void meridional_flux_thickness(
    const Box&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<Real> const&,
    Real,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    bool,
    bool,
    OceanOBC*,
    Array4<const Real> const&,
    Array4<const Real> const&);

/// @brief Zonal continuity update — advances layer thickness by the
/// convergence of the zonal thickness flux
void continuity_zonal_convergence(
    const Box&,
    Array4<Real> const&,
    Array4<const Real> const&,
    Real,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Real);

/// @brief Meridional continuity update — advances layer thickness by the
///  convergence of the meridional thickness flux
void continuity_meridional_convergence(
    const Box&,
    Array4<Real> const&,
    Array4<const Real> const&,
    Real,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Real);

/**
 * @brief Time steps the layer thicknesses, using a monotonically limit, directionally split PPM scheme,
 * based on Lin (1994).
 */
void continuity_PPM(
    Array4<const Real> const& u,                      //!< Zonal velocity [L T-1 ~> m s-1]
    Array4<const Real> const& v,                      //!< Meridional velocity [L T-1 ~> m s-1]
    Array4<const Real> const& hin,                    //!< Initial layer thickness [H ~> m or kg m-2]
    Array4<Real> const& h,                            //!< Final layer thickness [H ~> m or kg m-2]
    Array4<Real> const& uh,                           //!< Zonal volume flux, u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1]
    Array4<Real> const& vh,                           //!< Meridional volume flux, v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1]
    Real dt,                                          //!< Time increment [T ~> s]
    const Box& bx0,                                   //!< The core (unstencilled) iteration box
    int stencil,                                      //!< The continuity solver stencil width
    bool x_first,                                     //!< If true, advect zonally before meridionally
    Array4<const Real> const& mask2dT,                //!< Cell land/ocean mask [nondim]
    Array4<const Real> const& dy_Cu,                  //!< The grid cell's unblocked lengths of the u-faces of the
                                                      //!< h-cell [L ~> m]
    Array4<const Real> const& IareaT,                 //!< The grid cell's 1/areaT [L-2 ~> m-2]
    Array4<const Real> const& IdxT,                   //!< The grid cell's 1/dxT [L-1 ~> m-1]
    Array4<const Real> const& areaT,                  //!< The area of the h-cell [L2 ~> m2]
    Array4<const Real> const& dxT,                    //!< The x-extent of the h-cell [L ~> m]
    Array4<const Real> const& mask2dCu,               //!< 0 for land points, 1 for ocean points at u-locations [nondim]
    Array4<const Real> const& dxCu,                   //!< The grid cell's u-point x-extent [L ~> m]
    Array4<const Real> const& dx_Cv,                  //!< The grid cell's unblocked lengths of the v-faces of the
                                                      //!< h-cell [L ~> m]
    Array4<const Real> const& IdyT,                   //!< The grid cell's 1/dyT [L-1 ~> m-1]
    Array4<const Real> const& dyT,                    //!< The y-extent of the h-cell [L ~> m]
    Array4<const Real> const& mask2dCv,               //!< 0 for land points, 1 for ocean points at v-locations [nondim]
    Array4<const Real> const& dyCv,                   //!< The grid cell's v-point y-extent [L ~> m]
    int isd,                                          //!< Start i-index of the data domain (0-based)
    int ied,                                          //!< End i-index of the data domain (0-based)
    Real Angstrom_H,                                  //!< A one-Angstrom thickness [H ~> m or kg m-2]
    Real H_subroundoff,                               //!< A negligibly small thickness used to avoid division
                                                      //!< by zero [H ~> m or kg m-2]
    const reconstruction_CS_C& reconstruction_CS,     //!< Options controlling the edge-value reconstruction scheme
    const transport_adjust_CS_C& transport_adjust_CS, //!< Options controlling the transport adjustment and
                                                      //!< barotropic-consistency iteration
    OceanOBC* obc,                                    //!< Open boundary control structure
    Array4<const Real> const& por_face_areaU,         //!< Fractional open area of U-faces [nondim]
    Array4<const Real> const& por_face_areaV,         //!< Fractional open area of V-faces [nondim]
    Array4<const Real> const& uhbt,                   //!< Summed volume flux through zonal faces
                                                      //!< [H L2 T-1 ~> m3 s-1 or kg s-1];
                                                      //!< may be absent (.p == nullptr)
    Array4<const Real> const& vhbt,                   //!< Summed volume flux through meridional faces
                                                      //!< [H L2 T-1 ~> m3 s-1 or kg s-1];
                                                      //!< may be absent (.p == nullptr)
    Array4<const Real> const& visc_rem_u,             //!< Fraction of momentum/barotropic acceleration remaining
                                                      //!< after viscosity, zonal [nondim];
                                                      //!< may be absent (.p == nullptr)
    Array4<const Real> const& visc_rem_v,             //!< Fraction of momentum/barotropic acceleration remaining
                                                      //!< after viscosity, meridional [nondim];
                                                      //!< may be absent (.p == nullptr)
    Array4<Real> const& u_cor,                        //!< Zonal velocity with barotropic correction [L T-1 ~> m s-1];
                                                      //!< may be absent (.p == nullptr)
    Array4<Real> const& v_cor,                        //!< Meridional velocity with barotropic correction
                                                      //!< [L T-1 ~> m s-1];
                                                      //!< may be absent (.p == nullptr)
    Array4<Real> const& FA_u_W0,                      //!< Effective open face area, west, 0 transport;
                                                      //!< may be absent (.p == nullptr)
    Array4<Real> const& FA_u_E0,                      //!< Effective open face area, east, 0 transport;
                                                      //!< may be absent (.p == nullptr)
    Array4<Real> const& FA_u_WW,                      //!< Effective open face area, westerly test velocity;
                                                      //!< may be absent (.p == nullptr)
    Array4<Real> const& FA_u_EE,                      //!< Effective open face area, easterly test velocity;
                                                      //!< may be absent (.p == nullptr)
    Array4<Real> const& uBT_WW,                       //!< Westerly correction to the barotropic velocity;
                                                      //!< may be absent (.p == nullptr)
    Array4<Real> const& uBT_EE,                       //!< Easterly correction to the barotropic velocity;
                                                      //!< may be absent (.p == nullptr)
    Array4<Real> const& FA_v_S0,                      //!< Effective open face area, south, 0 transport;
                                                      //!< may be absent (.p == nullptr)
    Array4<Real> const& FA_v_N0,                      //!< Effective open face area, north, 0 transport;
                                                      //!< may be absent (.p == nullptr)
    Array4<Real> const& FA_v_SS,                      //!< Effective open face area, southerly test velocity;
                                                      //!< may be absent (.p == nullptr)
    Array4<Real> const& FA_v_NN,                      //!< Effective open face area, northerly test velocity;
                                                      //!< may be absent (.p == nullptr)
    Array4<Real> const& vBT_SS,                       //!< Southerly correction to the barotropic velocity;
                                                      //!< may be absent (.p == nullptr)
    Array4<Real> const& vBT_NN,                       //!< Northerly correction to the barotropic velocity;
                                                      //!< may be absent (.p == nullptr)
    Array4<Real> const& h_u,                          //!< Effective thickness at zonal faces [H ~> m or kg m-2];
                                                      //!< may be absent (.p == nullptr)
    Array4<Real> const& h_v,                          //!< Effective thickness at meridional faces [H ~> m or kg m-2];
                                                      //!< may be absent (.p == nullptr)
    Array4<Real> const& du_cor,                       //!< Zonal velocity increment from u that gives uhbt as the
                                                      //!< depth-integrated transport [L T-1 ~> m s-1];
                                                      //!< may be absent (.p == nullptr)
    Array4<Real> const& dv_cor);                      //!< Meridional velocity increment from v that gives vhbt as the
                                                      //!< depth-integrated transport [L T-1 ~> m s-1];
                                                      //!< may be absent (.p == nullptr)

/**
 * @brief Find the vertical sum of the thickness fluxes from the continuity solver without actually
 * updating the layer thicknesses.  Because the fluxes in the two directions are calculated
 * based on the input thicknesses, which are not updated between the directions, the fluxes
 * returned here are not the same as those that would be returned by a call to continuity.
 */
void continuity_PPM_2d_fluxes(
    Array4<const Real> const& u,                      //!< Zonal velocity [L T-1 ~> m s-1]
    Array4<const Real> const& v,                      //!< Meridional velocity [L T-1 ~> m s-1]
    Array4<const Real> const& h,                      //!< Layer thickness [H ~> m or kg m-2]
    Array4<Real> const& uhbt,                         //!< Vertically summed thickness flux through zonal
                                                      //!< faces [H L2 T-1 ~> m3 s-1 or kg s-1]
    Array4<Real> const& vhbt,                         //!< Vertically summed thickness flux through meridional
                                                      //!< faces [H L2 T-1 ~> m3 s-1 or kg s-1]
    Real dt,                                          //!< Time increment [T ~> s]
    const Box& bxC,                                   //!< Iteration box for continuity solver
    Array4<const Real> const& mask2dT,                //!< Cell land/ocean mask [nondim]
    Array4<const Real> const& dy_Cu,                  //!< The grid cell's unblocked lengths of the u-faces of the
                                                      //!< h-cell [L ~> m]
    Array4<const Real> const& IareaT,                 //!< The grid cell's 1/areaT [L-2 ~> m-2]
    Array4<const Real> const& IdxT,                   //!< The grid cell's 1/dxT [L-1 ~> m-1]
    Array4<const Real> const& dx_Cv,                  //!< The grid cell's unblocked lengths of the v-faces of the
                                                      //!< h-cell [L ~> m]
    Array4<const Real> const& IdyT,                   //!< The grid cell's 1/dyT [L-1 ~> m-1]
    Real Angstrom_H,                                  //!< A one-Angstrom thickness [H ~> m or kg m-2]
    const reconstruction_CS_C& reconstruction_CS,     //!< Options controlling the edge-value reconstruction scheme
    const transport_adjust_CS_C& transport_adjust_CS, //!< Options controlling the transport adjustment and
                                                      //!< barotropic-consistency iteration
    OceanOBC* obc,                                    //!< Open boundary control structure
    Array4<const Real> const& por_face_areaU,         //!< Fractional open area of U-faces [nondim]
    Array4<const Real> const& por_face_areaV);        //!< Fractional open area of V-faces [nondim]

/**
 * @brief Calculates the mass or volume fluxes through the zonal faces, and
 * other related quantities -- including, optionally, the barotropic
 * mass-flux correction (u_cor/du_cor) and the barotropic-consistency
 * face-area/velocity-correction diagnostics (FA_u_W0/E0/WW/EE, uBT_WW/EE)
 */
void zonal_mass_flux(
    const Box&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<Real> const&,
    Real,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Real,
    const transport_adjust_CS_C&,
    OceanOBC*,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<Real> const&);

/**
 * @brief Calculates the mass or volume fluxes through the meridional faces,
 * and other related quantities -- including, optionally, the barotropic
 * mass-flux correction (v_cor/dv_cor) and the barotropic-consistency
 * face-area/velocity-correction diagnostics (FA_v_S0/N0/SS/NN, vBT_SS/NN)
 */
void meridional_mass_flux(
    const Box&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<Real> const&,
    Real,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    int,
    int,
    Real,
    const transport_adjust_CS_C&,
    OceanOBC*,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<Real> const&);

/**
 * @brief Accumulates the vertically-summed zonal barotropic mass/volume
 * transport across the water column, for use as the barotropic solver's
 * target transport in the transport-adjustment iteration
 */
void zonal_BT_mass_flux(
    const Box&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<Real> const&,
    Real,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    const transport_adjust_CS_C&,
    OceanOBC*,
    Array4<const Real> const&);

/**
 * @brief Accumulates the vertically-summed meridional barotropic mass/volume
 * transport across the water column, for use as the barotropic solver's
 * target transport in the transport-adjustment iteration
 */
void meridional_BT_mass_flux(
    const Box&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<Real> const&,
    Real,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    const transport_adjust_CS_C&,
    OceanOBC*,
    Array4<const Real> const&);

/**
 * @brief Newton-iterates a barotropic velocity correction per zonal face so
 * that the vertically-summed zonal mass/volume transport matches the target
 * barotropic transport, to within the transport-adjustment iteration's
 * tolerance
 */
void zonal_flux_adjust(
    const Box&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Real,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    const transport_adjust_CS_C&,
    Array4<const Real> const&,
    Array4<const int> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<Real> const&,
    OceanOBC*);

/**
 * @brief Newton-iterates a barotropic velocity correction per meridional
 * face so that the vertically-summed meridional mass/volume transport
 * matches the target barotropic transport, to within the
 * transport-adjustment iteration's tolerance
 */
void meridional_flux_adjust(
    const Box&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Real,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    const transport_adjust_CS_C&,
    Array4<const Real> const&,
    Array4<const int> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<Real> const&,
    OceanOBC*);

/**
 * @brief Sets the effective open face areas and barotropic velocity
 * corrections at zonal faces that reproduce the summed layer
 * transports for three test barotropic velocities, for use in the
 * barotropic-consistency iteration
 */
void set_zonal_BT_cont(
    const Box&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<const Real> const&,
    Real,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    const transport_adjust_CS_C&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const int> const&,
    Array4<const Real> const&);

/**
 * @brief Sets the effective open face areas and barotropic velocity
 * corrections at meridional faces that reproduce the summed layer
 * transports for three test barotropic velocities, for use in the
 * barotropic-consistency iteration
 */
void set_merid_BT_cont(
    const Box&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<const Real> const&,
    Real,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    const transport_adjust_CS_C&,
    Array4<const Real> const&,
    Array4<const Real> const&,
    Array4<const int> const&,
    Array4<const Real> const&);

}
