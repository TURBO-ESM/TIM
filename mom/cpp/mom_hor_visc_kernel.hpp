// SKILLS: 0.3.3
#pragma once
///
/// @file mom_hor_visc_kernel.hpp
/// @brief Per-point device primitives for the MOM6 horizontal viscosity kernels.
///
#include <AMReX_Algorithm.H>
#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

#include <cmath>

namespace MOM {
using amrex::Real;

///
/// @brief Adds the MEKE-based (EY24_EBT_BS) backscatter contribution to the diagonal term of
///        the stress tensor at one h-point.
///
/// @tparam use_kh_struct     If true, shape the backscatter coefficient with VarMix%BS_struct
/// @param str_xx             The diagonal term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2]
/// @param Kh_BS              Laplacian antiviscosity [L2 T-1 ~> m2 s-1]
/// @param visc_limit_h_flag  determines whether backscatter is shut off [nondim]
/// @param MEKE_Ku            The MEKE-derived lateral viscosity coefficient [L2 T-1 ~> m2 s-1]
/// @param BS_struct          Vertical structure function used in backscatter [nondim]
/// @param sh_xx              horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
///
template <bool use_kh_struct>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void hor_visc_backscatter_h_point(Real& str_xx, Real& Kh_BS,
                                  Real const visc_limit_h_flag,
                                  Real const MEKE_Ku,
                                  Real const BS_struct,
                                  Real const sh_xx) noexcept
{
    using namespace amrex::literals;

    if (visc_limit_h_flag > 0.0_rt) {
        Kh_BS = 0._rt;
    } else {
        if constexpr (use_kh_struct) {
            Kh_BS = MEKE_Ku * BS_struct;
        } else {
            Kh_BS = MEKE_Ku;
        }
    }

    /// The diagonal term in the stress tensor due to backscatter [H L2 T-2 ~> m3 s-2 or kg s-2]
    Real const str_xx_BS = -Kh_BS * sh_xx;

    str_xx = str_xx + str_xx_BS;
}

///
/// @brief Sets the Leith+E biharmonic viscosity Ah and Ah_h at one h-point from the smoothed
///        square of Ah (CS%smooth_Ah true).
///
/// @param Ah        biharmonic viscosity (h or q) [L4 T-1 ~> m4 s-1]
/// @param Ah_h      biharmonic viscosity at thickness points [L4 T-1 ~> m4 s-1]
/// @param Ah_bg_xx  The background biharmonic viscosity at h points [L4 T-1 ~> m4 s-1]
/// @param Ah_sq     The square of the biharmonic viscosity [L8 T-2 ~> m8 s-2]
///
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void hor_visc_leithy_ah_h_smooth_point(Real& Ah, Real& Ah_h,
                                       Real const Ah_bg_xx,
                                       Real const Ah_sq) noexcept
{
    using namespace amrex::literals;

    Ah_h = amrex::max(Ah_bg_xx, std::sqrt(amrex::max(0._rt, Ah_sq)));
    Ah   = Ah_h;
}

///
/// @brief Sets the biharmonic viscosity Ah_h at one h-point from Ah (CS%smooth_Ah false).
///
/// @param Ah_h  biharmonic viscosity at thickness points [L4 T-1 ~> m4 s-1]
/// @param Ah    biharmonic viscosity (h or q) [L4 T-1 ~> m4 s-1]
///
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void hor_visc_leithy_ah_h_nosmooth_point(Real& Ah_h,
                                         Real const Ah) noexcept
{
    Ah_h = Ah;
}

} // namespace MOM
