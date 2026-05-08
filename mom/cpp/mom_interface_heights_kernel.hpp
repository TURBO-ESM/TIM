#pragma once

#include <AMReX_Box.H>
#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

namespace MOM {
using amrex::Real;

/**
 * @brief Per-cell conversion (Boussinesq path): dz = H_to_Z * h.
 *
 *  @param dz      [out] Geometric layer thickness [Z ~> m]
 *  @param h       [in]  Layer thickness [H ~> m or kg m-2]
 *  @param h_to_z  [in]  Conversion factor H → Z [Z H-1]
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void thickness_to_dz_3d_boussinesq_point(Real& dz,
                                         Real const h,
                                         Real const h_to_z) noexcept
{
    dz = h_to_z * h;
}

/**
 * @brief Per-cell conversion (non-Boussinesq path): dz = H_to_RZ * h * SpV_avg.
 *
 *  Called only when has_spv is true (i.e. !Boussinesq && allocated(tv%SpV_avg)).
 *
 *  @param dz      [out] Geometric layer thickness [Z ~> m]
 *  @param h       [in]  Layer thickness [H ~> m or kg m-2]
 *  @param spv     [in]  Layer-mean specific volume [R-1 ~> m3 kg-1]
 *  @param h_to_rz [in]  Conversion factor H → R*Z [R Z H-1]
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void thickness_to_dz_3d_nonboussinesq_point(Real& dz,
                                            Real const h,
                                            Real const spv,
                                            Real const h_to_rz) noexcept
{
    dz = h_to_rz * h * spv;
}

}
