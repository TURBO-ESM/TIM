// mom_interface_heights.hpp
#pragma once

#include "mom_interface_heights_kernel.hpp"

namespace MOM {
using amrex::Box;
using amrex::Array4;

/**
 * @brief Convert layer thickness in thickness units (H) to geometric
 *        layer thickness in height units (Z).
 *
 *  Box-level AMReX kernel for the MOM6
 *  MOM_interface_heights::thickness_to_dz_3d subroutine. Checks the
 *  Boussinesq/non-Boussinesq mode once, then dispatches to either
 *  thickness_to_dz_3d_point_boussinesq or
 *  thickness_to_dz_3d_point_nonboussinesq for every cell — no per-cell
 *  branch.
 *
 *  @param bx          Iteration domain (Fortran 1-based indices already
 *                     converted to 0-based by the bridge).
 *  @param h           Input layer thickness [H ~> m or kg m-2].
 *  @param dz          Output geometric layer thickness [Z ~> m].
 *  @param spv_avg     Layer-mean specific volume; read only in the
 *                     non-Boussinesq path (has_spv true).
 *  @param boussinesq  True if running in Boussinesq mode.
 *  @param h_to_z      Conversion factor from H to Z [Z H-1].
 *  @param h_to_rz     Conversion factor from H to R*Z [R Z H-1].
 *  @param has_spv     True iff non-Boussinesq AND tv%SpV_avg is allocated.
 */
void thickness_to_dz_3d(const Box& bx,
                        Array4<const Real> const& h,
                        Array4<Real>       const& dz,
                        Array4<const Real> const& spv_avg,
                        const bool boussinesq,
                        const Real h_to_z,
                        const Real h_to_rz,
                        const bool has_spv);
}
