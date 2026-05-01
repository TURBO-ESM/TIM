#pragma once

#include <AMReX_Box.H>
#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

namespace MOM {
using amrex::Real;

/**
 * @brief Per-cell conversion of thickness to geometric layer thickness.
 *
 *  Direct port of the inner loop body of MOM6
 *  MOM_interface_heights::thickness_to_dz_3d_fortran. Selects between the
 *  Boussinesq form (multiply by H_to_Z) and the non-Boussinesq form
 *  (multiply by H_to_RZ * SpV_avg) using the same predicate the Fortran
 *  shim sets: has_spv == ((.not.Boussinesq) .and. allocated(tv%SpV_avg)).
 *
 *  @param dz       [out] Geometric layer thickness [Z ~> m]
 *  @param h        [in]  Layer thickness in thickness units [H ~> m or kg m-2]
 *  @param spv      [in]  Layer-mean specific volume; only read when
 *                        has_spv is true
 *  @param h_to_z   [in]  Conversion factor from H to Z [Z H-1]
 *  @param h_to_rz  [in]  Conversion factor from H to R*Z [R Z H-1]
 *  @param has_spv  [in]  True iff non-Boussinesq AND tv%SpV_avg is allocated
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void thickness_to_dz_3d_point(Real& dz,
                              Real const h,
                              Real const spv,
                              Real const h_to_z,
                              Real const h_to_rz,
                              bool const has_spv) noexcept
{
    if (has_spv) {
        dz = h_to_rz * h * spv;
    } else {
        dz = h_to_z * h;
    }
}

}
