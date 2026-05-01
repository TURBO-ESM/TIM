// turbotmp_mom_interface_heights_bridge.cpp
//
// Marshalling bridge between the MOM6 Fortran shim
// MOM_interface_heights::thickness_to_dz_3d (mode = AMREX) and the AMReX
// kernel MOM::thickness_to_dz_3d. Pure host↔device transfer, layout
// transpose, and Box reconstruction — no math, no value-dependent
// branches (lessons.md §3).

#include "mom_interface_heights.hpp"
#include "turbotmp_helper.hpp"
#include "turbotmp_mom_interface_heights_bridge.h"

using namespace amrex;

/**
 * @brief Bridge for the MOM6 thickness_to_dz_3d subroutine.
 *
 * @param bx_HOST       Iteration domain (Fortran 1-based indices).
 * @param h_HOST        Input layer thickness [H ~> m or kg m-2]
 *                      (host, Fortran order).
 * @param dz_HOST       Output geometric layer thickness [Z ~> m]
 *                      (host, Fortran order). Declared inout in Fortran;
 *                      the kernel writes every cell in the iteration Box,
 *                      but the bridge still copies it host->device on entry
 *                      to preserve any halo values the caller relies on.
 * @param spv_avg_HOST  Layer-mean specific volume (host, Fortran order).
 *                      Read only when @p has_spv is true. The Fortran shim
 *                      always supplies a valid descriptor (a 1x1x1 placeholder
 *                      when has_spv is false), so the bridge can copy it
 *                      unconditionally.
 * @param boussinesq    True if running in Boussinesq mode.
 * @param h_to_z        Conversion factor from H to Z [Z H-1].
 * @param h_to_rz       Conversion factor from H to R*Z [R Z H-1].
 * @param has_spv       True iff non-Boussinesq AND tv%SpV_avg is allocated.
 */
void turbotmp_thickness_to_dz_3d_bridge(const Box_C* bx_HOST,
                                        const RealArray_C* h_HOST,
                                        RealArray_C* dz_HOST,
                                        const RealArray_C* spv_avg_HOST,
                                        const bool boussinesq,
                                        const double h_to_z,
                                        const double h_to_rz,
                                        const bool has_spv)
{
    /// Reconstruct the iteration Box (Fortran 1-based -> AMReX 0-based)
    Box bx(IntVect(bx_HOST->idxS[0]-1, bx_HOST->idxS[1]-1, bx_HOST->idxS[2]-1),
           IntVect(bx_HOST->idxE[0]-1, bx_HOST->idxE[1]-1, bx_HOST->idxE[2]-1));

    /// Allocate one A4Box per array, sized by RealArray_C::shape[]
    auto h_DEV       = turbotmp::make_array4(h_HOST->shape[0],
                                             h_HOST->shape[1],
                                             h_HOST->shape[2], 1);
    auto dz_DEV      = turbotmp::make_array4(dz_HOST->shape[0],
                                             dz_HOST->shape[1],
                                             dz_HOST->shape[2], 1);
    auto spv_avg_DEV = turbotmp::make_array4(spv_avg_HOST->shape[0],
                                             spv_avg_HOST->shape[1],
                                             spv_avg_HOST->shape[2], 1);

    /// Copy host -> device for every array, including the inout dz so its
    /// halo values are preserved across the call (lessons.md §7 #1).
    turbotmp::copy_FortranHost_to_array4(h_HOST->data,       h_DEV);
    turbotmp::copy_FortranHost_to_array4(dz_HOST->data,      dz_DEV);
    turbotmp::copy_FortranHost_to_array4(spv_avg_HOST->data, spv_avg_DEV);

    ///-------------------------------------------------
    /// Execute kernel
    ///-------------------------------------------------
    MOM::thickness_to_dz_3d(bx,
                            h_DEV.arr,
                            dz_DEV.arr,
                            spv_avg_DEV.arr,
                            boussinesq,
                            h_to_z,
                            h_to_rz,
                            has_spv);

    /// Ensure kernel is done before copying back
    Gpu::synchronize();

    /// Copy device -> host for outputs/inouts only
    turbotmp::copy_array4_to_FortranHost(dz_DEV, dz_HOST->data);

    /// Free everything we allocated
    turbotmp::free_array4(h_DEV);
    turbotmp::free_array4(dz_DEV);
    turbotmp::free_array4(spv_avg_DEV);
}
