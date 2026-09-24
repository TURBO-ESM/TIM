/**
 * @file tim_coms_infra_C_API.cpp
 * @brief C-API translation unit of the TIM checksum service.
 */

#include "tim_coms_infra_C_API.h"
#include "tim_coms_infra.hpp"

#include <AMReX_Arena.H>
#include <AMReX_Gpu.H>

int64_t tim_chksum_c(const RealArray_C* field_HOST, double* mask_val)
{
    // Default box dimmensions to 1 unless shape overrides.
    int dim_size[3] = {1, 1, 1};
    for (int d = 0; d < field_HOST->dim && d < 3; ++d)
      dim_size[d] = field_HOST->shape[d];

    // Default 4th dimmension to 1 element unless in 4D case
    int ncomp = 1;
    if (field_HOST->dim == 4)
      ncomp = field_HOST->shape[3];

    amrex::Box bx({0, 0, 0},
                  {dim_size[0]-1, dim_size[1]-1, dim_size[2]-1});

    size_t num_points = bx.numPts() * ncomp;

    amrex::Real* device_array = static_cast<amrex::Real*>(amrex::The_Arena()->alloc(num_points * sizeof(amrex::Real)));
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, field_HOST->data, field_HOST->data+num_points, device_array);

    amrex::BaseFab<amrex::Real> non_owning_fab(bx, ncomp, device_array);
    auto array4 = non_owning_fab.array();

    ///-------------------------------------------------
    /// Execute checksum
    ///-------------------------------------------------
    int64_t chksum = mask_val ? TIM::checksum(bx, array4, *mask_val)
                              : TIM::checksum(bx, array4);

    amrex::The_Arena()->free(device_array);

    return chksum;
}
