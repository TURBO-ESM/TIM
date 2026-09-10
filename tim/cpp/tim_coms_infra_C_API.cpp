/**
 * @file tim_coms_infra_C_API.cpp
 * @brief C-API translation unit of the TIM checksum service.
 */

#include "tim_coms_infra_C_API.h"
#include "tim_coms_infra.hpp"

int64_t tim_chksum_c(const RealArray_C* field_HOST, double* mask_val, bool global_chksum)
{
    amrex::Box bx({0, 0, 0},
                  {field_HOST->shape[0], field_HOST->shape[1], field_HOST->shape[2]});

    size_t num_points = bx.numPts();

    amrex::Real* device_array = static_cast<amrex::Real*>(TheArena()->alloc(num_points * sizeof(amrex::Real)));
    Gpu::copy(Gpu::hostToDevice, field_HOST->data, field_HOST.data+num_points, device_array);

    amrex::BaseFab<amrex::Real> non_owning_fab(bx, num_points, device_array);
    auto array_1d = non_owning_fab.array();

    ///-------------------------------------------------
    /// Execute checksum
    ///-------------------------------------------------
    int64_t chksum = mask_val ? TIM::checksum(bx, array_1d, *mask_val)
                              : TIM::checksum(bx, array_1d);

    TheArena()->free(device_array);

    return chksum;
}
