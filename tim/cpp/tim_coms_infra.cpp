/**
 * @file tim_coms_infra.cpp
 * @brief Implementation of the TIM checksum service.
 */

// #include <bit>

#include <AMReX_Gpu.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Reduce.H>

#include "tim_coms_infra.hpp"

namespace TIM {

int64_t checksum(double* field, size_t field_size, double* mask_val)
{
    // Cast needed due to lack of C++20 support in icpc needed to do:
    //   auto checksum = std::bit_cast<int64_t>(var);
    // See https://en.wikipedia.org/wiki/Fast_inverse_square_root#Overview_of_the_code
    // for logic behind casting.
    // Also need C++20 to use std::view::filter for mask.
    amrex::Long checksum = 0;
    if(mask_val)
    {
        amrex::Long mask_bytes = * ( amrex::Long * ) mask_val;
        checksum = amrex::Reduce::Sum<amrex::Long>(field_size,
            [=] AMREX_GPU_DEVICE (size_t i) -> amrex::Long {
                amrex::Long field_checksum = * ( amrex::Long * ) &field[i];
                return (field_checksum == mask_bytes) ? 0 : field_checksum;
            });
    }
    else
        checksum = amrex::Reduce::Sum<amrex::Long>(field_size,
            [=] AMREX_GPU_DEVICE (size_t i) -> amrex::Long {
                return * ( amrex::Long * ) &field[i];
            });
    amrex::ParallelDescriptor::ReduceLongSum(checksum);
    return checksum;
}

} // namespace TIM
