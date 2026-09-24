/**
 * @file tim_coms_infra.cpp
 * @brief Implementation of the TIM checksum service.
 */

// #include <bit>

#include <AMReX_Gpu.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Reduce.H>

#include "tim_coms_infra.hpp"

namespace {

// Per-rank checksum of a single box: the piece shared by every public
// overload, kept apart from the cross-rank reduction so callers that need to
// sum several boxes (e.g. a MultiFab's local boxes) can do so with a single
// ParallelDescriptor::ReduceLongSum instead of one per box.
amrex::Long local_checksum(amrex::Box const& bx, amrex::Array4<const amrex::Real> const& arr,
                           std::optional<amrex::Real> mask)
{
    amrex::Reducer<amrex::ReduceOpSum, amrex::Long> reducer;
    using Result_t = typename decltype(reducer)::Result_t;
    const int ncomp = arr.nComp();
    
    if (mask) {
        amrex::Long mask_bytes;
        std::memcpy(&mask_bytes, &(*mask), sizeof(mask_bytes));

        reducer.eval(bx, ncomp
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) -> Result_t
            {
                amrex::Long bits;
                std::memcpy(&bits, &arr(i,j,k, n), sizeof(bits));
                return { (bits == mask_bytes) ? 0 : bits };
            });
    } else {
        reducer.eval(bx, ncomp
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) -> Result_t
            {
                amrex::Long bits;
                std::memcpy(&bits, &arr(i,j,k,n), sizeof(bits));
                return bits;
            });
    }
    return amrex::get<0>(reducer.getResult());
}

} // namespace

namespace TIM {

amrex::Long checksum(amrex::Box const& bx,
                     amrex::Array4<const amrex::Real> const& arr,
                     std::optional<amrex::Real> mask)
{
    amrex::Long checksum = local_checksum(bx, arr, mask);
    amrex::ParallelDescriptor::ReduceLongSum(checksum);
    return checksum;
}

amrex::Long checksum(amrex::MultiFab const& mf, std::optional<amrex::Real> mask)
{
    amrex::Long checksum = 0;
    for (amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        checksum += local_checksum(mfi.growntilebox(), mf.const_array(mfi), mask);
    }
    amrex::ParallelDescriptor::ReduceLongSum(checksum);
    return checksum;
}

} // namespace TIM
