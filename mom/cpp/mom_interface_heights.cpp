// mom_interface_heights.cpp
#include "mom_interface_heights.hpp"

using namespace amrex;

namespace MOM {

void thickness_to_dz_3d(const Box& bx,
                        Array4<const Real> const& h,
                        Array4<Real>       const& dz,
                        Array4<const Real> const& spv_avg,
                        const bool boussinesq,
                        const Real h_to_z,
                        const Real h_to_rz,
                        const bool has_spv)
{
    if ((!boussinesq) && has_spv) {
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            thickness_to_dz_3d_nonboussinesq_point(dz(i,j,k),
                                                   h(i,j,k),
                                                   spv_avg(i,j,k),
                                                   h_to_rz);
        });
    } else {
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            thickness_to_dz_3d_boussinesq_point(dz(i,j,k),
                                                h(i,j,k),
                                                h_to_z);
        });
    }
}

}
