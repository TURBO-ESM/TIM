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
    // The Fortran kernel uses ((.not.Boussinesq) .and. allocated(tv%SpV_avg))
    // as its branch predicate. The Fortran shim has already evaluated this
    // and exposes it as has_spv; preserve the explicit conjunction here so
    // the C++ branch matches the Fortran one literally.
    const bool use_spv = (!boussinesq) && has_spv;

    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        thickness_to_dz_3d_point(dz(i,j,k),
                                 h(i,j,k),
                                 use_spv ? spv_avg(i,j,k) : Real(0),
                                 h_to_z,
                                 h_to_rz,
                                 use_spv);
    });
}

}
