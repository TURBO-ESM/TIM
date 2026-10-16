#include <AMReX_Arena.H>
#include <AMReX_Gpu.H>
#include <AMReX_Array4.H>
#include <AMReX_MultiFab.H>

#include "turbotmp_helper.hpp"

namespace turbotmp {

namespace {
// Extents and lower bounds of a 2D or 3D descriptor, padded to 3D with 1.
void descriptor_dims(int dim, const int* shape, const int* lb, int n[3], int l[3])
{
    if (dim != 2 && dim != 3) {
        amrex::Abort("turbotmp: array descriptor must be 2D or 3D");
    }
    for (int d = 0; d < 3; ++d) {
        n[d] = (d < dim) ? shape[d] : 1;
        l[d] = (d < dim) ? lb[d]    : 1;
    }
}
}

A4Box make_array4(const RealArray_C* a)
{
    using namespace amrex;

    int n[3], l[3];
    descriptor_dims(a->dim, a->shape, a->lb, n, l);
    const int nx = n[0], ny = n[1], nz = n[2];
    const int lbx = l[0], lby = l[1], lbz = l[2];
    const int ncomp = 1;

    A4Box a4;

    a4.nx = nx;
    a4.ny = ny;
    a4.nz = nz;
    a4.ncomp = ncomp;

    // lbx/lby/lbz are Fortran index-space lower bounds; position the box at the
    // array's true absolute location (matching the bridge's own -1 conversion
    // for Box_C), not just [0, shape-1] -- staggered-grid arrays legitimately
    // have lb != 1.
    a4.bx = Box(IntVect(lbx-1, lby-1, lbz-1),
               IntVect(lbx-1 + nx-1, lby-1 + ny-1, lbz-1 + nz-1));

    const Long npts = ncomp*a4.bx.numPts();

    // Allocate in AMReX arena (GPU-aware)
    a4.data   = (Real*) The_Arena()->alloc(npts * sizeof(Real));
    a4.data_f = (Real*) The_Arena()->alloc(npts * sizeof(Real));

    // setup AMReX views
    a4.arr = Array4<Real>(a4.data, begin(a4.bx), end(a4.bx),ncomp);

    return a4;
}

void free_array4(A4Box& a4)
{
    using namespace amrex;

    if (a4.data)
    {
        The_Arena()->free(a4.data);
        a4.data = nullptr;
    }
    if(a4.data_f) {
	The_Arena()->free(a4.data_f);
	a4.data_f = nullptr;
    }
}

void copy_FortranHost_to_array4(const double* f, A4Box& a4)
{
    using namespace amrex;

    int nx = a4.nx;
    int ny = a4.ny;
    int nz = a4.nz;
    int ncomp = a4.ncomp;

    auto arr = a4.arr;
    Real* d_f = a4.data_f;
    Box bx = a4.bx;

    auto lo = lbound(bx);

    // Totl number of elements
    Long n = static_cast<Long>(nx) * ny * nz * ncomp;


    // Copy host -> device
    Gpu::copy(Gpu::hostToDevice, f, f+n, d_f);


    // Transpose array from Fortran to C
    ParallelFor(bx, ncomp,  [=] AMREX_GPU_DEVICE (int i, int j, int k, int nc)
    {
        int ii = i - lo.x;
	int jj = j - lo.y;
	int kk = k - lo.z;
        int idx = ii + nx*(jj + ny*(kk + nz*nc));  // Fortran layout with component stride
        arr(i,j,k,nc) = d_f[idx];
    });
    Gpu::streamSynchronize();
}

void copy_array4_to_FortranHost(const A4Box& a4, double* f)
{
    using namespace amrex;

    int nx = a4.nx;
    int ny = a4.ny;
    int nz = a4.nz;
    int ncomp = a4.ncomp;

    auto arr = a4.arr;
    Real* d_f = a4.data_f;
    Box bx = a4.bx;

    auto lo = lbound(bx);

    Long n = static_cast<Long>(nx) * ny * nz * ncomp;

    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int nc)
    {
        int ii = i - lo.x;
	int jj = j - lo.y;
	int kk = k - lo.z;
        int idx = ii + nx*(jj + ny*(kk + nz*nc));  // Fortran layout with component stride
        d_f[idx] = arr(i,j,k,nc);
    });


    // Copy device -> Host
    Gpu::copy(Gpu::deviceToHost, d_f, d_f+n,f);
    Gpu::streamSynchronize();

}

IntA4Box make_int_array4(const LogicalArray_C* a)
{
    using namespace amrex;

    int n[3], l[3];
    descriptor_dims(a->dim, a->shape, a->lb, n, l);
    const int nx = n[0], ny = n[1], nz = n[2];
    const int lbx = l[0], lby = l[1], lbz = l[2];
    const int ncomp = 1;

    IntA4Box a4;

    a4.nx = nx;
    a4.ny = ny;
    a4.nz = nz;
    a4.ncomp = ncomp;

    //  See make_array() comment -- identical convention
    a4.bx = Box(IntVect(lbx-1, lby-1, lbz-1),
               IntVect(lbx-1 + nx-1, lby-1 + ny-1, lbz-1 + nz-1));

    const Long npts = ncomp*a4.bx.numPts();

    // Allocate in AMReX arena (GPU-aware)
    a4.data   = (int*) The_Arena()->alloc(npts * sizeof(int));
    a4.data_f = (int*) The_Arena()->alloc(npts * sizeof(int));

    // setup AMReX views
    a4.arr = Array4<int>(a4.data, begin(a4.bx), end(a4.bx), ncomp);

    return a4;
}

void free_int_array4(IntA4Box& a4)
{
    using namespace amrex;

    if (a4.data)
    {
        The_Arena()->free(a4.data);
        a4.data = nullptr;
    }
    if (a4.data_f)
    {
        The_Arena()->free(a4.data_f);
        a4.data_f = nullptr;
    }
}

void copy_FortranHost_to_int_array4(const int* f, IntA4Box& a4)
{
    using namespace amrex;

    int nx = a4.nx;
    int ny = a4.ny;
    int nz = a4.nz;
    int ncomp = a4.ncomp;

    auto arr = a4.arr;
    int* d_f = a4.data_f;
    Box bx = a4.bx;

    auto lo = lbound(bx);

    // Total number of elements
    Long n = static_cast<Long>(nx) * ny * nz * ncomp;

    // Copy host -> device
    Gpu::copy(Gpu::hostToDevice, f, f+n, d_f);

    // Transpose array from Fortran to C
    ParallelFor(bx, ncomp,  [=] AMREX_GPU_DEVICE (int i, int j, int k, int nc)
    {
        int ii = i - lo.x;
	int jj = j - lo.y;
	int kk = k - lo.z;
        int idx = ii + nx*(jj + ny*(kk + nz*nc));  // Fortran layout with component stride
        arr(i,j,k,nc) = d_f[idx];
    });
    Gpu::streamSynchronize();
}

void copy_int_array4_to_FortranHost(const IntA4Box& a4, int* f)
{
    using namespace amrex;

    int nx = a4.nx;
    int ny = a4.ny;
    int nz = a4.nz;
    int ncomp = a4.ncomp;

    auto arr = a4.arr;
    int* d_f = a4.data_f;
    Box bx = a4.bx;

    auto lo = lbound(bx);

    Long n = static_cast<Long>(nx) * ny * nz * ncomp;

    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int nc)
    {
        int ii = i - lo.x;
	int jj = j - lo.y;
	int kk = k - lo.z;
        int idx = ii + nx*(jj + ny*(kk + nz*nc));  // Fortran layout with component stride
        d_f[idx] = arr(i,j,k,nc);
    });

    // Copy device -> Host
    Gpu::copy(Gpu::deviceToHost, d_f, d_f+n, f);
    Gpu::streamSynchronize();
}

amrex::Box to_box(const Box_C* b)
{
    return amrex::Box(amrex::IntVect(b->idxS[0]-1, b->idxS[1]-1, b->idxS[2]-1),
                      amrex::IntVect(b->idxE[0]-1, b->idxE[1]-1, b->idxE[2]-1));
}

}
