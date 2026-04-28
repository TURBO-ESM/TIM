#include <AMReX_Arena.H>
#include <AMReX_Gpu.H>
#include <AMReX_Array4.H>
#include <AMReX_MultiFab.H>

#include "turbotmp_helper.hpp"

namespace turbotmp {

A4Box make_array4(int nx, int ny, int nz, int ncomp)
{
    using namespace amrex;

    A4Box a4;

    a4.nx = nx;
    a4.ny = ny;
    a4.nz = nz;
    a4.ncomp = ncomp;

    a4.bx = Box(IntVect(0,0,0), IntVect(nx-1, ny-1, nz-1));

    const Long npts = ncomp*a4.bx.numPts();

    // Allocate in AMReX arena (GPU-aware)
    a4.data   = (Real*) The_Arena()->alloc(npts * sizeof(Real));
    a4.data_f = (Real*) The_Arena()->alloc(npts * sizeof(Real));

    // setup AMReX views
    a4.arr = Array4<Real>(a4.data, lbound(a4.bx), ubound(a4.bx),ncomp);

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


    // Tranpose array from Fortran to C
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

}
