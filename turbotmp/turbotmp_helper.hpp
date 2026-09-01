#pragma once
#include <AMReX_Array4.H>
#include <AMReX_MultiFab.H>

namespace turbotmp {

struct A4Box
{
    amrex::Box bx;
    amrex::Real* data = nullptr;     // Raw C-pointer on device C order
    amrex::Array4<amrex::Real> arr;  //  AMReX view C order

    amrex::Real* data_f = nullptr;  // Raw C-pointer on device Fortran order

    int nx, ny, nz, ncomp;
};

   // lbx/lby/lbz are the array's Fortran index-space lower bounds (i.e., the
   // values returned by lbound() on the Fortran side). For a faked-out 2D
   // dimension, pass 1 (matching how nz=1 is passed for 2D arrays).
   // The returned A4Box's box -- and therefore a4.arr's valid index range --
   // is positioned at the array's true absolute location (lb-1 .. lb-1+n-1),
   // not just [0, n-1]. Staggered-grid arrays (u/v-point fields) legitimately
   // have lb != 1, and shared iteration boxes require each array4 to be
   // positioned at its own real offset so indexing hits the correct element.
   A4Box make_array4(int nx, int ny, int nz, int ncomp, int lbx, int lby, int lbz);
   void free_array4(A4Box& a4);
   void copy_FortranHost_to_array4(const double* f, A4Box& a4);
   void copy_array4_to_FortranHost(const A4Box& a4, double* f);

// IntA4Box mirrors A4Box for integer-typed data. Serves both LogicalArray_C
// (integer-encoded 0/1) and IntArray_C -- both are plain `int*` on the C side.
struct IntA4Box
{
    amrex::Box bx;
    int* data = nullptr;     // Raw C-pointer on device C order
    amrex::Array4<int> arr;  //  AMReX view C order

    int* data_f = nullptr;  // Raw C-pointer on device Fortran order

    int nx, ny, nz, ncomp;
};

   // See make_array4()'s lbx/lby/lbz comment -- identical convention.
   IntA4Box make_int_array4(int nx, int ny, int nz, int ncomp, int lbx, int lby, int lbz);
   void free_int_array4(IntA4Box& a4);
   void copy_FortranHost_to_int_array4(const int* f, IntA4Box& a4);
   void copy_int_array4_to_FortranHost(const IntA4Box& a4, int* f);
}
