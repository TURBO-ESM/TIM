#pragma once
/**
 * @file turbotmp_bridge_c_types.h
 * @brief Plain C structs used to pass Fortran arrays and boxes across the
 *        temporary turbotmp bridge.
 */

#include <stdint.h>
#include <stddef.h>

/// @brief C view of a Fortran multidimensional real array passed across the bridge.
struct RealArray_C {
    double* data;   ///< Pointer to the multidimensional array data.
    int* shape;     ///< Array of per-dimension extents.
    int* lb;        ///< Per-dimension lower bounds.
    int* ub;        ///< Per-dimension upper bounds.
    int dim;        ///< Number of dimensions.
};

/// @brief C view of a MOM6 index box (start/end indices per dimension).
struct Box_C {
    int* idxS;      ///< Per-dimension start indices.
    int* idxE;      ///< Per-dimension end indices.
};

/* C has no implicit `struct` elision, so a C translation unit cannot spell these
 * as bare `RealArray_C` / `Box_C` the way C++ can.  The typedefs make the plain
 * names usable from both languages; they are redundant-but-legal in C++. */
typedef struct RealArray_C RealArray_C;
typedef struct Box_C       Box_C;
