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

/// @brief C view of a Fortran multidimensional logical array passed across the bridge.
/// @details The Fortran side stores a `logical` array natively but exposes it here through
/// an integer-encoded (0/1) shadow buffer (`LogicalArray_t%data_c`, filled by `%to_c()` and
/// read back by `%from_c()`), since Fortran `logical` is not itself C-interoperable. `data`
/// therefore always points at plain `int`s, never at genuine Fortran logicals.
struct LogicalArray_C {
    int* data;      ///< Pointer to the integer-encoded (0/1) array data.
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
