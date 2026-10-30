// SKILLS: 0.3.2
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

/// @brief C view of a Fortran multidimensional integer array passed across the bridge.
struct IntArray_C {
    int* data;      ///< Pointer to the multidimensional array data.
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

/// @brief A C-interoperable mirror of grid_core_type, with its members in the same order.
struct grid_core_C {
    struct RealArray_C mask2dT;     ///< 0 for land points and 1 for ocean points on the h-grid
                                    ///  [nondim].
    struct RealArray_C areaT;       ///< The area of an h-cell [L2 ~> m2].
    struct RealArray_C IareaT;      ///< 1/areaT [L-2 ~> m-2].
    struct RealArray_C dxT;         ///< dxT is delta x at h points [L ~> m].
    struct RealArray_C dyT;         ///< dyT is delta y at h points [L ~> m].
    struct RealArray_C IdxT;        ///< 1/dxT [L-1 ~> m-1].
    struct RealArray_C IdyT;        ///< IdyT is 1/dyT [L-1 ~> m-1].
    struct RealArray_C bathyT;      ///< Ocean bottom depth, referenced to Z_ref at tracer points.
                                    ///  bathyT is in depth units and positive *below* Z_ref
                                    ///  [Z ~> m].
    struct RealArray_C mask2dCu;    ///< 0 for boundary points and 1 for ocean points on the u grid
                                    ///  [nondim].
    struct RealArray_C dy_Cu;       ///< The unblocked lengths of the u-faces of the h-cell [L ~> m].
    struct RealArray_C IdxCu;       ///< 1/dxCu [L-1 ~> m-1].
    struct RealArray_C dxCu;        ///< dxCu is delta x at u points [L ~> m].
    struct RealArray_C mask2dCv;    ///< 0 for boundary points and 1 for ocean points on the v grid
                                    ///  [nondim].
    struct RealArray_C dx_Cv;       ///< The unblocked lengths of the v-faces of the h-cell [L ~> m].
    struct RealArray_C IdyCv;       ///< 1/dyCv [L-1 ~> m-1].
    struct RealArray_C dyCv;        ///< dyCv is delta y at v points [L ~> m].
    struct RealArray_C CoriolisBu;  ///< The Coriolis parameter at corner points [T-1 ~> s-1].
    struct RealArray_C meanSL;      ///< Spatially varying time mean sea level, referenced to Z_ref
                                    ///  at tracer points. meanSL is in height units and positive
                                    ///  *above* Z_ref. It is used a) as the height where p = p_atm
                                    ///  or zero; b) to calculate time mean thickness of the water
                                    ///  column, where mean thickness = max(meanSL + bathyT, 0.0).
                                    ///  meanSL is 2D for the consideration of a domain with
                                    ///  spatically varying mean height, e.g. the Great Lakes system
                                    ///  [Z ~> m].
    struct RealArray_C geoLatT;     ///< The geographic latitude at tracer (h) points
                                    ///  [degrees_N] or [km] or [m]
    struct RealArray_C IdyCu;       ///< 1/dyCu [L-1 ~> m-1].
    struct RealArray_C areaCu;      ///< The areas of the u-grid cells [L2 ~> m2].
    struct RealArray_C IdxCv;       ///< 1/dxCv [L-1 ~> m-1].
    struct RealArray_C areaCv;      ///< The areas of the v-grid cells [L2 ~> m2].
    struct RealArray_C mask2dBu;    ///< 0 for boundary points and 1 for ocean points on the q grid
                                    ///  [nondim].
    struct RealArray_C dxBu;        ///< dxBu is delta x at q points [L ~> m].
    struct RealArray_C dyBu;        ///< dyBu is delta y at q points [L ~> m].
    struct RealArray_C IareaBu;     ///< IareaBu = 1/areaBu [L-2 ~> m-2].
    struct RealArray_C geoLonT;     ///< The geographic longitude at tracer (h) points
                                    ///  [degrees_E] or [km] or [m]
    struct RealArray_C df_dx;       ///< Derivative d/dx f (Coriolis parameter) at h-points
                                    ///  [T-1 L-1 ~> s-1 m-1].
    struct RealArray_C df_dy;       ///< Derivative d/dy f (Coriolis parameter) at h-points
                                    ///  [T-1 L-1 ~> s-1 m-1].
    struct RealArray_C IareaCu;     ///< The masked inverse areas of u-grid cells [L-2 ~> m-2].
    struct RealArray_C IareaCv;     ///< The masked inverse areas of v-grid cells [L-2 ~> m-2].
    struct RealArray_C IdxBu;       ///< 1/dxBu [L-1 ~> m-1].
    struct RealArray_C IdyBu;       ///< 1/dyBu [L-1 ~> m-1].
    struct RealArray_C Coriolis2Bu; ///< The square of the Coriolis parameter at corner points
                                    ///  [T-2 ~> s-2].
};

/// @brief A C-interoperable mirror of grid_OBC_type, with its members in the same order.
struct grid_OBC_C {
    struct RealArray_C OBCmaskCu;     ///< 0 for boundary or OBC points and 1 for ocean points on the u
                                      ///  grid [nondim].
    struct RealArray_C IdxCu_OBCmask; ///< 1/dxCu or 0 at boundary or OBC points [L-1 ~> m-1].
    struct RealArray_C OBCmaskCv;     ///< 0 for boundary or OBC points and 1 for ocean points on the v
                                      ///  grid [nondim].
    struct RealArray_C IdyCv_OBCmask; ///< 1/dxCv or 0 at boundary or OBC points [L-1 ~> m-1].
};

/* C has no implicit `struct` elision, so a C translation unit cannot spell these
 * as bare `RealArray_C` / `Box_C` the way C++ can.  The typedefs make the plain
 * names usable from both languages; they are redundant-but-legal in C++. */

/// @brief Alias letting C callers spell the array view without the `struct` keyword.
typedef struct RealArray_C    RealArray_C;
/// @brief Alias letting C callers spell the integer array view without the `struct` keyword.
typedef struct IntArray_C     IntArray_C;
/// @brief Alias letting C callers spell the logical array view without the `struct` keyword.
typedef struct LogicalArray_C LogicalArray_C;
/// @brief Alias letting C callers spell the index box without the `struct` keyword.
typedef struct Box_C          Box_C;
/// @brief Alias letting C callers spell the grid core mirror without the `struct` keyword.
typedef struct grid_core_C    grid_core_C;
/// @brief Alias letting C callers spell the grid OBC mirror without the `struct` keyword.
typedef struct grid_OBC_C     grid_OBC_C;
