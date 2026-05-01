#pragma once

#include <stdint.h>
#include <stddef.h>

struct RealArray_C {
    double* data;   // Pointer to multidimensional array
    int* shape;     // An array of dimension extents
    int* lb;        // Lower bounds
    int* ub;        // Upper bounds
    int dim;        // The number of dimension
};
struct Box_C {
    int* idxS;
    int* idxE;
};

#ifdef __cplusplus
extern "C" {
#endif

void turbotmp_thickness_to_dz_3d_bridge(const Box_C* bx,
                                        const RealArray_C* h,
                                        RealArray_C* dz,
                                        const RealArray_C* spv_avg,
                                        const bool boussinesq,
                                        const double h_to_z,
                                        const double h_to_rz,
                                        const bool has_spv);

#ifdef __cplusplus
}
#endif
