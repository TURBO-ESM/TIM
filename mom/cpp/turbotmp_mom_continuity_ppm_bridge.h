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

struct OceanOBC;    // Undefined at the moment

#ifdef __cplusplus
extern "C" {
#endif

void turbotmp_ppm_limit_pos_bridge(const Box_C* bx_HOST, const RealArray_C* h_in_HOST, 
			RealArray_C* h_L_HOST, RealArray_C* h_R_HOST, const double h_min);
void turbotmp_ppm_limit_cw84_bridge(const Box_C* bx_HOST, const RealArray_C* h_in_HOST,
                        RealArray_C* h_L_HOST, RealArray_C* h_R_HOST);
void turbotmp_ppm_reconstruction_y_bridge(const Box_C* bx_HOST, const RealArray_C* h_in_HOST,
                        RealArray_C* h_S_HOST, RealArray_C* h_N_HOST, const RealArray_C* mask2dT_HOST,
                        const double h_min, const bool monotonic, const bool simple_2nd, OceanOBC* obc);

#ifdef __cplusplus
}
#endif
