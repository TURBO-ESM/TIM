#pragma once

#include "turbotmp_bridge_c_types.h"

#ifdef __cplusplus
extern "C" {
#endif

void turbotmp_thickness_to_dz_3d_bridge(const Box_C* bx_HOST,
                                        const RealArray_C* h_HOST,
                                        RealArray_C* dz_HOST,
                                        const RealArray_C* spv_avg_HOST,
                                        const bool boussinesq,
                                        const double h_to_z,
                                        const double h_to_rz,
                                        const bool has_spv);

#ifdef __cplusplus
}
#endif
