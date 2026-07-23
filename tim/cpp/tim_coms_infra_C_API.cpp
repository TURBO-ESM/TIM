/**
 * @file tim_coms_infra_C_API.cpp
 * @brief C-API translation unit of the TIM checksum service.
 */

#include "tim_coms_infra_C_API.h"
#include "tim_coms_infra.hpp"

int64_t tim_chksum_c(double* field, size_t field_size, double* mask_val) {
    return TIM::checksum(field, field_size, mask_val);
}
