#pragma once
/**
 * @file tim_coms_infra_C_API.h
 * @brief C API of the TIM checksum service (bind(C) surface for Fortran).
 */

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/// @brief C entry point for TIM::checksum (FMS mpp_chksum replacement).
/// @param field      Per-rank field data.
/// @param field_size Number of elements in @p field on this rank.
/// @param mask_val   Value marking masked elements (compared bitwise);
///                   pass NULL for an unmasked checksum.
/// @return The global checksum (identical on every rank).
int64_t tim_chksum_c(double* field, size_t field_size, double* mask_val);

#ifdef __cplusplus
}
#endif
