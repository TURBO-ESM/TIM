#pragma once
/**
 * @file tim_coms_infra.hpp
 * @brief Checksum service of the TIM communication infrastructure.
 */

#include <cstdint>
#include <cstddef>

/// @brief TURBO Infrastructure for MOM (TIM) — the C++ infrastructure layer.
namespace TIM {

/// @brief Bitwise checksum of a distributed field.
/// @param field      Per-rank field data (host or device memory).
/// @param field_size Number of elements in @p field on this rank.
/// @param mask_val   Value marking masked elements (compared bitwise);
///                   pass nullptr for an unmasked checksum.
/// @return The global checksum (identical on every rank).
int64_t checksum(double* field, size_t field_size, double* mask_val);

}  // namespace TIM
