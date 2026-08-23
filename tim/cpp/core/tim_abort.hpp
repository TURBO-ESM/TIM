#pragma once
/**
 * @file tim_abort.hpp
 * @brief TIM's fatal-error path.
 */

namespace TIM {

/// @brief Abort the run with the diagnostic @p msg. Safe in any MPI/AMReX state
/// @param msg The diagnostic message to print before aborting.
[[noreturn]] void abort(const char* msg);

}  // namespace TIM
