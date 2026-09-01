#pragma once
/**
 * @file tim_abort.hpp
 * @brief TIM's fatal-error path.
 */

#include <string_view>

namespace TIM {

/// @brief Abort the run with the diagnostic @p msg. Safe in any MPI/AMReX state
/// @param msg The diagnostic message to print before aborting.
[[noreturn]] void abort(std::string_view msg);

}  // namespace TIM
