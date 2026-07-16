#pragma once
/**
 * @file tim_profile.hpp
 * @brief Runtime activation of the AMReX profiling capability.
 */

namespace TIM {

#ifdef __cplusplus
extern "C" {
#endif

/// @brief Activates the AMReX profiling capability.
/// @param level Profiling detail level (1 = basic profiling).
void tim_set_profile(int level);

#ifdef __cplusplus
}
#endif

}  // namespace TIM
