/**
 * @file tim_profile.cpp
 * @brief Implementation of the AMReX profiling activation.
 */

#include <AMReX_ParmParse.H>

#include "tim_profile.hpp"

namespace TIM {

void tim_set_profile(int level)
{
    amrex::ParmParse pp("amrex");
    pp.add("profile", level);
}

}  // namespace TIM
