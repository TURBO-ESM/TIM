#include <AMReX_ParmParse.H>

#include "tim_profile.hpp"

namespace TIM {

/// Activates the AMReX profilling capability
void tim_set_profile(int level)  //!< Level of detail level = 1 basic profiling 
{
    amrex::ParmParse pp("amrex");
    pp.add("profile", level);
}

}
