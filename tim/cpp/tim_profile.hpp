#pragma once

#include <AMReX_ParmParse.H>

namespace TIM {

#ifdef __cplusplus
extern "C" {
#endif
	void tim_set_profile(int level); //!< Activates the AMReX profilling capability
#ifdef __cplusplus
}
#endif

}
