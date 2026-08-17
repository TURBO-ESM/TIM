#pragma once
/**
 * @file tim_stagger.hpp
 * @brief The staggering flags to specify where a field's values sit
 * within a horizontal grid cell.
 */

#include <AMReX_IntVect.H>

namespace TIM {

/// @brief The staggering (grid residency) of a field: the position of its
/// values within a horizontal grid cell. On an Arakawa C-grid, for instance,
/// Cell hosts the tracers, XFace the zonal velocity (u points), YFace the
/// meridional velocity (v points), and Node the vorticity and the Coriolis
/// parameter (q points).
enum class Stagger {
    Cell,   ///< Cell centers (h points).
    XFace,  ///< Cell faces normal to x (u points).
    YFace,  ///< Cell faces normal to y (v points).
    Node    ///< Cell corners (q points).
};

/// @brief The AMReX nodality of a staggering: 0 corresponds to cell centers
/// while 1 corresponds to points sitting on cell boundaries (nodes).
/// This class is used to specify horizontal staggering only
/// @param stagger The staggering.
/// @return The nodal flag vector (consumable by amrex::convert).
inline amrex::IntVect nodality(const Stagger stagger) {
    switch (stagger) {
        case Stagger::XFace: return amrex::IntVect(1, 0, 0);
        case Stagger::YFace: return amrex::IntVect(0, 1, 0);
        case Stagger::Node:  return amrex::IntVect(1, 1, 0);
        case Stagger::Cell:  break;
    }
    return amrex::IntVect(0, 0, 0);
}

}  // namespace TIM
