#pragma once
/**
 * @file tim_stagger.hpp
 * @brief The staggering flags to specify where a field's values sit
 * within a horizontal grid cell.
 */

namespace TIM {

/// @brief The staggering (grid residency) of a field: the position of its
/// values within a horizontal grid cell. On an Arakawa C-grid, for instance,
/// Cell hosts the tracers, XFace the zonal velocity (u points), YFace the
/// meridional velocity (v points), and Node the vorticity and the Coriolis
/// parameter (q points).
enum class Stagger {
    Cell,   ///< Cell centers (h points).
    XFace,  ///< Cell faces normal to x (u points; the west face).
    YFace,  ///< Cell faces normal to y (v points; the south face).
    Node    ///< Cell corners (q points; the southwest corner).
};

}  // namespace TIM
