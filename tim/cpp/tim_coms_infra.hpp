#pragma once

#include <AMReX_Array4.H>
#include <AMReX_Box.H>

namespace TIM {
    amrex::Long checksum(amrex::Box const& bx, amrex::Array4<amrex::Real> const& arr);
    amrex::Long checksum(amrex::Box const& bx, amrex::Array4<amrex::Real> const& arr, amrex::Real mask);
}
