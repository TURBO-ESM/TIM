# TURBO Infrastructure for MOM6

TIM (TURBO Infrastructure for MOM6) is a work-in-progress effort to build an
AMReX-based software infrastructure for the Modular Ocean Model version 6 (MOM6),
with the goal of enabling GPU acceleration, performance portability, and
ultra-high-resolution ocean simulations on current and emerging heterogeneous
architectures.

TIM is being developed incrementally. The initial phase extracts and refactors the
subset of FMS functionality required by MOM6, re-implementing it in modern C++ with
AMReX to provide infrastructure features including domain decomposition, MPI/GPU
parallelism, device-portable kernels, memory management, I/O, and diagnostics.

These pages document the C++ infrastructure layer under `tim/cpp` (namespace `TIM`).
Browse the [class list](annotated.html) to get started.
