#pragma once
/**
 * @file tim_runtime.hpp
 * @brief Startup/shutdown of the infrastructure runtime (MPI and AMReX).
 */

#include <mpi.h>

namespace TIM {

/// @brief Abort the run with the diagnostic @p msg. Safe in any MPI/AMReX state
/// @param msg The diagnostic message to print before aborting.
[[noreturn]] void abort(const char* msg);

/// @brief RAII encapsulation of the infrastructure runtime. Brings up MPI and
/// AMReX in order on construction, and shuts them down in reverse order on
/// destruction. (Future additional infrastructure dependencies, e.g. parallel I/O,
/// will be added here without affecting clients of this class.)
///
/// Whether TIM owns MPI is an explicit caller decision, made by picking a
/// constructor: (A) owner mode initializes and finalizes MPI (solo drivers);
/// (B) guest mode adopts a communicator from a caller that owns MPI (e.g. a
/// coupler) and never finalizes it. Either way AMReX runs as a guest on that
/// communicator.
///
/// Exactly one Runtime may exist per process lifetime (a copy would finalize
/// AMReX and MPI twice; AMReX's runtime state is global). Violations abort.
class Runtime {
public:
    /// @brief Owner mode: initialize MPI, then AMReX on MPI_COMM_WORLD.
    /// @param argc Command-line argument count.
    /// @param argv Command-line argument vector.
    ///
    /// Aborts if MPI is already initialized: Runtime owns the process in this
    /// mode, and a second would-be owner indicates a driver bug. Callers
    /// that receive their MPI from elsewhere use the guest-mode constructor.
    Runtime(int& argc, char**& argv);

    /// @brief Guest mode: adopt @p comm from a caller that owns MPI (e.g.,
    ///        a coupler) and initialize AMReX on it.
    /// @param comm The communicator TIM runs on. Must not be MPI_COMM_NULL
    ///             and must belong to an MPI that is initialized and not yet
    ///             finalized (aborts otherwise). In guest mode, MPI is not
    ///             freed or finalized by TIM.
    explicit Runtime(MPI_Comm comm);

    /// @brief Finalize AMReX. If in owner mode, finalize MPI too.
    ~Runtime();

    /// @brief The top-level communicator the infrastructure runs on.
    /// @return MPI_COMM_WORLD in owner mode. The adopted communicator in
    ///         guest mode.
    MPI_Comm comm() const { return comm_; }

    /// @brief Whether the infrastructure runtime is up. true if Runtime has been
    /// constructed and not yet destroyed (and hence MPI and AMReX are alive).
    /// @return true between a Runtime's construction and its destruction.
    static bool active();

    // Exactly one Runtime per process (see class comment).
    Runtime(const Runtime&) = delete;
    Runtime& operator=(const Runtime&) = delete;
    Runtime(Runtime&&) = delete;
    Runtime& operator=(Runtime&&) = delete;

private:
    /// @brief The MPI communicator TIM runs on.
    MPI_Comm comm_ = MPI_COMM_WORLD;
    /// @brief True in owner mode: the destructor also finalizes MPI.
    bool owns_mpi_ = false;
};

}  // namespace TIM
