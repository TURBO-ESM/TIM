/**
 * @file tim_runtime.cpp
 * @brief Implementation of TIM::Runtime, the RAII encapsulation of the
 * infrastructure runtime.
 */

#include <atomic>

#include <AMReX.H>

#include "tim_runtime.hpp"

namespace TIM {

namespace {

// Enforce that exactly one Runtime exists per process lifetime. If a second is
// constructed, abort with a message.
std::atomic<bool> runtime_exists{false};

void enforceSingleRuntime() {
    if (runtime_exists.exchange(true)) {
        amrex::Abort("TIM::Runtime: a Runtime already exists in this process; "
                     "exactly one is allowed per process lifetime.");
    }
}

}  // namespace

Runtime::Runtime(int& argc, char**& argv) : owns_mpi_(true) {
    enforceSingleRuntime();
    int mpi_already_initialized = 0;
    MPI_Initialized(&mpi_already_initialized);
    if (mpi_already_initialized) {
        amrex::Abort("TIM::Runtime (owner mode): MPI is already initialized, but this "
                     "constructor owns MPI startup/shutdown. Construct Runtime before "
                     "any other MPI use, or hand the existing communicator to the "
                     "guest-mode constructor.");
    }
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
        amrex::Abort("TIM::Runtime (owner mode): MPI_Init failed.");
    }
    // AMReX is a guest on the communicator Runtime owns: it still sets up its
    // internals (arenas, devices, ...) but does not own argc/argv or MPI.
    amrex::Initialize(comm_);
}

Runtime::Runtime(MPI_Comm comm) : comm_(comm), owns_mpi_(false) {
    enforceSingleRuntime();
    int mpi_already_initialized = 0;
    MPI_Initialized(&mpi_already_initialized);
    if (!mpi_already_initialized) {
        amrex::Abort("TIM::Runtime (guest mode): MPI is not initialized. This "
                     "constructor adopts a communicator from a caller that owns "
                     "MPI; use the owner-mode constructor if TIM should own MPI "
                     "startup/shutdown.");
    }
    amrex::Initialize(comm_);
}

Runtime::~Runtime() {
    amrex::Finalize();
    if (owns_mpi_) {
        MPI_Finalize();
    }
}

}  // namespace TIM
