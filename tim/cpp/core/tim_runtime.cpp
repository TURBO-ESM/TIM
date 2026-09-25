/**
 * @file tim_runtime.cpp
 * @brief Implementation of TIM::Runtime, the RAII encapsulation of the
 * infrastructure runtime.
 */

#include <atomic>

#include <AMReX.H>

#include "tim_runtime.hpp"
#include "tim_abort.hpp"

namespace TIM {

namespace {

// Enforce that exactly one Runtime exists per process lifetime. If a second is
// constructed, abort with a message. (Deliberately never reset: the
// one-per-lifetime rule outlives the Runtime itself.)
std::atomic<bool> runtime_exists{false};

// Whether the runtime is up right now: set once a constructor completes,
// cleared when destruction begins. Backs Runtime::active().
std::atomic<bool> runtime_active{false};

void enforceSingleRuntime() {
    if (runtime_exists.exchange(true)) {
        TIM::abort("TIM::Runtime: a Runtime already exists in this process; "
                   "exactly one is allowed per process lifetime.");
    }
}

}  // namespace

bool Runtime::active() {
    return runtime_active.load();
}

Runtime::Runtime(int& argc, char**& argv) : owns_mpi_(true) {
    enforceSingleRuntime();

    // Check if MPI is already finalized or initialized. If so, abort with a message.
    int mpi_already_finalized = 0;
    MPI_Finalized(&mpi_already_finalized);
    if (mpi_already_finalized) {
        TIM::abort("TIM::Runtime (owner mode): MPI has already been "
                   "finalized and cannot be initialized again.");
    }
    int mpi_already_initialized = 0;
    MPI_Initialized(&mpi_already_initialized);
    if (mpi_already_initialized) {
        TIM::abort("TIM::Runtime (owner mode): MPI is already initialized, but this "
                   "constructor owns MPI startup/shutdown. Construct Runtime before "
                   "any other MPI use, or hand the existing communicator to the "
                   "guest-mode constructor.");
    }
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
        TIM::abort("TIM::Runtime (owner mode): MPI_Init failed.");
    }

    // Having initialized MPI, initialize AMReX.
    amrex::Initialize(comm_);

    runtime_active.store(true);
}

Runtime::Runtime(MPI_Comm comm) : comm_(comm), owns_mpi_(false) {
    enforceSingleRuntime();

    // Check that the adopted communicator is valid and that MPI is initialized but not finalized.
    if (comm_ == MPI_COMM_NULL) {
        TIM::abort("TIM::Runtime (guest mode): the adopted communicator is "
                   "MPI_COMM_NULL; pass a valid communicator.");
    }
    int mpi_already_initialized = 0;
    MPI_Initialized(&mpi_already_initialized);
    if (!mpi_already_initialized) {
        TIM::abort("TIM::Runtime (guest mode): MPI is not initialized. This "
                   "constructor adopts a communicator from a caller that owns "
                   "MPI; use the owner-mode constructor if TIM should own MPI "
                   "startup/shutdown.");
    }
    int mpi_already_finalized = 0;
    MPI_Finalized(&mpi_already_finalized);
    if (mpi_already_finalized) {
        TIM::abort("TIM::Runtime (guest mode): MPI has already been "
                   "finalized; construct Runtime while MPI is alive.");
    }

    // Having verified that MPI is alive and the communicator is valid, initialize AMReX on it.
    amrex::Initialize(comm_);

    runtime_active.store(true);
}

Runtime::~Runtime() {
    runtime_active.store(false);
    amrex::Finalize();
    if (owns_mpi_) {
        MPI_Finalize();
    }
}

}  // namespace TIM
