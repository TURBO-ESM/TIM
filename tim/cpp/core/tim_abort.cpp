/**
 * @file tim_abort.cpp
 * @brief Implementation of TIM's fatal-error path.
 */

#include "tim_abort.hpp"

#include <cstdio>
#include <cstdlib>

#include <mpi.h>

namespace TIM {

void abort(const std::string_view msg) {
    // Not necessarily null-terminated, so bound the write by its size.
    std::fprintf(stderr, "%.*s\n", static_cast<int>(msg.size()), msg.data());
    std::fflush(stderr);
    int mpi_initialized = 0;
    MPI_Initialized(&mpi_initialized);
    if (mpi_initialized) {
        int mpi_finalized = 0;
        MPI_Finalized(&mpi_finalized);
        if (!mpi_finalized) {
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    std::abort();
}

}  // namespace TIM
