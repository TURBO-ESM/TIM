/**
 * @file tim_abort.cpp
 * @brief Implementation of TIM's fatal-error path.hpp).
 */

#include <cstdio>
#include <cstdlib>

#include <mpi.h>

#include "tim_abort.hpp"

namespace TIM {

void abort(const char* msg) {
    std::fprintf(stderr, "%s\n", msg);
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
