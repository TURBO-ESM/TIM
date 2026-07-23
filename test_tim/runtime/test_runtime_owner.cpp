// Owner-mode test for the RAII encapsulation of the infrastructure runtime (TIM::Runtime).
//
// In owner mode Runtime owns the MPI communicator and hands it to AMReX,
// which then runs as a guest on that communicator. These tests verify that,
// with a Runtime alive, direct MPI calls on Runtime's communicator and AMReX
// data structures/collectives work side by side in one process.

#include <optional>

#include <gtest/gtest.h>

#include <AMReX.H>
#include <AMReX_BoxArray.H>
#include <AMReX_MultiFab.H>

#include "core/tim_runtime.hpp"

namespace {

// The one Runtime of this test binary, owned here and constructed/destroyed
// by main() via emplace()/reset() (exactly one Runtime may exist per process).
std::optional<TIM::Runtime> g_runtime;

// The communicator Runtime owns is usable for direct MPI calls alongside
// AMReX: a collective sum of one contribution per rank equals the comm size.
TEST(RuntimeOwner, MPICollectiveOnOwnedCommunicator) {
    ASSERT_TRUE(g_runtime.has_value());
    int size = -1;
    ASSERT_EQ(MPI_Comm_size(g_runtime->comm(), &size), MPI_SUCCESS);
    int one = 1;
    int sum = 0;
    ASSERT_EQ(MPI_Allreduce(&one, &sum, 1, MPI_INT, MPI_SUM, g_runtime->comm()), MPI_SUCCESS);
    EXPECT_EQ(sum, size);
}

// AMReX works as a guest: build a MultiFab and run an AMReX collective
// reduction (MultiFab::sum reduces over the communicator AMReX was handed).
TEST(RuntimeOwner, MultiFabCollectiveReduction) {
    ASSERT_TRUE(amrex::Initialized());
    const int n = 4;
    const amrex::Box domain(amrex::IntVect(0), amrex::IntVect(n - 1));
    const amrex::BoxArray ba(domain);
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab mf(ba, dm, 1, 0);
    mf.setVal(1.5);
    EXPECT_DOUBLE_EQ(mf.sum(), 1.5 * n * n * n);
}

}  // namespace

// Test-binary entry point: initialize GTest, bring up the infrastructure
// layer in owner mode (Runtime calls MPI_Init itself), and run all tests.
int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    g_runtime.emplace(argc, argv);
    const int rc = RUN_ALL_TESTS();
    g_runtime.reset();  // finalize AMReX and MPI before static teardown
    return rc;
}
