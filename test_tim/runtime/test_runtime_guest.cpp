// Guest-mode test for the RAII encapsulation of the infrastructure runtime
// (TIM::Runtime).
//
// In guest mode the caller (here: this binary's main, representing a coupler
// or a cap) owns MPI and hands Runtime the communicator to run on. These tests
// verify that Runtime adopts the communicator and brings up AMReX on it, and
// after Runtime is destroyed inside main(), that MPI is still alive. 

#include <cstdio>
#include <cstdlib>
#include <optional>

#include <gtest/gtest.h>

#include <AMReX.H>
#include <AMReX_BoxArray.H>
#include <AMReX_MultiFab.H>

#include "core/tim_runtime.hpp"

namespace {

[[noreturn]] void fail(const char* msg) {
    std::fprintf(stderr, "test_runtime_guest: FAIL -- %s\n", msg);
    std::exit(1);
}

// The one Runtime of this test binary, owned here and constructed/destroyed
// by main() via emplace()/reset() (exactly one Runtime may exist per process).
std::optional<TIM::Runtime> g_runtime;

// Check that Runtime adopts exactly the communicator it was handed, and it is
// usable for direct MPI calls alongside AMReX.
TEST(RuntimeGuest, MPICollectiveOnAdoptedCommunicator) {
    ASSERT_TRUE(g_runtime.has_value());
    EXPECT_EQ(g_runtime->comm(), MPI_COMM_WORLD);
    int size = -1;
    ASSERT_EQ(MPI_Comm_size(g_runtime->comm(), &size), MPI_SUCCESS);
    int one = 1;
    int sum = 0;
    ASSERT_EQ(MPI_Allreduce(&one, &sum, 1, MPI_INT, MPI_SUM, g_runtime->comm()), MPI_SUCCESS);
    EXPECT_EQ(sum, size);
}

// Check that AMReX was brought up as a guest on the adopted communicator.
TEST(RuntimeGuest, MultiFabCollectiveReduction) {
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

// Test-binary entry point: own MPI like a coupler would, run the tests with
// a guest-mode Runtime alive, then verify the guest left MPI to its owner.
int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
        fail("MPI_Init failed");
    }
    g_runtime.emplace(MPI_COMM_WORLD);
    const int rc = RUN_ALL_TESTS();
    g_runtime.reset();
    // Runtime is destroyed; the guest must not have finalized MPI, and the
    // caller's MPI must still be usable.
    int finalized = 1;
    MPI_Finalized(&finalized);
    if (finalized) {
        fail("guest-mode Runtime finalized MPI");
    }
    int size = -1;
    if (MPI_Comm_size(MPI_COMM_WORLD, &size) != MPI_SUCCESS || size < 1) {
        fail("MPI unusable after Runtime teardown");
    }
    MPI_Finalize();
    return rc;
}
