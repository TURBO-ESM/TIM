// Custom GoogleTest entry point for the TIM C++ unit tests
//
// Initialize/finalize AMReX.

#include <AMReX.H>
#include <gtest/gtest.h>

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    amrex::Initialize(argc, argv);
    int rc = RUN_ALL_TESTS();
    amrex::Finalize();
    return rc;
}
