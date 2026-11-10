// Smoke test for the AMReX-backed core (TIM::tim): AMReX comes up through the
// custom main (see test_tim_main.cpp), so the harness can build AMReX-backed
// fixtures the way protoMOMxx's tests/test_mom does. The PIO link check lives
// with the io layer (io/test_tim_io_smoke.cpp).

#include <gtest/gtest.h>

#include <AMReX.H>

TEST(Smoke, AmrexIsInitialized) {
    EXPECT_TRUE(amrex::Initialized());
}
