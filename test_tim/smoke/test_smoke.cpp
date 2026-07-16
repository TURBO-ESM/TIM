// PIO linking test for the TIM C++ I/O and diagnostics implementation 
//
//   * Initialize AMReX via the custom main (see test_tim_main.cpp), so the
//     harness can build AMReX-backed fixtures the way test_mom does.
//   * Confirm that the ParallelIO C library links and its header resolves, so the
//     io/ layer can #include <pio.h> and call libpioc.

#include <gtest/gtest.h>

#include <AMReX.H>

#include <pio.h>

TEST(Smoke, AmrexIsInitialized) {
    EXPECT_TRUE(amrex::Initialized());
}

TEST(Smoke, ParallelIoLinks) {
    EXPECT_EQ(static_cast<int>(PIO_IOTYPE_NETCDF), 2);
    void* pio_entry = reinterpret_cast<void*>(&PIOc_Init_Intracomm);
    EXPECT_NE(pio_entry, nullptr);
}
