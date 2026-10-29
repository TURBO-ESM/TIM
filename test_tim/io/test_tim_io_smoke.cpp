// Smoke test for the AMReX-free I/O layer (TIM::tim_io)

// Build-enforced AMReX isolation: this TU links tim_io only, so AMReX headers
// can become visible here only through a usage requirement leaking from tim_io
// or tim_base. (Assumes AMReX lives under its own prefix, as in every TIM build,
// not a system include directory.)
#if __has_include(<AMReX.H>)
#error "AMReX headers are visible to a TIM::tim_io consumer; tim_io must not link AMReX"
#endif

#include <gtest/gtest.h>
#include <mpi.h>
#include <pio.h>

TEST(IoSmoke, PioSystemInitializesAndFrees) {
    int iosysid = -1;
    ASSERT_EQ(PIOc_Init_Intracomm(MPI_COMM_WORLD, /*num_iotasks=*/1, /*stride=*/1, /*base=*/0,
                                  PIO_REARR_BOX, &iosysid),
              PIO_NOERR);
    EXPECT_EQ(PIOc_free_iosystem(iosysid), PIO_NOERR);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);
    const int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
