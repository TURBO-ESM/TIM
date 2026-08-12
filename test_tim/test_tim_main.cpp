// Custom GoogleTest entry point for the TIM C++ unit tests
//
// Bring up the infrastructure layer in owner mode (TIM::Runtime initializes
// and finalizes MPI and AMReX), so every linked test exercises the real
// startup path and tests needing a live runtime (e.g. Domain) get one for
// free. The runtime tests keep their own mains: they test construction and
// teardown itself and inspect state after teardown, which a shared main
// cannot express.

#include <gtest/gtest.h>

#include "core/tim_runtime.hpp"

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    const TIM::Runtime runtime(argc, argv);
    return RUN_ALL_TESTS();
}
