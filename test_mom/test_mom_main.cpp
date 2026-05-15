// Custom GoogleTest entry point for the MOM C++ unit tests.
//
// Responsibilities:
//   * Strip GoogleTest's own flags via ::testing::InitGoogleTest.
//   * Parse --data-dir=PATH from the remaining argv. This flag is required.
//   * Initialize/Finalize AMReX.
//
// Tests retrieve the captured-data directory via test_mom::data_dir.

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <string_view>

#include <gtest/gtest.h>

#include <AMReX.H>

#include "data_dir.hpp"

namespace {

constexpr std::string_view kDataDirFlag = "--data-dir=";

// Pulls --data-dir=PATH out of argv. Returns the value, or an empty string if
// the flag is not present. Removes the matched argv slot so AMReX does not
// trip on an unknown option.
std::string extract_data_dir(int& argc, char** argv) {
    for (int i = 1; i < argc; ++i) {
        std::string_view a = argv[i];
        if (a.rfind(kDataDirFlag, 0) == 0) {
            std::string value(a.substr(kDataDirFlag.size()));
            for (int j = i; j < argc - 1; ++j) {
                argv[j] = argv[j + 1];
            }
            argv[argc - 1] = nullptr;
            --argc;
            return value;
        }
    }
    return {};
}

} // namespace

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);

    std::string path = extract_data_dir(argc, argv);

    // Build-time discovery (gtest_discover_tests POST_BUILD) invokes the
    // binary with --gtest_list_tests and no extra args. Honour that mode
    // without requiring --data-dir so test names can be enumerated.
    const bool list_only = ::testing::GTEST_FLAG(list_tests);
    if (!list_only) {
        if (path.empty()) {
            std::cerr
                << "ERROR: --data-dir=PATH is required.\n"
                   "       Point it at the directory containing the captured\n"
                   "       Fortran I/O .bin/.meta files (e.g. ppm_limit_pos.bin).\n";
            return EXIT_FAILURE;
        }
        if (!std::filesystem::is_directory(path)) {
            std::cerr << "ERROR: --data-dir='" << path
                      << "' is not a directory.\n";
            return EXIT_FAILURE;
        }
        test_mom::data_dir = path;
    }

    amrex::Initialize(argc, argv);
    int rc = RUN_ALL_TESTS();
    amrex::Finalize();
    return rc;
}
