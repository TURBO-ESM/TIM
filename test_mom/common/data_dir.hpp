// Directory holding the captured Fortran I/O fixtures (.bin/.meta files).
// Set once at startup by the test main; tests read it to build per-fixture
// paths (e.g. data_dir / "ppm_limit_pos").
#pragma once

#include <filesystem>

namespace test_mom {

inline std::filesystem::path data_dir;

} // namespace test_mom
