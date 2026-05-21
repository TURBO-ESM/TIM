// Reader for the captured Fortran I/O fixtures (.meta + big-endian .bin).
//
// Format
// ------
//   <basename>.meta is a small text file. Each line has three whitespace-
//   separated tokens:
//
//       <name>  <type>  <offset>
//
//   where <offset> is a 1-based byte offset into <basename>.bin. <type> is
//   one of: Box_t, RealArray_t, real64, logical.
//
//   The .bin payload is big-endian and uses these layouts:
//
//     Box_t       : int32 ndim
//                   int32 idxS[ndim]
//                   int32 idxE[ndim]
//
//     RealArray_t : int32 ndim
//                   int32 shape[ndim]
//                   (int32 lb, int32 ub)[ndim]    // interleaved
//                   int32 nelem
//                   float64 data[nelem]           // Fortran column-major
//
//     real64      : float64
//     logical     : int32 (0 = false, non-zero = true)
//
// Index conventions
// -----------------
//   * idxS / idxE / lb / ub are Fortran 1-based. box() subtracts 1 to convert
//     into AMReX's 0-based world (matching the turbotmp bridge).
//   * Captured arrays in this project happen to have lb == 1 for every
//     dimension. fab() asserts that and builds a storage box that spans
//     [0, shape-1] in each dim.
//
// AMReX dependency
// ----------------
//   fab_host() allocates on The_Pinned_Arena(); fab_device() allocates on
//   The_Arena() (device memory on GPU builds, host on CPU builds) and stages
//   the payload through a pinned-host buffer.

#pragma once

#include <cstddef>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>

namespace test_mom {

// Loads a .meta / .bin pair given the base path (no extension). The .meta
// directory is parsed eagerly; payload values are decoded on each accessor
// call. fab_host()/fab_device() return a fresh FArrayBox by move, so the
// caller owns it and may mutate it (e.g. pass to a kernel as in/out storage).
class CapturedFile {
public:
    explicit CapturedFile(const std::filesystem::path& base);

    amrex::Box       box       (const std::string& name) const;
    amrex::FArrayBox fab_host  (const std::string& name) const;
    amrex::FArrayBox fab_device(const std::string& name) const;
    double           real64    (const std::string& name) const;
    bool             logical   (const std::string& name) const;

private:
    struct Entry {
        std::string type;
        std::size_t offset; // 1-based into bin_
    };

    const Entry& lookup(const std::string& name, const char* want_type) const;

    std::unordered_map<std::string, Entry> entries_;
    std::vector<std::byte>                 bin_;
    std::filesystem::path                  base_;
};

} // namespace test_mom
