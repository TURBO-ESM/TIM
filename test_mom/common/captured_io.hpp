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
//   one of: Box_t, RealArray_t, LogicalArray_t, real64, logical.
//
//   The .bin payload is big-endian and uses these layouts:
//
//     Box_t          : int32 ndim
//                      int32 idxS[ndim]
//                      int32 idxE[ndim]
//
//     RealArray_t    : int32 ndim
//                      int32 shape[ndim]
//                      (int32 lb, int32 ub)[ndim]    // interleaved
//                      int32 nelem
//                      float64 data[nelem]           // Fortran column-major
//
//     LogicalArray_t : same layout as RealArray_t, but the payload is
//                      int32 data[nelem] (0 = false, non-zero = true)
//                      instead of float64 -- mirrors the Fortran side's
//                      write_binaryLogical, which mirrors write_binary
//                      (RealArray_t) field-for-field.
//
//   A RealArray_t/LogicalArray_t entry can also be null-encoded: just
//   int32 ndim == -1, no shape/lb/ub/nelem/payload -- this is how the
//   Fortran side's write_binary(Logical) records an unassociated
//   (Fortran-optional-and-absent) array container. Fields the Fortran
//   shim captures unconditionally (regardless of %associated()) can
//   still be null-encoded this way; check is_associated() before calling
//   fab_host()/fab_device()/int_fab_host()/int_fab_device() on such a
//   field, since those throw on ndim == -1 (a real 2D/3D array is required).
//   A field can also be missing from the capture entirely -- when it is
//   itself a member of an optionally-associated derived type (e.g.
//   BT_cont%FA_u_W0), the whole capture block for it is skipped when the
//   parent container is unassociated, since rec%add() can't dereference an
//   unassociated parent. is_associated() treats a missing entry the same
//   as a null-encoded one: both mean "not associated," not a stale fixture.
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
//   the payload through a pinned-host buffer. int_fab_host()/int_fab_device()
//   are the LogicalArray_t counterparts, returning amrex::IArrayBox.

// SKILLS: 0.3.1

#pragma once

#include <cstddef>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>

namespace test_mom {

// Loads a .meta / .bin pair given the base path (no extension). The .meta
// directory is parsed eagerly; payload values are decoded on each accessor
// call. fab_host()/fab_device()/int_fab_host()/int_fab_device() return a
// fresh FArrayBox/IArrayBox by move, so the caller owns it and may mutate it
// (e.g. pass to a kernel as in/out storage).
class CapturedFile {
public:
    explicit CapturedFile(const std::filesystem::path& base);

    amrex::Box       box           (const std::string& name) const;
    amrex::FArrayBox fab_host      (const std::string& name) const;
    amrex::FArrayBox fab_device    (const std::string& name) const;
    amrex::IArrayBox int_fab_host  (const std::string& name) const;
    amrex::IArrayBox int_fab_device(const std::string& name) const;
    double           real64       (const std::string& name) const;
    bool             logical      (const std::string& name) const;
    int              integer      (const std::string& name) const;
    // True unless a RealArray_t/LogicalArray_t entry is null-encoded
    // (ndim == -1, i.e. the Fortran container was unassociated at capture
    // time). Check this before calling fab_host()/fab_device()/
    // int_fab_host()/int_fab_device() on a field the Fortran shim captures
    // unconditionally.
    bool             is_associated(const std::string& name) const;

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
