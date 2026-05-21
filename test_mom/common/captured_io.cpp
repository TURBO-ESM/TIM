#include "captured_io.hpp"

#include <cstdint>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <AMReX_Arena.H>
#include <AMReX_Gpu.H>

namespace test_mom {

namespace {

// --------------------------------------------------------------------------
// Big-endian primitive readers from an in-memory buffer
// --------------------------------------------------------------------------

uint32_t read_be_u32(const std::vector<std::byte>& buf, std::size_t& off) {
    if (off + 4 > buf.size()) {
        throw std::runtime_error("captured_io: truncated bin (u32 at " +
                                 std::to_string(off) + ")");
    }
    auto p = reinterpret_cast<const unsigned char*>(buf.data()) + off;
    off += 4;
    return (uint32_t(p[0]) << 24) | (uint32_t(p[1]) << 16) |
           (uint32_t(p[2]) << 8)  | (uint32_t(p[3]));
}

int32_t read_be_i32(const std::vector<std::byte>& buf, std::size_t& off) {
    return static_cast<int32_t>(read_be_u32(buf, off));
}

double read_be_f64(const std::vector<std::byte>& buf, std::size_t& off) {
    if (off + 8 > buf.size()) {
        throw std::runtime_error("captured_io: truncated bin (f64 at " +
                                 std::to_string(off) + ")");
    }
    auto p = reinterpret_cast<const unsigned char*>(buf.data()) + off;
    unsigned char swapped[8];
    for (int i = 0; i < 8; ++i) swapped[i] = p[7 - i];
    double v;
    std::memcpy(&v, swapped, 8);
    off += 8;
    return v;
}

std::vector<std::byte> slurp(const std::filesystem::path& p) {
    std::ifstream in(p, std::ios::binary | std::ios::ate);
    if (!in) {
        throw std::runtime_error("captured_io: cannot open " + p.string());
    }
    auto sz = in.tellg();
    in.seekg(0, std::ios::beg);
    std::vector<std::byte> buf(static_cast<std::size_t>(sz));
    in.read(reinterpret_cast<char*>(buf.data()), sz);
    return buf;
}

} // namespace

// --------------------------------------------------------------------------
// CapturedFile
// --------------------------------------------------------------------------

CapturedFile::CapturedFile(const std::filesystem::path& base) : base_(base) {
    std::ifstream meta(base.string() + ".meta");
    if (!meta) {
        throw std::runtime_error("captured_io: cannot open " +
                                 base.string() + ".meta");
    }
    std::string line;
    while (std::getline(meta, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string name, type;
        std::size_t offset;
        if (!(ss >> name >> type >> offset)) {
            throw std::runtime_error("captured_io: bad meta line: '" + line +
                                     "' in " + base.string() + ".meta");
        }
        entries_.emplace(std::move(name), Entry{std::move(type), offset});
    }
    bin_ = slurp(base.string() + ".bin");
}

const CapturedFile::Entry&
CapturedFile::lookup(const std::string& name, const char* want_type) const {
    auto it = entries_.find(name);
    if (it == entries_.end()) {
        throw std::runtime_error("captured_io: missing entry '" + name +
                                 "' in " + base_.string());
    }
    if (it->second.type != want_type) {
        throw std::runtime_error("captured_io: entry '" + name +
                                 "' has type '" + it->second.type +
                                 "', expected '" + want_type + "'");
    }
    return it->second;
}

amrex::Box CapturedFile::box(const std::string& name) const {
    const auto& e = lookup(name, "Box_t");
    std::size_t off = e.offset - 1;
    int ndim = read_be_i32(bin_, off);
    if (ndim < 1 || ndim > 3) {
        throw std::runtime_error("captured_io: Box_t ndim=" +
                                 std::to_string(ndim) + " for '" + name + "'");
    }
    // Pad unused dims with 1..1 so the resulting AMReX Box is always 3D.
    int idxS[3] = {1, 1, 1};
    int idxE[3] = {1, 1, 1};
    for (int i = 0; i < ndim; ++i) idxS[i] = read_be_i32(bin_, off);
    for (int i = 0; i < ndim; ++i) idxE[i] = read_be_i32(bin_, off);
    return amrex::Box(amrex::IntVect(idxS[0] - 1, idxS[1] - 1, idxS[2] - 1),
                      amrex::IntVect(idxE[0] - 1, idxE[1] - 1, idxE[2] - 1));
}

amrex::FArrayBox CapturedFile::fab_host(const std::string& name) const {
    const auto& e = lookup(name, "RealArray_t");
    std::size_t off = e.offset - 1;
    int ndim = read_be_i32(bin_, off);
    if (ndim != 2 && ndim != 3) {
        throw std::runtime_error(
            "captured_io: fab expects 2D or 3D arrays (got ndim=" +
            std::to_string(ndim) + ") for '" + name + "'");
    }
    int shape[3] = {1, 1, 1};
    for (int i = 0; i < ndim; ++i) shape[i] = read_be_i32(bin_, off);
    for (int i = 0; i < ndim; ++i) {
        int lb = read_be_i32(bin_, off);
        int ub = read_be_i32(bin_, off);
        if (lb != 1) {
            throw std::runtime_error(
                "captured_io: fab assumes Fortran lb==1 (got " +
                std::to_string(lb) + " on dim " + std::to_string(i) +
                ") for '" + name + "'");
        }
        if (ub - lb + 1 != shape[i]) {
            throw std::runtime_error(
                "captured_io: fab bounds inconsistent with shape on dim " +
                std::to_string(i) + " (lb=" + std::to_string(lb) + ", ub=" +
                std::to_string(ub) + ", shape=" + std::to_string(shape[i]) +
                ") for '" + name + "'");
        }
    }
    int32_t nelem = read_be_i32(bin_, off);
    int64_t expected = int64_t{shape[0]} * shape[1] * shape[2];
    if (nelem != expected) {
        throw std::runtime_error(
            "captured_io: RealArray_t nelem mismatch for '" + name + "' (got " +
            std::to_string(nelem) + ", expected " + std::to_string(expected) +
            ")");
    }
    amrex::Box sbx(amrex::IntVect(0, 0, 0),
                   amrex::IntVect(shape[0] - 1, shape[1] - 1, shape[2] - 1));
    amrex::FArrayBox fab(sbx, 1, amrex::The_Pinned_Arena());
    auto arr = fab.array();
    // Fortran column-major in the file: i fastest, then j, then k.
    for (int k = 0; k < shape[2]; ++k)
        for (int j = 0; j < shape[1]; ++j)
            for (int i = 0; i < shape[0]; ++i)
                arr(i, j, k) = read_be_f64(bin_, off);
    return fab;
}

amrex::FArrayBox CapturedFile::fab_device(const std::string& name) const {
    amrex::FArrayBox host = fab_host(name);
    amrex::FArrayBox dev(host.box(), host.nComp(), amrex::The_Arena());
    const std::size_t n = static_cast<std::size_t>(host.box().numPts()) *
                          static_cast<std::size_t>(host.nComp());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                     host.dataPtr(), host.dataPtr() + n,
                     dev.dataPtr());
    amrex::Gpu::streamSynchronize();
    return dev;
}

double CapturedFile::real64(const std::string& name) const {
    const auto& e = lookup(name, "real64");
    std::size_t off = e.offset - 1;
    return read_be_f64(bin_, off);
}

bool CapturedFile::logical(const std::string& name) const {
    const auto& e = lookup(name, "logical");
    std::size_t off = e.offset - 1;
    return read_be_i32(bin_, off) != 0;
}

} // namespace test_mom
