# Lessons from TIM PR #8: AMReX/C++ implementation of the PPM bridges

Source: [TURBO-ESM/TIM PR #8](https://github.com/TURBO-ESM/TIM/pull/8) — merged commit `2374014`.

**Pre-condition:** The `generate_amrex_code` skill operates on a
pre-existing TURBO-ESM/TIM checkout. The work directory must already
contain the source tree (i.e. have `mom/cpp/` and `turbotmp/`) and must
be on (or rebased onto) the `main` branch. Cloning is not performed —
use `git clone -b main git@github.com:TURBO-ESM/TIM.git <dir>` once to
set up the directory, then pass it as `<work-directory>` on every
subsequent skill invocation.

This document distills the design, logic, and patterns used on the **C++/TIM
side** of the PPM bridge work. PR #8 adds 7 files that implement the AMReX
versions of three Fortran subroutines whose Fortran-side wrapper lives in
the MOM6 tree (`PPM_limit_pos`, `PPM_limit_cw84`, `PPM_reconstruction_y`).
The Fortran side is out of scope here; this file covers only what gets
written in the TIM repository.

The implementation is a strict **three-tier separation**:

```
   Fortran shim (in MOM6)         <-- bind(C) call across the boundary
        │
        ▼
   <prefix>_FOO_bridge   (extern "C", marshalling only, no math)
        │  copies Fortran-host arrays → device, transposes layout
        ▼
   MOM::FOO              (AMReX kernel: Box + Array4, ParallelFor)
        │
        ▼
   FOO_point             (per-cell device function, AMREX_GPU_DEVICE inline)
```

The bridge has no math; the box-level kernel has no Fortran or pointer
juggling; the pointwise primitive has no AMReX iteration logic. Each layer
can be unit-tested or replaced independently.

---

## 1. C++ file layout in TIM

| File | Role |
|---|---|
| `turbotmp/turbotmp_helper.{hpp,cpp}` | generic `A4Box`/`IntA4Box` containers + GPU arena alloc/free + Fortran↔C transpose. Reusable across all bridges. |
| `mom/cpp/turbotmp_<module>_bridge.h` | `extern "C"` prototypes of the bridge entry points; mirror C structs for `RealArray_C`, `LogicalArray_C`, `Box_C`. |
| `mom/cpp/turbotmp_<module>_bridge.cpp` | bridge implementations: marshalling only, no math. |
| `mom/cpp/<module>.hpp` | public AMReX kernel API in namespace `MOM`. |
| `mom/cpp/<module>.cpp` | box-level AMReX kernels (`ParallelFor` lambdas). |
| `mom/cpp/<module>_kernel.hpp` | per-cell device primitives `*_point`. |

Naming is symmetric with the Fortran side: bridges are co-located with the
kernels they wrap and prefixed `turbotmp_` to mark them as scaffolding —
the prefix advertises that these symbols are temporary and may be retired
once full ports are validated.

For PR #8 the concrete files were
`turbotmp_mom_continuity_ppm_bridge.{h,cpp}`,
`mom_continuity_ppm.{hpp,cpp}`, and `mom_continuity_ppm_kernel.hpp`.

---

## 2. The `turbotmp::A4Box` helper

A small struct that owns a GPU-arena allocation plus an `Array4<Real>`
view. It carries **two** device buffers per array:

- `data` — AMReX C-order, exposed via `arr` (the `Array4` consumed by
  kernels).
- `data_f` — Fortran-order staging, written/read by host↔device transfers.

```cpp
struct A4Box {
    amrex::Box                  bx;
    amrex::Real*                data   = nullptr;  // device, C order
    amrex::Array4<amrex::Real>  arr;               // view onto `data`
    amrex::Real*                data_f = nullptr;  // device, Fortran order
    int nx, ny, nz, ncomp;
};

A4Box make_array4(const RealArray_C* a);
void  free_array4(A4Box& a4);
void  copy_FortranHost_to_array4(const double* f, A4Box& a4);
void  copy_array4_to_FortranHost(const A4Box& a4, double* f);
```

Memory comes from `amrex::The_Arena()` so the same code runs on CPU and
GPU. The Fortran↔C transpose happens in a `ParallelFor` lambda *on the
device* — there is no host-side reordering. `make_array4` takes the
extents and Fortran 1-based lower bounds from the descriptor itself
(`RealArray_C::shape[]`/`lb[]`; a 2D descriptor, `dim == 2`, gets
`nz = 1` and `lbz = 1`; `ncomp` is always 1) and positions `a4.bx`,
and therefore `a4.arr`'s valid index range, at the array's **true
absolute location** (`lb-1 .. lb-1+n-1`), not at `[0, n-1]`. This
matters because the kernel is launched over the bridge's *iteration*
box (§3 step 1, built from `Box_C::idxS/idxE`, itself the array's real
absolute position) and indexes every `Array4` argument with those same
absolute coordinates — if two arrays sharing one call have different
`lb` (e.g. an h-point array vs. a staggered u/v-point array in the same
kernel), each one's `Array4` must be positioned at its own real offset
for that shared indexing to land on the correct element in both. Get
this wrong and every array with `lb == 1` is silently fine while every
staggered array reads/writes the wrong element (or out of bounds) —
see §7 #12. The staging buffer `data_f` keeps the input layout intact
while the kernel-facing `data` is reordered.

### 2.1 The `turbotmp::IntA4Box` helper (for `LogicalArray_C` / `IntArray_C`)

Byte-for-byte the same shape as `A4Box`, with `int` in place of `Real`
throughout (`amrex::Array4<int> arr`, `int* data`/`data_f`):

```cpp
IntA4Box make_int_array4(const LogicalArray_C* a);
void      free_int_array4(IntA4Box& a4);
void      copy_FortranHost_to_int_array4(const int* f, IntA4Box& a4);
void      copy_int_array4_to_FortranHost(const IntA4Box& a4, int* f);
```

Named generically (`Int`, not `Logical`) because it serves **two**
Fortran-side types that share an identical C-side layout:

- **`LogicalArray_C`** (from a Fortran `LogicalArray_t`) — the Fortran
  `logical` data is not itself C-interoperable, so `LogicalArray_t` keeps
  a separate integer-encoded (0/1) shadow buffer (`%data_c`); `%to_c()`
  fills it and returns a `LogicalArray_C` pointing at it, `%from_c()`
  reads it back into genuine `logical`s after a C/AMReX call mutates it.
  The bridge only ever sees the `int*` shadow — `LogicalArray_C.data` is
  never a pointer to real Fortran logicals.
- **`IntArray_C`** (from a Fortran `IntArray_t`) — not yet used by any
  ported kernel, but already defined Fortran-side with the same
  `data`/`shape`/`lb`/`ub`/`rank` layout, so `IntA4Box` covers it too
  without any new helper being needed when the first `IntArray_C`
  bridge parameter shows up.

At the kernel layer, a parameter fed by either type is
`Array4<const int> const&` (read-only) or `Array4<int> const&`
(writable) — 0 means false (or the integer value 0), any nonzero value
means true (or that integer value). There is no `Array4<bool>`; C++
`bool` never appears at this boundary.

### 2.2 Plain `bind(C)` value structs (no helper needed)

Some Fortran `bind(C)` types are a flat aggregate of scalars — no
`c_ptr`/shape/rank fields, nothing array-shaped — e.g.
`transport_adjust_CS_C` (3 `real(c_double)` + 5 `logical(c_bool)`,
passed `intent(in)` without `value`, so it crosses as a pointer:
`const transport_adjust_CS_C*`). These need **no** `A4Box`/`IntA4Box`-
style helper at all: there is nothing to allocate on the device and
nothing to transpose, since a `bind(C)` struct of plain scalars already
has the same memory layout the C++ side needs. The bridge just
dereferences the host pointer once and passes the resulting struct
straight through.

Because there is no bridge-vs-kernel translation step (unlike
`Box_C`→`Box` or `RealArray_C`→`Array4`), the **same struct type**
serves both layers — there is no separate "_C" bridge type and
"native" kernel type. Placement follows the existing `OceanOBC`
pattern for a module-specific type: forward-declared in the module's
`turbotmp_<module>_bridge.h` (a pointer parameter only needs a forward
declaration), with the **full** definition living in the module's
kernel header (`mom/cpp/<module>.hpp`), where kernels read its fields
directly (e.g. `CS.vol_CFL`) — never in the shared
`turbotmp_bridge_c_types.h`, which is for universal, reusable-across-
every-module types (`RealArray_C`, `Box_C`, `LogicalArray_C`), not a
single Fortran module's own config type. At the kernel layer this is
`const transport_adjust_CS_C&` (an ordinary aggregate parameter, no
different from taking `const Box&`).

---

## 3. Bridge function recipe (the 7-step marshalling pattern)

Every `<prefix>_FOO_bridge` body looks the same:

```cpp
void turbotmp_ppm_limit_pos_bridge(const Box_C* bx_HOST,
                                   const RealArray_C* h_in_HOST,
                                   RealArray_C* h_L_HOST,
                                   RealArray_C* h_R_HOST,
                                   const double h_min)
{
    // 1. Reconstruct the iteration Box (Fortran 1-based -> AMReX 0-based)
    Box bx = turbotmp::to_box(bx_HOST);

    // 2. Allocate one A4Box per array, sized by RealArray_C::shape[] and
    //    positioned by RealArray_C::lb[] (make_array4 reads both from the
    //    descriptor -- see §2).
    //    Positioning matters even for a single-array-type kernel like this
    //    one: get it wrong and it stays invisible here (h_in/h_L/h_R all
    //    share the same lb) but silently corrupts any kernel that mixes
    //    differently-staggered arrays over one iteration box (§7 #12).
    auto h_in_DEV = turbotmp::make_array4(h_in_HOST);
    auto h_L_DEV  = turbotmp::make_array4(h_L_HOST);
    auto h_R_DEV  = turbotmp::make_array4(h_R_HOST);

    // 3. Copy host -> device for every array (including inout, kernels may read first)
    turbotmp::copy_FortranHost_to_array4(h_in_HOST->data, h_in_DEV);
    turbotmp::copy_FortranHost_to_array4(h_L_HOST->data,  h_L_DEV);
    turbotmp::copy_FortranHost_to_array4(h_R_HOST->data,  h_R_DEV);

    // 4. Call the AMReX kernel
    MOM::ppm_limit_pos(bx, h_in_DEV.arr, h_L_DEV.arr, h_R_DEV.arr, h_min);

    // 5. Synchronize before reading back
    Gpu::synchronize();

    // 6. Copy device -> host for every output / inout (skip pure inputs)
    turbotmp::copy_array4_to_FortranHost(h_L_DEV, h_L_HOST->data);
    turbotmp::copy_array4_to_FortranHost(h_R_DEV, h_R_HOST->data);

    // 7. Free everything we allocated
    turbotmp::free_array4(h_in_DEV);
    turbotmp::free_array4(h_L_DEV);
    turbotmp::free_array4(h_R_DEV);
}
```

The bridge is pure marshalling — no math, no value-dependent branches.
Signature conventions:

- `const Box_C*` for the iteration domain (always the first arg).
- `const RealArray_C*` for inputs, `RealArray_C*` for outputs/inouts.
- Scalars by value: `const double h_min`, `const bool monotonic`.
- Opaque types as raw pointers: `OceanOBC* obc`.
- **`_HOST` suffix on all `Box_C*` and `RealArray_C*` parameter names** — in
  both the `.h` declaration and the `.cpp` implementation — to advertise
  that these point into Fortran host (CPU) memory, not device memory.
  Example: `const Box_C* bx_HOST`, `const RealArray_C* h_HOST`,
  `RealArray_C* dz_HOST`.

`RealArray_C` and `Box_C` are defined once in `mom/cpp/turbotmp_bridge_c_types.h`.
Bridge headers `#include "turbotmp_bridge_c_types.h"` instead of re-declaring
the structs. The bridge header then looks like:

```cpp
#pragma once
#include "turbotmp_bridge_c_types.h"

struct OceanOBC;    // forward-declared, defined elsewhere

#ifdef __cplusplus
extern "C" {
#endif
void turbotmp_ppm_limit_pos_bridge(const Box_C* bx_HOST,
                                   const RealArray_C* h_in_HOST,
                                   RealArray_C* h_L_HOST,
                                   RealArray_C* h_R_HOST,
                                   const double h_min);
// ... other prototypes ...
#ifdef __cplusplus
}
#endif
```

The struct member order in `turbotmp_bridge_c_types.h` **must** match the
Fortran-side `RealArray_C` / `Box_C` declarations field-for-field. The
Fortran field is `int :: rank`; the C field is `int dim` — only the names
differ, the layout matches.

### 3.1 Mapping a Fortran `bind(C)` interface to a C prototype

When deriving a C bridge prototype from the Fortran-side `interface …
end interface` block, every dummy maps mechanically:

| Fortran dummy declaration                       | C prototype parameter        |
|---|---|
| `type(Box_C),       intent(in)`                 | `const Box_C*`               |
| `type(RealArray_C), intent(in)`                 | `const RealArray_C*`         |
| `type(RealArray_C), intent(inout)` / `out`      | `RealArray_C*`               |
| `type(IntArray_C),  intent(in)`                 | `const IntArray_C*`          |
| `type(IntArray_C),  intent(inout)` / `out`      | `IntArray_C*`                |
| `type(LogicalArray_C), intent(in)`              | `const LogicalArray_C*`      |
| `type(LogicalArray_C), intent(inout)` / `out`   | `LogicalArray_C*`            |
| `type(<Config>_C),  intent(in)` (flat scalar struct, no `value`) | `const <Config>_C*` |
| `real(c_double),    intent(in), value`          | `const double`               |
| `integer(c_int),    intent(in), value`          | `const int`                  |
| `logical(c_bool),   intent(in), value`          | `const bool`                 |
| `type(c_ptr),       intent(in), value`          | forward-declared `<Type>*`   |

Argument order is preserved one-for-one. Names may differ; layout and
const-qualification cannot. Both `IntArray_C` and `LogicalArray_C`
marshal through the same `turbotmp::IntA4Box` helper — see §2.1. A flat
scalar-only config struct needs no helper at all and is forward-declared/
defined per the pattern in §2.2. If a `RealArray_C` arrives without a clear
intent (a writable pointer used for both input and output), it must be
copied host→device on entry — see §7 #1.

---

## 4. The AMReX kernel layer (`MOM::` namespace)

The kernel layer is the seam between AMReX and the model. Conventions:

- Public functions live in `namespace MOM` (one model-wide namespace per
  module — `MOM::ppm_limit_pos`, `MOM::PPM_reconstruction_y`).
- Iteration domain is `const Box&`.
- Read-only arrays are `Array4<const Real> const&`; writable arrays are
  `Array4<Real> const&`. A logical-mask or integer array
  (`LogicalArray_C`/`IntArray_C` on the bridge side) is
  `Array4<const int> const&` / `Array4<int> const&` — never
  `Array4<bool>` (see §2.1).
- Scalars are `const Real h_min`, `bool monotonic`, etc. — by value.
- A flat scalar-only config struct (e.g. `transport_adjust_CS_C`) is
  `const <Config>_C&` — an ordinary aggregate parameter, fields read
  directly (`CS.vol_CFL`) — never split into individual bool/double
  parameters (see §2.2).
- Opaque types (`OceanOBC*`) are forward-declared in the same header
  (`struct OceanOBC;`) and never dereferenced unless implemented. Until
  then, abort with `AMREX_ABORT_LOC("...not yet implemented")` if a
  non-null pointer arrives.
- Bodies are AMReX `ParallelFor` lambdas with
  `[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept` that delegate to
  `*_point` primitives or write inline math.
- For multi-stage kernels, intermediate buffers are
  `amrex::FArrayBox tmp_fab(box, ncomp)` allocated *inside the kernel*,
  never returned to the bridge.
- Import AMReX names with `using amrex::Box; using amrex::Array4;` at
  file scope, or `using namespace amrex;` inside the .cpp — never
  sprinkle `amrex::` everywhere, and never `using namespace` in headers.

Example shape:

```cpp
namespace MOM {
using amrex::Box;
using amrex::Array4;

void ppm_limit_pos(const Box& bx,
                   Array4<const Real> const& h_in,
                   Array4<Real>       const& h_L,
                   Array4<Real>       const& h_R,
                   const Real h_min)
{
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        ppm_limit_pos_point(h_L(i,j,k), h_R(i,j,k), h_in(i,j,k), h_min);
    });
}
}
```

---

## 5. Pointwise primitives (`*_point` functions)

The innermost, stencil-free per-cell math is factored out into a
`<module>_kernel.hpp`. Conventions:

- One function per pointwise rule, named `<kernel>_point`.
- **`_point` is always the final suffix.** For mode-split primitives use
  `<kernel>_<mode>_point` (e.g. `thickness_to_dz_3d_boussinesq_point`,
  `thickness_to_dz_3d_nonboussinesq_point`) — not `<kernel>_point_<mode>`.
  This is consistent with `ppm_limit_pos_point` / `ppm_limit_cw84_point`.
- Marked `AMREX_GPU_DEVICE AMREX_FORCE_INLINE noexcept`. Write the two
  macros together on one line, above the return type:
  `AMREX_GPU_DEVICE AMREX_FORCE_INLINE`. The same applies to
  `AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE` and to file-local device
  helpers in the `.cpp`.
- Inputs by value (`Real const h_in`), outputs by reference (`Real& h_L`).
- No `Array4`, no `Box`, no indexing — the caller supplies the cell value.
- Body is a literal port of the Fortran inner loop body.

```cpp
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void ppm_limit_pos_point(Real& h_L, Real& h_R,
                         Real const h_in, Real const h_min) noexcept
{
    Real const curv = 3.0 * ((h_L + h_R) - 2.0 * h_in);
    if (curv > 0.0) {
        Real const dh = h_R - h_L;
        if (amrex::Math::abs(dh) < curv) {
            if (h_in <= h_min) { h_L = h_in; h_R = h_in; }
            else if (12.0*curv*(h_in - h_min) < (curv*curv + 3.0*dh*dh)) {
                Real const scale = 12.0*curv*(h_in - h_min) / (curv*curv + 3.0*dh*dh);
                h_L = h_in + scale*(h_L - h_in);
                h_R = h_in + scale*(h_R - h_in);
            }
        }
    }
}
```

This split lets the same per-cell math be reused by a future CPU code
path or a unit-test driver without going through AMReX, and keeps the
device-inlined hot path clearly separated from the iteration scheduling
in `ppm_limit_pos`.

Where the original Fortran takes a stencil (e.g.,
`PPM_reconstruction_y`), the kernel keeps the AMReX `ParallelFor` (which
*is* the stencil iteration) and inlines the math directly — no `_point`
primitive is necessary unless the inner expression is used in multiple
places.

### 5.1 Fortran per-point helper procedures → one device function each

A separate case: the Fortran loop body *calls* per-point helper
procedures. They are `pure`/`elemental` subroutines or functions taking
scalars or small fixed-size stencil arrays, never grid arrays: `flux_elem`
in continuity, and CorAdCalc's 23 WENO/UP3 helpers
(`weno_seven_h_weight_reconstruction(q8, h8, u8, h_tiny, u, qr, …)` and
the weight/reconstruction helpers it calls, `fac_fn`,
`UP3_reconstruction`). These were never bridged on the Fortran side
(`convert_calltree` Step 2), so the `generate_amrex_code` invocation for
the kernel that calls them ports each one as its own `_point` device
function in `<module>_kernel.hpp`, alongside that kernel.
- **One-to-one with the Fortran helpers.** Keep the call structure
  (driver → weight/reconstruction helpers → `fac_fn`) and each helper's
  arithmetic in the same order. Bit-for-bit comparison against the
  Fortran reference depends on the order of operations.
- **Stencil inputs.** A Fortran `q8(8)` that the caller filled from a
  grid slice (`abs_vort(I,J-4:J+3,k)`) becomes either:
  - a small local array, `Real q8[8]` or `amrex::GpuArray<Real,8>`,
    filled in the `ParallelFor` lambda. This is the most literal port and
    easiest to check against the Fortran;
  - or the `Array4` plus the point's indices, reading neighbours directly.
    This avoids the copies. Choose once per helper family, before
    porting, and record it in the Step 12 report.
- **Reversed-slice stencils** (`q8(5:2:-1)`, the upwind direction) must
  keep exactly the Fortran element order: index the local array
  backwards, or pass an explicit stride. Never "simplify" them by
  symmetry.
- **Optional stencil inputs** used only under a run-time flag (the WENO
  velocity stencils `u4/u6/u8`, read only when `weno_velocity_smooth`
  is true): form and pass them only when the flag is set, testing it
  once outside the `ParallelFor` (§7 #9).
- **Naming:** `<helper>_point`, following §5's `_point`-suffix rule
  (`weno_seven_h_weight_reconstruction_point`, `fac_fn_point`).

---

## 6. Index-base, const, and pointer conventions

- **1-based → 0-based.** Fortran `idxS[0] == 1` at the leftmost cell;
  AMReX `IntVect` is 0-based. Build the bridge's `Box` with
  `turbotmp::to_box(bx_HOST)`, which subtracts 1 from every index. Do
  this in the bridge, never in the kernel.
- **The array, not just the iteration box, needs positioning.** The
  bridge's iteration `Box` is built at the array's true absolute
  location (`idxS-1`..`idxE-1`), but `RealArray_C`/`LogicalArray_C`
  arguments have their **own** `lb`/`ub`, independent of the iteration
  box and of each other — a staggered u/v-point array legitimately has
  a different `lb` than the h-point arrays it's called alongside.
  `make_array4`/`make_int_array4` take each array's own descriptor (§2)
  and position its `Array4` at that array's real offset; a kernel
  indexing several such arrays with the same shared iteration
  coordinate only lands on the correct element in each when every one
  is positioned correctly. Never rebase an array to `[0, shape-1]`
  regardless of its real `lb` — see §7 #12.
- **Const-correctness.** Bridge inputs are `const RealArray_C*` (and the
  underlying data is `const double*` once dereferenced); kernel inputs
  are `Array4<const Real> const&`. The `const`-qualification is part of
  the signature contract.
- **Pass scalars by value.** `bool` (not C99 `_Bool`), `double`, `int`
  — never references for plain scalars at the bridge boundary. This
  matches the Fortran-side `intent(in), value`.
- **Opaque pointers.** Forward-declare with `struct Foo;` in both the
  bridge header and the kernel header. Until the type is implemented,
  `AMREX_ABORT_LOC(...)` if a non-null arrives, and keep the future
  implementation behind a `/* … */` block rather than deleting it — the
  Fortran side already passes the pointer; merely accepting and
  ignoring it would silently change behaviour.
- **Headers use `#pragma once`.** No `#ifndef … #define … #endif` guards.
- **`extern "C"` only on the bridge header.** Helper and kernel headers
  are normal C++.

---

## 7. C++-side pitfalls observed in the commit log

PR #8's ~30 commits surface the friction worth flagging:

1. **Forgetting to copy inout arrays in.** A `RealArray_C*` (non-const)
   used as both input and output must still be copied host→device before
   the kernel runs — the limiter kernels read `h_L`/`h_R` *before*
   writing them. Skipping the inbound copy silently corrupts results.
2. **`bool` vs `_Bool`.** The bridge declares `const bool monotonic` to
   match the Fortran `logical(c_bool)`. C++ `bool` and Fortran
   `c_bool` are ABI-compatible; using C99 `_Bool` directly is not
   guaranteed and was specifically corrected in the commit log
   ("All scalars are now passed by value, and use bool").
3. **Namespace pollution.** Adding `using namespace amrex;` at .cpp
   scope is fine and was deliberately introduced ("Reduced the number
   of amrex:: by using namespace"), but never put `using namespace` in
   headers — it leaks into every translation unit that includes them.
4. **`#pragma once` typos.** A `#pragama once` (sic) silently disables
   the guard. Verify with `grep -n pragma` on every new header (the PR
   has a "Typo in pragma" commit).
5. **Memory leaks via early return.** Every `make_array4` must reach a
   `free_array4`. If a guard clause aborts early (e.g., OBC not
   implemented), free already-allocated containers first or — better —
   restructure to allocate only after the guards.
6. **`Gpu::synchronize` placement.** Call it after the kernel and
   before the device→host copy. Without it, the copy may race the
   kernel on async streams.
7. **OBC-not-implemented landmines.** Forward-declared `OceanOBC*` is
   ergonomic, but a non-null pointer arriving at the bridge will
   silently pass through. Guard with `AMREX_ABORT_LOC` to surface the
   case explicitly until OBC is implemented.
   `AMREX_ABORT_LOC` is **not** an AMReX macro: each kernel `.cpp` that
   uses it defines it locally, after its includes (with `<string>`):
   ```cpp
   /// @brief Abort with @p msg annotated by the source file and line number.
   /// @param msg Message to report before aborting.
   #define AMREX_ABORT_LOC(msg) \
   	amrex::Abort(std::string(msg) + " [" + __FILE__ + ":" + std::to_string(__LINE__) + "]")
   ```
   The first kernel in a module to use it (an OBC guard or a halo-update
   guard, #18) must add this block if the `.cpp` lacks it, or the build
   fails with `'AMREX_ABORT_LOC' was not declared in this scope`.
8. **AMReX include hygiene.** Bridges only need
   `<AMReX_Array4.H>`, `<AMReX_MultiFab.H>`, `<AMReX_Gpu.H>`,
   `<AMReX_FArrayBox.H>` (the last only when intermediate `FArrayBox`
   buffers are allocated). Pulling in `<AMReX.H>` umbrella headers
   slows builds without benefit.
9. **Per-cell mode flags on a call-constant bool.** If a kernel
   switches between two math paths based on a `bool` parameter that is
   constant for the entire kernel call (e.g. Boussinesq vs.
   non-Boussinesq), do **not** pass the flag inside the `ParallelFor`
   lambda. Instead create two separate `*_point` primitives — one per
   mode — and select between two separate `ParallelFor` calls with a
   single `if/else` at the box level before the loop. This eliminates
   the divergent branch from every GPU thread. (Lesson from PR #11.)
10. **Comments document physics, never the porting process.** Every doc
    comment on a ported function or parameter states what the quantity
    *physically is* — meaning, units, sign convention — normally lifted
    straight from the Fortran `!<`/`!!` doc comments or the `bind(C)`
    interface's own comments. Never add a comment that narrates *how or
    why the C++ port was done*: no explaining the 1-based→0-based index
    remapping inline (§6 already covers the mechanics; the code doing it
    doesn't also need a comment defending it), no notes on a Fortran
    naming quirk (e.g. a loop variable written as `J` in one place and
    `j` in another purely as MOM6 house style — pick one C++ name and
    move on, silently), no essay on whether a Fortran dummy is `optional`
    vs. merely guarded by `%associated()` on a mandatory argument. A
    parameter whose Fortran container may be unassociated should be
    documented as what the C++ code does about it now (e.g. "may be
    unallocated; pass a default-constructed `Array4` when unallocated"),
    not as a discussion of Fortran's optional-argument semantics or the
    porting decision behind the guard. Say "unallocated"/"allocated",
    never "absent"/"present"/"optional": the `bind(C)` interfaces have
    no `optional` dummies, only containers (`RealArray_t` arguments or
    derived-type components) that may not be allocated, which
    `%to_c()` passes as `data == nullptr`. If a disabled/deferred block
    needs a comment, reuse the file's existing convention for that (see
    the `NOTE: OBC support temporarily disabled` comment already used by
    `PPM_reconstruction_x`/`_y`) rather than inventing a new explanatory
    note about why this particular instance is deferred. Anything about
    *why a choice was made during this port* belongs in the commit
    message, not in the shipped source.
    The inverse also holds: the Fortran's own comments must not be lost.
    That includes comments in the executable body (carried over at the
    same place) and the comments on local-variable declarations
    (`real :: hvel ! hvel is the thickness used at a velocity grid
    point [H ~> m or kg m-2]`). When a Fortran local becomes a C++
    local, whether a scalar, an `FArrayBox` scratch field or the
    `Array4` view of one, its description and units go in a `///`
    comment directly above the C++ declaration; for a scratch field put
    it above the `FArrayBox`. Locals with no C++ counterpart (`%view`
    pointers, loop indices) drop their comments with them. A port that
    keeps every body comment but declares `hvel`, `dz_harm` and
    `a_cpl_gl90` bare loses the units the reader most needs (seen in
    the first `mom_vert_friction.cpp` port).
11. **Disabled/deferred blocks get real code, one level of disabling.**
    When Fortran logic beyond an `AMREX_ABORT_LOC` guard (typically an
    OBC segment loop) is deferred rather than ported live, wrap it in a
    single `/* ... */` block containing **genuine, faithfully-ported
    C++** — a real `ParallelFor` over a real derived `Box`, real
    assignment statements — exactly like `PPM_reconstruction_y`'s
    existing disabled OBC block (the one overriding `h_S`/`h_N` near
    the end of `mom_continuity_ppm.cpp`). Do not additionally prefix
    the inner statements with `//`: the enclosing `/* */` already
    disables compilation, so a second layer of line-commenting is
    redundant, and it silently degrades a real port into an inert
    stub that reads as unfinished even after the real translation was
    done. (A `for`/`if` skeleton with commented-out assignments inside
    it was exactly this mistake, caught during PR review of the
    `meridional_flux_thickness`/`zonal_flux_thickness` ports — both
    had to be fixed to match `PPM_reconstruction_y`'s one-level
    pattern.) One level of disabling only, real code inside it.
12. **`make_array4`/`make_int_array4` must be positioned by the array's
    own `lb`, not just sized by `shape[]`.** Discovered when 8 kernels
    in the `continuity_PPM` chain — all ones combining h-point and
    staggered u/v-point arrays in one call — crashed or produced
    non-bit-for-bit results in real MOM6 integration while their
    GoogleTest unit tests all passed. The unit tests call `MOM::kernel`
    directly with correctly-positioned `Array4`s (built from
    `CapturedFile`, which does honor real `lb`), so they never exercise
    the bridge's marshalling at all — this bug is structurally
    invisible to a kernel-level unit test and only shows up once the
    real Fortran caller drives the actual `turbotmp_*_bridge` entry
    point. If a kernel-level test suite is fully green but real
    integration crashes or drifts non-bit-for-bit on a kernel that
    takes both h-point and staggered (u/v-point) arrays, suspect this
    exact class of bug first: check whether `make_array4`/
    `make_int_array4` are being called with the array's real `lb[]`
    (§2, §6), not just `shape[]`.

---
13. **A Fortran output with no bridge channel.** `zonal_mass_flux`'s
    Fortran body ends by writing `BT_cont%h_u` via
    `zonal_flux_thickness`, but its `bind(C)` interface carried only
    the six `FA_u_*`/`uBT_*` BT_cont fields. The C++ port noticed and
    left a source comment calling `h_u` "out of scope" because the
    bridge "does not currently expose an h_u channel" — consistent with
    this skill's no-Fortran-edits rule, but the gap was never reported
    back, so AMReX mode silently returned a stale `h_u`. A missing
    channel is a bridge-contract bug on the MOM6 side: stop and report
    it (SKILL.md Step 3) rather than documenting it in code.
14. **An orchestrator skipping a shim's argument transformation.**
    `continuity_PPM_fortran` calls `zonal_edge_thickness(..., Angstrom_H,
    ...)`; that shim computes `h_min = 2.0 * Angstrom_H` and passes the
    doubled value to its bridge and `_fortran` worker. The C++
    `MOM::continuity_PPM` called `MOM::zonal_edge_thickness` directly
    with `Angstrom_H`, so the PPM reconstruction ran with half the
    minimum thickness (`continuity_PPM_2d_fluxes` had the same bug).
    Every leaf kernel matched Fortran bit-for-bit when run individually
    — each shim still did its own doubling — and only
    `CONTINUITY_PPM_MODE=AMREX` diverged, and only where an edge value
    reached the thickness floor. Fix: for every callee the C++
    orchestrator invokes directly, read the callee's Fortran shim and
    reproduce each scaled, derived, or defaulted argument (SKILL.md
    Step 6).
15. **Comments use the C++ Doxygen form `///`, never `/** ... */`.**
    Every comment written into a generated C++ file — the `@file`/
    `@brief` header, function and parameter doc blocks (`@param`),
    struct-member docs (`///<`), and in-body step comments — uses
    `///`. Model new files on `mom_continuity_ppm*`, which use this
    style throughout. Do not copy the `/** ... */` blocks in
    `mom_interface_heights*` or in the pre-existing header of
    `turbotmp_bridge_c_types.h`. Leave those existing comments as
    they are; the rule covers only text this skill writes. The gradKE
    bridge was first written with `/** */` blocks copied from the
    interface_heights bridge and had to be converted. The one
    allowed `/* ... */` is the single disabling wrapper around a
    deferred block of real code (#11), which is not a comment.
16. **Array4 bounds come from `Box(arr)`, never from `arr.begin.x`.**
    Since AMReX 26.x, `Array4::begin`/`end` are `IntVectND<4>`, not
    `Dim3`. They have no `.x`/`.y`/`.z` members, so code like
    `IntVect(q.begin.x, q.begin.y, k)` no longer compiles. To size a
    scratch `FArrayBox` from an argument array, use
    `Box(arr).smallEnd(d)` / `Box(arr).bigEnd(d)`. When the whole box
    is wanted, use `Box(arr)` itself, as `mom_continuity_ppm.cpp` does
    (commit 28b0798). `Box(arr).bigEnd(d)` already equals the old
    `arr.end.x - 1`, so do not subtract 1 again. Both the PPM
    continuity port and the CorAdCalc port (`coradv_setup`,
    `coradv_arakawa`, `coradv_sadourny`, `coradcalc_tr`) hit this
    build error. The
    `Array4(p, Dim3 begin, Dim3 end, ncomp)` constructor, and the
    `.p` member, are still available in 26.x.
17. **Every box-level kernel starts with `BL_PROFILE`.** The first
    statement of each public `MOM::<kernel>` in `mom/cpp/<module>.cpp`
    is `BL_PROFILE("<kernel>");`, with the kernel's own C++ name as the
    label, followed by a blank line. It goes before any OBC guard, so
    the timer covers the whole call. Include `<AMReX.H>` in the
    `.cpp`. The `_point` device helpers, file-local templates in the
    `.cpp` and the bridge functions get no profiling call. This is the
    convention in `mom_continuity_ppm.cpp` (e.g.
    `BL_PROFILE("ppm_limit_pos")`) and `mom_coriolisadv.cpp` (e.g.
    `BL_PROFILE("coradv_setup")`). The `mom_vert_friction` port first
    shipped without these calls because this checklist did not ask
    for them, which left its kernels out of the AMReX TinyProfiler
    output.
18. **Halo updates (`pass_var` / `pass_vector`) abort, everything else is
    ported.** There is no tested C++ halo exchange yet. In MOM6 these calls
    always sit inside a run-time flag (`if (CS%...) then ... call
    pass_var(...) ... end if`), and the initial test configurations never
    take that branch. Port the `if` faithfully and make the first
    statement of its body
    `AMREX_ABORT_LOC("<kernel>: halo update (pass_var) not yet implemented");`
    (name `pass_vector` when that is the call). Port the rest of the branch
    as live, compiled C++. In place of each halo-update call, put the
    original Fortran call as a single `///` comment, e.g.
    `/// call pass_var(u, G%Domain)`. This comment is a deliberate
    exception to #10 and #11: it marks exactly where the missing exchange
    belongs. Reaching the abort at run time means it's time to implement
    the exchange. If a `pass_var`/`pass_vector` call is ever found
    **outside** a flag, stop and report it rather than guarding it.
19. **Floating-point literals carry the `_rt` suffix.** Write every
    floating-point literal in kernel code as an `amrex::Real` literal:
    `0.25_rt`, `1.0_rt`, `1.e-30_rt`, `0._rt`. Leave integer literals
    used for indices and counts alone.
    - **Why it matters:** in a double-precision build, `Real` is `double`, so
      this changes nothing and results stay bit-for-bit identical. In a
      single-precision AMReX build, a bare `0.25` is a `double`. That
      promotes the expression to double, and mixed calls such as
      `amrex::max(x, 0.0)` fail template deduction and do not compile.
    - **How to enable it:** in a `.cpp`, add
      `using namespace amrex::literals;` after `using namespace amrex;`. In
      a header (`<module>_kernel.hpp`), put that using-directive inside each
      function body that needs it, never at file or namespace scope. The
      hard rule against `using namespace` in headers still holds.
      `mom_continuity_ppm` set the convention, although its
      `mom_continuity_ppm_kernel.hpp` uses a file-scope directive that this
      rule now forbids.
    - **Where it slipped:** the CorAdCalc, vertvisc and interface-heights
      ports used bare literals because this list did not ask for `_rt`.
      `mom_hor_visc` was converted afterwards.
20. **Plain fills and copies use the `BaseFab` member functions, not
    hand-written loops.** Wrap the `Array4` in a non-owning
    `FArrayBox(arr)` and call the AMReX member function. The wrapper
    neither allocates nor copies; it just views the same memory.

    | Instead of a `ParallelFor` that … | Write |
    |---|---|
    | sets every point of `arr` to `v` | `FArrayBox(arr).setVal<RunOn::Device>(v);` |
    | sets the points of `bx` to `v` | `FArrayBox(arr).setVal<RunOn::Device>(v, bx);` |
    | copies all of `src` into a scratch FAB on `Box(src)` | `dst_fab.copy<RunOn::Device>(FArrayBox(src));` |

    - **Scope:** this is only for loops whose body is a single assignment of
      one value, or a straight copy.
    - **Keep the box:** pass the same box the loop would have iterated over,
      so results are unchanged. That box must lie inside the array's extent.
    - **Comments:** any Fortran comment on the original assignment goes on a
      `///` line above the call.
    - **When to keep the loop:** a loop that sets a value and then computes
      in the same pass stays a `ParallelFor`, for example
      `visc_bound_rem = 1.0_rt` followed by the bounding logic. So does a
      copy between differently shaped boxes, such as one layer of a 3-D
      block into a 2-D buffer.
    - **Where it slipped:** `mom_hor_visc` first shipped 21 hand-written
      zero-fill loops and a `copy_prev_field` helper. Both were converted
      afterwards.

---

## 8. C++ verification harness (replay against captures)

The capture files written by the Fortran shim — `capture/<kernel>.bin`
and `capture/<kernel>.meta` produced when the Fortran side runs in its
CAPTURE mode — are the C++ side's regression input. A C++ test driver:

1. Reads `capture/<kernel>.{bin,meta}`.
2. Materialises Fortran-layout host buffers from the recorded inputs
   (`_<argname>` / `_<argname>_before`).
3. Calls `<prefix>_FOO_bridge(...)` exactly as the Fortran shim would.
4. Diffs the resulting buffers against the recorded
   `_<argname>_after` arrays.

Bit-identity is the goal for the limiter kernels (no FP reordering);
for `PPM_reconstruction_y` the OBC branches are currently disabled, so
any comparison run must exclude OBC-active configurations until the
forward-declared `OceanOBC` type is implemented in C++.

## 9. Backlog

**Orchestrator kernels calling other already-ported kernels on
internal scratch (e.g. `zonal_mass_flux`/`meridional_mass_flux` calling
`zonal_flux_adjust`/`set_zonal_BT_cont`/`zonal_flux_thickness`):**
currently ported as separate kernel launches over locally-allocated
scratch `FArrayBox`/`IArrayBox` buffers (mirroring the Fortran
orchestrator's own `allocView`/`free` calls) — the direct, low-risk
translation of "allocate local arrays, call other subroutines." A
faster but more invasive alternative, deferred for now: extract the
Newton-iteration body of `zonal_flux_adjust`/`meridional_flux_adjust`
into a shared per-thread device function (like `flux_elem_point`) that
the orchestrator's own per-thread body calls directly, keeping
everything as per-thread scalars and avoiding scratch allocation and
the extra kernel launches entirely. Revisit once the option-1 ports are
working and there's a reason to care about the extra launches (e.g. a
profiling result), since option 2 means refactoring
`zonal_flux_adjust`/`meridional_flux_adjust`'s already-approved kernel
bodies rather than just adding new code.

The same shared-per-thread-device-function pattern is now the planned
approach for per-point helper families (§5.1). CorAdCalc's WENO/UP3
helpers were decided that way from the start, so the Fortran helpers
are never bridged individually.
