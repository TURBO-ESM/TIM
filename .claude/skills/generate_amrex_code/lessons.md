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

A4Box make_array4(int nx, int ny, int nz, int ncomp,
                  int lbx, int lby, int lbz);
void  free_array4(A4Box& a4);
void  copy_FortranHost_to_array4(const double* f, A4Box& a4);
void  copy_array4_to_FortranHost(const A4Box& a4, double* f);
```

Memory comes from `amrex::The_Arena()` so the same code runs on CPU and
GPU. The Fortran↔C transpose happens in a `ParallelFor` lambda *on the
device* — there is no host-side reordering. `lbx`/`lby`/`lbz` are the
array's Fortran 1-based lower bounds (`RealArray_C::lb[]`/
`LogicalArray_C::lb[]`, pass `1` for a faked-out 2D dimension, matching
how `nz=1` is passed for 2D arrays) — `make_array4` positions `a4.bx`,
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
IntA4Box make_int_array4(int nx, int ny, int nz, int ncomp,
                         int lbx, int lby, int lbz);
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
    Box bx(IntVect(bx_HOST->idxS[0]-1, bx_HOST->idxS[1]-1, bx_HOST->idxS[2]-1),
           IntVect(bx_HOST->idxE[0]-1, bx_HOST->idxE[1]-1, bx_HOST->idxE[2]-1));

    // 2. Allocate one A4Box per array, sized by RealArray_C::shape[] and
    //    positioned by RealArray_C::lb[] (pass 1 for a faked-out 2D
    //    dimension, matching how nz=1 is passed for 2D arrays -- see §2).
    //    Positioning matters even for a single-array-type kernel like this
    //    one: get it wrong and it stays invisible here (h_in/h_L/h_R all
    //    share the same lb) but silently corrupts any kernel that mixes
    //    differently-staggered arrays over one iteration box (§7 #12).
    auto h_in_DEV = turbotmp::make_array4(h_in_HOST->shape[0], h_in_HOST->shape[1], h_in_HOST->shape[2], 1,
                                          h_in_HOST->lb[0], h_in_HOST->lb[1], h_in_HOST->lb[2]);
    auto h_L_DEV  = turbotmp::make_array4(h_L_HOST->shape[0],  h_L_HOST->shape[1],  h_L_HOST->shape[2],  1,
                                          h_L_HOST->lb[0], h_L_HOST->lb[1], h_L_HOST->lb[2]);
    auto h_R_DEV  = turbotmp::make_array4(h_R_HOST->shape[0],  h_R_HOST->shape[1],  h_R_HOST->shape[2],  1,
                                          h_R_HOST->lb[0], h_R_HOST->lb[1], h_R_HOST->lb[2]);

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
- Marked `AMREX_GPU_DEVICE AMREX_FORCE_INLINE noexcept`.
- Inputs by value (`Real const h_in`), outputs by reference (`Real& h_L`).
- No `Array4`, no `Box`, no indexing — the caller supplies the cell value.
- Body is a literal port of the Fortran inner loop body.

```cpp
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
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

---

## 6. Index-base, const, and pointer conventions

- **1-based → 0-based.** Fortran `idxS[0] == 1` at the leftmost cell;
  AMReX `IntVect` is 0-based. Subtract 1 from every index when
  constructing the `Box` in the bridge. Do this in the bridge, never in
  the kernel.
- **The array, not just the iteration box, needs positioning.** The
  bridge's iteration `Box` is built at the array's true absolute
  location (`idxS-1`..`idxE-1`), but `RealArray_C`/`LogicalArray_C`
  arguments have their **own** `lb`/`ub`, independent of the iteration
  box and of each other — a staggered u/v-point array legitimately has
  a different `lb` than the h-point arrays it's called alongside. Pass
  every array's own `lb[]` into `make_array4`/`make_int_array4` (§2) so
  its `Array4` is positioned at that array's real offset; a kernel
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
    absent; pass a default-constructed `Array4` to represent absence"),
    not as a discussion of Fortran's optional-argument semantics or the
    porting decision behind the guard. If a disabled/deferred block
    needs a comment, reuse the file's existing convention for that (see
    the `NOTE: OBC support temporarily disabled` comment already used by
    `PPM_reconstruction_x`/`_y`) rather than inventing a new explanatory
    note about why this particular instance is deferred. Anything about
    *why a choice was made during this port* belongs in the commit
    message, not in the shipped source.
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
