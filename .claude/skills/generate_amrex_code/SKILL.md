---
name: generate_cpp_bridge
description: Produce the C++/AMReX side of a bridge for a MOM6 Fortran subroutine. Adds a three-tier implementation inside a TURBO-ESM/TIM checkout — an extern "C" marshalling bridge, an AMReX kernel in namespace MOM, and (when the kernel is stencil-free per cell) a pointwise device primitive — and reuses the existing turbotmp::A4Box helper for host↔device transfer and Fortran↔C layout transpose. Assumes the Fortran-side bind(C) interface and capture-mode regression input already exist; this skill does not modify any Fortran source. Mirrors the pattern established in TURBO-ESM/TIM PR #8.
user-invocable: true
argument-hint: <work-directory> <function-name>
---

# Generate C++ bridge for a MOM6 Fortran subroutine

This skill is the **execution checklist**. All templates, type tables,
conventions, and pitfalls live in [lessons.md](lessons.md) — read it
once at the start of every run (Step 1 enforces this) and refer back
to its numbered sections from each step below. Do not reproduce
templates here.

The Fortran-side wrapper (dispatcher shim, `bind(C)` interface,
caller rewrite, capture-mode plumbing) is **out of scope**. It is
assumed already merged on the MOM6 side; this skill only produces
the AMReX implementation that the Fortran shim calls into.

## Help message

If `$ARGUMENTS` is empty, or equals `help`, or equals `--help`, or
equals `-h`, do NOT run any steps. Print the following help message
verbatim and stop:

```
Usage: /generate_cpp_bridge <work-directory> <function-name>

Produce the C++/AMReX side of a bridge for a MOM6 Fortran subroutine
inside a TURBO-ESM/TIM checkout. Writes (or extends) three layers:

  - mom/cpp/<module>_kernel.hpp           pointwise device primitive
                                          (stencil-free kernels only)
  - mom/cpp/<module>.{hpp,cpp}            box-level AMReX kernel in
                                          namespace MOM
  - mom/cpp/turbotmp_<module>_bridge.{h,cpp}
                                          extern "C" marshalling bridge

Reuses turbotmp/turbotmp_helper.{hpp,cpp}; extends it only if a needed
helper is missing.

Arguments:
  <work-directory>   Absolute path to a TURBO-ESM/TIM checkout. Must
                     already exist; Step 1 populates it from `main` if
                     empty.
  <function-name>    Name of the original Fortran subroutine being
                     bridged (e.g. PPM_limit_pos). Used to derive the
                     kernel symbol, the bridge symbol, and the capture
                     filename.

Example:
  /generate_cpp_bridge /glade/derecho/scratch/sunjian/TIM PPM_limit_pos
```

## Step 0 — validate inputs

`$0` = work-directory (TIM), `$1` = function-name. Stop on the first
failure with a one-line, actionable error. Do not retry, do not assume
defaults, do not create anything.

1. **Argument count.** If `$0` empty OR `$1` empty → stop:
   `Error: missing arguments. Run "/generate_cpp_bridge --help" for usage.`
2. **Work directory exists.** If `$0` is not an existing directory → stop:
   `Error: work directory "<value>" does not exist.`
   It **may be empty** — Step 1 will clone into it.

Layout / lessons.md / plan-confirmation checks run in Step 1 after
the clone populates the tree.

## Settle these decisions (ask if not obvious from the tree)

1. **Bridge prefix** — default to whatever existing
   `turbotmp_*_bridge` declarations under `$0/mom/cpp` use; otherwise
   `turbotmp_`.
2. **Host C++ module name** (`<module>`) — lowercase the Fortran host
   module that owns `$1` (e.g. `MOM_continuity_PPM` →
   `mom_continuity_ppm`). If a bridge file for `<module>` already
   exists, **append** to it; do not create a second file.
3. **Pointwise factor** — yes if the kernel is stencil-free per cell
   (its inner loop body reads only `(i,j,k)`); no if it uses a
   stencil (`(i, j±1, k)` etc.). Controls whether Step 4 produces a
   `*_point` primitive.

If the user already specified any of these, take their values as-is.

## Procedure

Each step is one action with a pointer to the lessons.md section that
holds the template or rationale.

### 1. Clone TIM and validate
   Target branch: `main`. Apply this policy on `$0`:

   - empty → `git clone -b main git@github.com:TURBO-ESM/TIM.git $0`.
   - already a TURBO-ESM/TIM checkout (`git -C $0 remote -v`) →
     `git -C $0 fetch origin main`, then stop and surface to the user
     if any of: working tree dirty (`git status --porcelain` non-empty),
     HEAD not at `origin/main`. Do not stash, discard, or switch
     branches silently.
   - non-empty but not TURBO-ESM/TIM → stop and surface; do not
     overwrite.

   Then validate the tree:
   - `$0/mom/cpp` and `$0/turbotmp` both missing → stop:
     `Error: "<value>" is not a TURBO-ESM/TIM checkout (missing mom/cpp/ and turbotmp/).`
     One missing is OK — Steps 6–8 will create the sibling.
   - `$0/.claude/skills/generate_cpp_bridge/lessons.md` missing → stop:
     `Error: lessons.md not found at <work-directory>/.claude/skills/generate_cpp_bridge/lessons.md.`
   - **Read `lessons.md` in full.** Authoritative for naming,
     signature conventions, and the three-tier layout. Prefer it over
     this file on conflict and report the discrepancy.

   **Confirm the plan.** Print one paragraph naming: kernel symbol
   (`MOM::<lowercased_$1>`), bridge symbol
   (`<prefix>_<lowercased_$1>_bridge`), target file paths, whether
   each will be created or appended, and whether a `*_point` factor
   will be produced.

### 2. Gather the bridge contract
   Ask the user for the Fortran-side `bind(C)` interface, either as
   the literal `interface … end interface` block or as a path to a
   MOM6 source file containing
   `<prefix>_<lowercased_$1>_bridge) bind(C)` (the skill `grep`s the
   signature out).

   Derive the C prototype mechanically from the table in
   lessons.md §3.1. Print the derived C signature and pause for user
   confirmation. If any field is ambiguous, ask before proceeding.

### 3. Locate the capture file (regression input)
   Ask where `capture/<lowercased_$1>.bin` and `…meta` live. Default
   search: the user's CWD, then `$0/capture`. If not found, surface
   the absence but continue — Step 10 will skip the replay and
   record "AMReX-mode replay deferred" in Step 11.

### 4. Pointwise primitive in `mom/cpp/<module>_kernel.hpp`
   Skip if the kernel is stencil-using. Otherwise add
   `MOM::<lowercased_$1>_point(...)` per lessons.md §5. If the file
   exists, append inside the existing `namespace MOM`.

### 5. Box-level AMReX kernel in `mom/cpp/<module>.{hpp,cpp}`
   Add `MOM::<lowercased_$1>(...)` per lessons.md §4. Forward-declare
   any opaque types in the header; in the body, guard non-null
   opaque pointers with
   `AMREX_ABORT_LOC("...not yet implemented")` until the type is
   implemented. If the module file exists, append the new
   declaration in the header and the new definition in the .cpp.

### 6. Bridge header in `mom/cpp/turbotmp_<module>_bridge.h`
   First-time creation: include the mirror C structs (`RealArray_C`,
   `Box_C`), forward-declare any opaque types, wrap prototypes in
   `#ifdef __cplusplus extern "C" { … } #endif` — per lessons.md §3.
   Append the new prototype derived in Step 2. Do not re-declare the
   structs or guards if the file already exists.

### 7. Bridge implementation in `mom/cpp/turbotmp_<module>_bridge.cpp`
   Implement the new prototype following the **7-step marshalling
   pattern** in lessons.md §3 (build `Box` with index −1, `make_array4`
   each array, copy host→device for inputs **and** inouts, call the
   `MOM::` kernel, `Gpu::synchronize`, copy device→host for
   outputs/inouts only, free everything). No math, no value-dependent
   branches.

### 8. Extend `turbotmp/turbotmp_helper.{hpp,cpp}` only if needed
   Reuse existing helpers (lessons.md §2). Only add a new helper if
   the new bridge needs a layout/type that does not yet exist (e.g.
   2-D, `IntArray_C`).

### 9. Wire into the build
   Add the new sources alongside the existing `mom_continuity_ppm.cpp`
   in `CMakeLists.txt` / `Makefile.am`. Run a clean build; verify the
   resulting library exports `<prefix>_<lowercased_$1>_bridge` with
   C linkage (`nm -D --defined-only` + `grep`). On failure, surface
   the diagnostic; do not patch by changing kernel signatures or
   removing `const` without checking lessons.md §6.

### 10. Verify (replay against capture)
   If Step 3 located a capture file, run the replay procedure in
   lessons.md §8. Expected outcome:
   - **stencil-free limiter kernels** → bit-identity. Any nonzero
     diff is a bug.
   - **stencil/multi-stage kernels** → bit-identity for OBC-inactive
     configs only; OBC-active configs excluded until the opaque type
     is implemented in C++.

   If no capture file: skip and note "replay deferred" for Step 11.

### 11. Commit and push
   Branch: `claude_<lowercased_$1>_bridge` based on `main`. Stage all
   newly created and modified files. Commit message names the bridge
   symbol, the kernel symbol, whether a pointwise primitive was
   produced, and the replay result. Use a HEREDOC so the trailer is
   preserved verbatim:

   ```
   BRANCH="claude_$(echo "$1" | tr '[:upper:]' '[:lower:]')_bridge"
   git -C "$0" checkout -B "$BRANCH"
   git -C "$0" add -A
   git -C "$0" commit -m "$(cat <<'EOF'
   <one-line summary of the C++ implementation for $1>

   <1–3 line body: bridge symbol, kernel symbol,
   files added/extended, replay result>

   Co-authored-by: Claude <noreply@anthropic.com>
   EOF
   )"
   git -C "$0" push -u origin "$BRANCH"
   ```

   If the push is rejected because the branch already exists upstream
   with unrelated history, stop and surface to the user; do not
   force-push.

## Hard rules

- Bridge functions contain marshalling only — no math, no
  value-dependent branches (lessons.md §3).
- Convert Fortran 1-based indices to AMReX 0-based **in the bridge**
  (subtract 1 per axis), never in the kernel (lessons.md §6).
- Pass scalars by value with `const`; use C++ `bool` to mate with
  Fortran `logical(c_bool)`. Never C99 `_Bool`.
- Headers use `#pragma once`. No `using namespace` in headers.
- `extern "C"` lives only in the bridge `.h`. Kernel and helper
  headers are normal C++.
- Every `make_array4` reaches a `free_array4` on every path —
  including early-return / abort paths (lessons.md §7 #5).
- Reuse `turbotmp::A4Box` and its helpers; never duplicate the
  Fortran↔C transpose logic.
- Forward-declared opaque pointer types stay forward-declared until
  the model implements them; guard non-null arrivals with
  `AMREX_ABORT_LOC` (lessons.md §7 #7).
- Do not modify any Fortran source, the existing
  `turbotmp_helper.{hpp,cpp}` core API, or the mirror C struct
  layout (`RealArray_C`, `Box_C`).

If something not covered here comes up, consult lessons.md §7
(C++-side pitfalls) before improvising.

## Output to the user on success

Report:

1. Bridge symbol implemented and the kernel symbol it dispatches to.
2. List of TIM files created or extended.
3. Build result (clean, with the bridge symbol exported).
4. Replay-against-capture result, or "deferred" with the reason.
5. Branch pushed to `origin`.
6. What still needs to happen, if anything (e.g. `OceanOBC`
   implementation, OBC-active replay once the opaque type lands).
