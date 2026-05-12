---
name: generate_amrex_code
description: Produce the C++/AMReX side of a bridge for a MOM6 Fortran subroutine. Adds a three-tier implementation inside a TURBO-ESM/TIM checkout — an extern "C" marshalling bridge, an AMReX kernel in namespace MOM, and (when the kernel is stencil-free per cell) a pointwise device primitive — and reuses the existing turbotmp::A4Box helper for host↔device transfer and Fortran↔C layout transpose. Grounds the AMReX port in the original Fortran source, located either via an optional MOM6 checkout path or by user paste. Assumes the Fortran-side bind(C) interface and capture-mode regression input already exist; this skill does not modify any Fortran source. Mirrors the pattern established in TURBO-ESM/TIM PR #8.
user-invocable: true
argument-hint: <work-directory> <function-name> [<mom6-directory>] [--enable_src_validate] [--enable_git_commit]
---

# Generate C++/AMReX implementation for a MOM6 Fortran subroutine

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
Usage: /generate_amrex_code <work-directory> <function-name> [<mom6-directory>] [--enable_src_validate] [--enable_git_commit]

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
  <work-directory>     Absolute path to an existing TURBO-ESM/TIM checkout
                       (must contain mom/cpp/ and turbotmp/, and should be
                       on the main branch). Cloning is not performed.
  <function-name>      Name of the original Fortran subroutine being
                       bridged (e.g. PPM_limit_pos). Used to derive
                       the kernel symbol, the bridge symbol, and the
                       capture filename.
  <mom6-directory>     OPTIONAL. Absolute path to a TURBO-ESM/MOM6
                       checkout. When given, Step 3 grep's the
                       Fortran subroutine body out of src/ and
                       config_src/ to ground the AMReX port. When
                       omitted, Step 3 asks the user to paste the
                       Fortran source body.
  --enable_src_validate  (optional) Run Step 1: verify the work directory is a
                         TURBO-ESM/TIM checkout on the main branch and that the
                         tree layout is valid. Off by default.
  --enable_git_commit    (optional) Run Step 12: create branch
                         claude_<function-name>_bridge, commit all changes, and
                         push to origin. Off by default.

Example:
  /generate_amrex_code /glade/derecho/scratch/sunjian/TIM PPM_limit_pos \
                       /glade/derecho/scratch/sunjian/MOM6
```

## Step 0 — validate inputs

`$0` = work-directory (TIM), `$1` = function-name, `$2` =
mom6-directory (optional). Stop on the first failure with a
one-line, actionable error. Do not retry, do not assume defaults,
do not create anything.

1. **Argument count.** If `$0` empty OR `$1` empty → stop:
   `Error: missing arguments. Run "/generate_amrex_code --help" for usage.`
2. **TIM work directory is an existing checkout.** If `$0` is not an
   existing directory → stop:
   `Error: work directory "<value>" does not exist.`
   The directory must already contain a TURBO-ESM/TIM checkout —
   cloning is not performed by this skill. Step 1 verifies the checkout
   identity and validates the tree layout.
3. **MOM6 directory (optional).** If `$2` was provided but is not
   an existing directory → stop:
   `Error: MOM6 directory "<value>" does not exist.` If `$2` was
   omitted, record `mom6_mode=paste` and Step 3 will ask the user
   for the Fortran source body inline.
4. **Parse optional flags.** Scan remaining arguments for `--enable_src_validate`
   and `--enable_git_commit`. Store as boolean flags (default: off). Any
   unrecognised argument that starts with `--` → stop:
   `Error: unknown option "<value>". Run "/generate_amrex_code --help" for usage.`

Step 1 runs only when `--enable_src_validate` is set; Step 12 runs only when
`--enable_git_commit` is set. Layout / lessons.md / plan-confirmation checks
run in Step 1.

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
   stencil (`(i, j±1, k)` etc.). Controls whether Step 5 produces a
   `*_point` primitive. The Fortran source captured in Step 3 is
   what determines this — confirm or correct after Step 3.

If the user already specified any of these, take their values as-is.

## Procedure

Each step is one action with a pointer to the lessons.md section that
holds the template or rationale.

### 1. Validate the existing checkout *(runs only when `--enable_src_validate` is passed; skip otherwise and proceed to Step 2)*
   If `--enable_src_validate` was not supplied, skip this entire step and go to Step 2.

   **Verify the work directory.** All bridge work in this skill is based
   on the `main` branch — that is the agreed base and what the Fortran
   side expects to merge against.

   Verify `$0` is a TURBO-ESM/TIM checkout by running
   `git -C $0 remote -v` and confirming the output contains
   `TURBO-ESM/TIM`. If not → stop:
   `Error: "<value>" is not a TURBO-ESM/TIM checkout.`

   Then run `git -C $0 fetch origin main` and check the branch state:
   - If the working tree is dirty (`git -C $0 status --porcelain` is
     non-empty) → stop and surface to the user; do not stash or discard.
   - If HEAD is not already at `origin/main` (compare
     `git -C $0 rev-parse HEAD` and `git -C $0 rev-parse origin/main`)
     → stop and ask the user whether to `git checkout main` before
     continuing. Do not switch branches silently.

   Then validate the tree:
   - `$0/mom/cpp` and `$0/turbotmp` both missing → stop:
     `Error: "<value>" is not a TURBO-ESM/TIM checkout (missing mom/cpp/ and turbotmp/).`
     One missing is OK — Steps 7–9 will create the sibling.
   - `$0/.claude/skills/generate_amrex_code/lessons.md` missing → stop:
     `Error: lessons.md not found at <work-directory>/.claude/skills/generate_amrex_code/lessons.md.`
   - **Read `lessons.md` in full.** Authoritative for naming,
     signature conventions, and the three-tier layout. Prefer it over
     this file on conflict and report the discrepancy.

   **Confirm the plan.** Print one paragraph naming: kernel symbol
   (`MOM::<lowercased_$1>`), bridge symbol
   (`<prefix>_<lowercased_$1>_bridge`), target file paths, whether
   each will be created or appended, and whether a `*_point` factor
   is the working assumption (Step 3 may override it).

### 2. Gather the bridge contract
   Ask the user for the Fortran-side `bind(C)` interface, either as
   the literal `interface … end interface` block or as a path to a
   MOM6 source file containing
   `<prefix>_<lowercased_$1>_bridge) bind(C)` (the skill `grep`s the
   signature out — if `$2` was supplied, default the search root to
   `$2/src` and `$2/config_src`).

   Derive the C prototype mechanically from the table in
   lessons.md §3.1. Print the derived C signature and pause for user
   confirmation. If any field is ambiguous, ask before proceeding.

### 3. Locate and read the Fortran kernel source
   The AMReX kernel body in Steps 5–6 must be a faithful port of the
   Fortran subroutine body, so the skill must read that body before
   it writes any C++.

   - **If `$2` was supplied**, run
     `grep -irln "^[[:space:]]*subroutine[[:space:]]\+\(<function-name>_fortran\|<function-name>\)\b" $2/src $2/config_src`.
     - 1 match → open the file, locate the
       `subroutine … end subroutine` block, and capture the body for
       Steps 5–6 along with its file path and line range.
     - 0 matches → fall through to user paste below.
     - >1 matches → list the candidate file paths and ask the user
       which to use. Do not pick silently.
   - **If `$2` was omitted, or the grep found nothing**, ask the
     user to paste the Fortran subroutine body inline. Either the
     `<function-name>_fortran` form (post-wrap) or the original
     pre-wrap form is acceptable; record which.

   After capture, print a one-line acknowledgement: e.g.
   `Loaded body from <path>:<lineS>-<lineE> (NN lines).` or
   `Loaded body from user paste (NN lines).`

   Then re-evaluate "Settle these decisions" #3 (pointwise factor):
   if the captured body reads any neighbour index (`j-1`, `j+1`,
   etc.), set the factor to `no` and record the change. If the body
   is purely `(i,j,k)`-local, keep `yes`.

   Do not proceed to Step 4 until the body has been captured.

### 4. Locate the capture file (regression input)
   Ask where `capture/<lowercased_$1>.bin` and `…meta` live. Default
   search: the user's CWD, then `$0/capture`. If not found, surface
   the absence but continue — Step 11 will skip the replay and
   record "AMReX-mode replay deferred" in Step 12.

### 5. Pointwise primitive in `mom/cpp/<module>_kernel.hpp`
   Skip if "Settle these decisions" #3 (as confirmed at the end of
   Step 3) said "no pointwise factor". Otherwise add
   `MOM::<lowercased_$1>_point(...)` per lessons.md §5, porting the
   Fortran inner loop body captured in Step 3 verbatim (preserve
   parenthesisation to keep floating-point evaluation order). If the
   file exists, append inside the existing `namespace MOM`.

   **Mode flags:** If the Fortran inner loop contains an `if`-test on
   a bool that is constant for the entire loop (a mode flag such as
   Boussinesq), create **two** separate `_point` functions — one per
   mode — rather than a single function with a flag parameter. Name
   them `<kernel>_point_<mode>` (e.g.
   `thickness_to_dz_3d_point_boussinesq` /
   `thickness_to_dz_3d_point_nonboussinesq`). See lessons.md §7 #9.

### 6. Box-level AMReX kernel in `mom/cpp/<module>.{hpp,cpp}`
   Add `MOM::<lowercased_$1>(...)` per lessons.md §4. The body of the
   `ParallelFor` lambda is either a one-line call to the `*_point`
   primitive from Step 5 (stencil-free case) or a direct port of the
   Fortran loop nest from Step 3 (stencil case). Forward-declare any
   opaque types in the header; in the body, guard non-null opaque
   pointers with `AMREX_ABORT_LOC("...not yet implemented")` until
   the type is implemented. If the module file exists, append the
   new declaration in the header and the new definition in the .cpp.

   **Mode dispatch:** When two `_point` primitives exist for distinct
   modes (Step 5), select between them with a single `if/else`
   **before** the `ParallelFor`, never inside the lambda. Each lambda
   calls exactly one primitive with no conditional. See lessons.md §7 #9.

### 7. Bridge header in `mom/cpp/turbotmp_<module>_bridge.h`
   First-time creation: start with `#include "turbotmp_bridge_c_types.h"`
   (do **not** re-define `RealArray_C` / `Box_C` — they live in that shared
   file). Forward-declare any opaque types, wrap prototypes in
   `#ifdef __cplusplus extern "C" { … } #endif` — per lessons.md §3.
   Append the new prototype derived in Step 2. Do not re-declare the
   structs or guards if the file already exists.
   All `Box_C*` and `RealArray_C*` parameter names must carry the `_HOST`
   suffix in the declaration here, matching the `.cpp` implementation —
   see lessons.md §3.

### 8. Bridge implementation in `mom/cpp/turbotmp_<module>_bridge.cpp`
   Implement the new prototype following the **7-step marshalling
   pattern** in lessons.md §3 (build `Box` with index −1, `make_array4`
   each array, copy host→device for inputs **and** inouts, call the
   `MOM::` kernel, `Gpu::synchronize`, copy device→host for
   outputs/inouts only, free everything). No math, no value-dependent
   branches.

### 9. Extend `turbotmp/turbotmp_helper.{hpp,cpp}` only if needed
   Reuse existing helpers (lessons.md §2). Only add a new helper if
   the new bridge needs a layout/type that does not yet exist (e.g.
   2-D, `IntArray_C`).

### 10. Wire into the build
   Add the new sources alongside the existing `mom_continuity_ppm.cpp`
   in `CMakeLists.txt` / `Makefile.am`. Run a clean build; verify the
   resulting library exports `<prefix>_<lowercased_$1>_bridge` with
   C linkage (`nm -D --defined-only` + `grep`). On failure, surface
   the diagnostic; do not patch by changing kernel signatures or
   removing `const` without checking lessons.md §6.

### 11. Verify (replay against capture)
   If Step 4 located a capture file, run the replay procedure in
   lessons.md §8. Expected outcome:
   - **stencil-free limiter kernels** → bit-identity. Any nonzero
     diff is a bug.
   - **stencil/multi-stage kernels** → bit-identity for OBC-inactive
     configs only; OBC-active configs excluded until the opaque type
     is implemented in C++.

   If no capture file: skip and note "replay deferred" for Step 12.

### 12. Commit and push *(runs only when `--enable_git_commit` is passed; skip otherwise)*
   If `--enable_git_commit` was not supplied, skip this entire step and report the
   files that were modified so the user can commit manually.

   Branch: `claude_<lowercased_$1>_bridge` based on `main`. Stage all
   newly created and modified files. Commit message names the bridge
   symbol, the kernel symbol, whether a pointwise primitive was
   produced, the source-of-truth for the port (file path + line
   range, or "user paste"), and the replay result. Use a HEREDOC so
   the trailer is preserved verbatim:

   ```
   BRANCH="claude_$(echo "$1" | tr '[:upper:]' '[:lower:]')_bridge"
   git -C "$0" checkout -B "$BRANCH"
   git -C "$0" add -A
   git -C "$0" commit -m "$(cat <<'EOF'
   <one-line summary of the C++ implementation for $1>

   <1–3 line body: bridge symbol, kernel symbol,
   files added/extended, source-of-truth, replay result>

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
- Do not invent kernel math. The AMReX body must come from the
  Fortran source captured in Step 3, ported verbatim.

If something not covered here comes up, consult lessons.md §7
(C++-side pitfalls) before improvising.

## Output to the user on success

Report:

1. Bridge symbol implemented and the kernel symbol it dispatches to.
2. Source-of-truth for the port: file path + line range (if grep'd
   from `$2`) or "user paste".
3. List of TIM files created or extended.
4. Build result (clean, with the bridge symbol exported).
5. Replay-against-capture result, or "deferred" with the reason.
6. Branch pushed to `origin`.
7. What still needs to happen, if anything (e.g. `OceanOBC`
   implementation, OBC-active replay once the opaque type lands).
