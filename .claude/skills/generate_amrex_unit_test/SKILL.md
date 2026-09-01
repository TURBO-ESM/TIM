---
name: generate_amrex_unit_test
version: "0.3.1"
description: Add (or extend) GoogleTest coverage under test_mom/ for an already-implemented MOM::<kernel> AMReX function inside a TURBO-ESM/TIM checkout, replaying it against the Fortran-captured regression fixture (capture/<kernel>.{bin,meta}) via the existing test_mom::CapturedFile / expect_arrays_equal harness. Assumes the box-level AMReX kernel (mom/cpp/<module>.{hpp,cpp}) already exists -- produced by generate_amrex_code; this skill does not write kernel or bridge code, only test_mom/ test files and CMake wiring. When no capture fixture exists yet, grounds the field mapping in the Fortran shim's capture-mode rec%add(...) calls (via an optional MOM6 checkout path or by user paste) so a real assertion can still be written, falling back to a GTEST_SKIP() stub only when neither a fixture nor Fortran source is available. Mirrors the pattern established by the existing test_mom/mom_continuity_ppm/test_mom_continuity_ppm.cpp tests.
user-invocable: true
argument-hint: <work-directory> <function-name> [<mom6-directory>] [--enable_src_validate] [--enable_build] [--enable_git_commit] [--disable_git_commit]
---

# Generate a GoogleTest unit test for a MOM AMReX kernel

This skill is the **execution checklist**; there is no companion
lessons file — the reference material (the `test_mom/` layout, the
`CapturedFile` API, the naming conventions) is small enough to live
inline in each step below.

**Out of scope:** writing or modifying the kernel itself
(`mom/cpp/<module>.{hpp,cpp}`), the pointwise primitive, the bridge
(`turbotmp_<module>_bridge.{h,cpp}`), or any Fortran source. All of
that is `generate_amrex_code`'s job and is assumed already merged —
Step 2 stops if the kernel isn't there yet. This skill only adds to
`test_mom/`: a `TEST(...)` case, and CMake wiring if the module has no
test executable yet.

**These tests bypass the bridge.** Every `TEST(...)` this skill writes
calls `MOM::<kernel>(...)` directly, in-process, with `Array4`s built
straight from `CapturedFile` — never through
`turbotmp_<module>_bridge.{h,cpp}`'s `extern "C"` marshalling (the
`Box_C`/`RealArray_C` conversion, `make_array4` sizing/positioning,
host↔device copies) that the real Fortran caller actually goes
through. A green suite here proves the kernel logic matches the
captured Fortran output; it proves nothing about the marshalling
layer. Concretely: `generate_amrex_code`'s lessons.md §7 #12 documents
a real case where 8 kernels' bridges mis-positioned marshalled arrays
(sized by `shape[]` but not positioned by `lb[]`), causing crashes and
non-bit-for-bit results in live MOM6 integration, while every one of
these kernel-level unit tests stayed green throughout — because the
bug lived entirely in code this skill's tests never call. Do not treat
a passing suite here as evidence the bridge is correct; that requires
either a real integration run or a test that goes through
`turbotmp_<module>_bridge` itself (not something this skill currently
produces).

A MOM6 checkout is **read-only reference material**, used only when no
capture fixture exists yet (Step 3) to ground the field-name mapping in
the Fortran shim's `rec%add(...)` calls — never to modify Fortran
source, and never required when a fixture is already on disk.

## Help message

If `$ARGUMENTS` is empty, or equals `help`, or equals `--help`, or
equals `-h`, do NOT run any steps. Print the following help message
verbatim and stop:

```
Usage: /generate_amrex_unit_test <work-directory> <function-name> [<mom6-directory>] [--enable_src_validate] [--enable_build] [--enable_git_commit] [--disable_git_commit]

Add GoogleTest coverage under test_mom/ for an already-implemented
MOM::<kernel> AMReX function, replayed against its Fortran-captured
regression fixture. Does not write kernel or bridge code.

Arguments:
  <work-directory>     Absolute path to an existing TURBO-ESM/TIM checkout
                       (must contain mom/cpp/ and test_mom/, and should be
                       on the main branch). Cloning is not performed.
  <function-name>      Name of the original Fortran subroutine whose
                       AMReX kernel (MOM::<lowercased_name>) is being
                       tested (e.g. PPM_limit_pos). Used to derive the
                       kernel symbol, the test name, and the capture
                       fixture's base filename.
  <mom6-directory>     OPTIONAL. Absolute path to a TURBO-ESM/MOM6
                       checkout. Only consulted when Step 3 finds no
                       capture fixture on disk: Step 3 then greps the
                       Fortran shim's CAPTURE-mode rec%add(...) calls
                       out of src/ and config_src/ to ground the field
                       mapping (read-only -- no Fortran source is ever
                       modified). When omitted and no fixture exists,
                       Step 3 asks the user to paste that block instead.
  --enable_src_validate  (optional) Run Step 1: verify the work directory is a
                         TURBO-ESM/TIM checkout on the main branch and that the
                         test_mom/ scaffold (common/, test_mom_main.cpp) is
                         present. Off by default.
  --enable_build          (optional) Run Steps 9-10: configure, build (from
                         source, including AMReX if no install is found --
                         this can take a long time on a first build), and
                         ctest the new test. Off by default -- without it,
                         this skill stops after writing the harness (Step 8)
                         and reports the build as not run.
  --enable_git_commit    (optional) Force Step 11 to run this invocation
                         only, overriding the global git_commit_and_push
                         preference in ~/.claude/preferences.json.
                         Mutually exclusive with --disable_git_commit.
  --disable_git_commit   (optional) Force Step 11 to be skipped this
                         invocation only, overriding the global
                         preference. Mutually exclusive with
                         --enable_git_commit.

When neither --enable_git_commit nor --disable_git_commit is passed, Step 11
follows ~/.claude/preferences.json's "git_commit_and_push" key ("auto" or
"manual"; treated as "manual" when the file is missing, the key is absent,
or the JSON fails to parse). There is no equivalent global preference for
--enable_build -- it is off by default every run, no exceptions, because a
first build compiles AMReX from source and can take a long time.

Example:
  /generate_amrex_unit_test /glade/derecho/scratch/sunjian/TIM PPM_limit_pos \
                            /glade/derecho/scratch/sunjian/MOM6
```

## Step 0 — validate inputs

`$0` = work-directory (TIM), `$1` = function-name, `$2` =
mom6-directory (optional). Stop on the first failure with a one-line,
actionable error. Do not retry, do not assume defaults, do not create
anything.

1. **Argument count.** If `$0` empty OR `$1` empty → stop:
   `Error: missing arguments. Run "/generate_amrex_unit_test --help" for usage.`
2. **TIM work directory is an existing checkout.** If `$0` is not an
   existing directory → stop:
   `Error: work directory "<value>" does not exist.`
   Cloning is not performed by this skill. Step 1 verifies checkout
   identity and validates the tree layout.
3. **MOM6 directory (optional).** If `$2` was provided but is not an
   existing directory → stop:
   `Error: MOM6 directory "<value>" does not exist.` If `$2` was
   omitted, record `mom6_mode=paste` — Step 3 will ask the user to
   paste the Fortran shim's capture block inline if it turns out to be
   needed (no fixture found). If a fixture is found, `$2` is never
   consulted either way.
4. **Parse optional flags.** Scan remaining arguments for
   `--enable_src_validate`, `--enable_build`, `--enable_git_commit`,
   and `--disable_git_commit`. Store as boolean flags (default: off).
   If both `--enable_git_commit` and `--disable_git_commit` are passed
   → stop: `Error: --enable_git_commit and --disable_git_commit are
   mutually exclusive.` Any unrecognised argument that starts with
   `--` → stop:
   `Error: unknown option "<value>". Run "/generate_amrex_unit_test --help" for usage.`

Step 1 runs only when `--enable_src_validate` is set. Steps 9-10 run
only when `--enable_build` is set — there is no global-preference
fallback for it; omitting the flag always means "stop after writing
the harness." Step 11's behavior is decided by
`--enable_git_commit`/`--disable_git_commit` if passed, or otherwise
by the global preference described in Step 11.

## Settle these decisions (ask if not obvious from the tree)

1. **Host C++ module (`<module>`)** — discovered in Step 2 by grepping
   `mom/cpp/*.hpp`, not assumed from any Fortran module name.
2. **Test target** — `test_mom_<module>`, file
   `test_mom/<module>/test_<module>.cpp`. Append to it if it exists;
   create it (and the CMake block) if not.
3. **Capture-fixture location** — default search: the user's CWD, then
   `$0/capture`. Ask if not found in either.
4. **Field-mapping ground truth** — the captured `.meta` file if a
   fixture was found in #3; otherwise the Fortran shim's CAPTURE-mode
   `rec%add(...)` block, from `$2` if supplied (Step 3 greps for it)
   or from a user paste; otherwise none, in which case Step 7 writes a
   `GTEST_SKIP()` stub rather than fabricate one.

If the user already specified any of these, take their values as-is.

## Procedure

### 1. Validate the existing checkout *(runs only when `--enable_src_validate` is passed; skip otherwise and proceed to Step 2)*
   All test work in this skill is based on the `main` branch — the
   same base `generate_amrex_code` uses.

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

   Then validate the tree — stop on the first failure:
   - `$0/mom/cpp` missing → stop:
     `Error: "<value>"/mom/cpp not found — the AMReX kernel must exist before it can be tested. Run generate_amrex_code first.`
   - `$0/test_mom/common` or `$0/test_mom/test_mom_main.cpp` missing →
     stop:
     `Error: test_mom/ scaffold not found under <work-directory> (common/, test_mom_main.cpp). This skill extends that scaffold; it does not create it from scratch.`

   **Confirm the plan.** Print one paragraph naming: the function name,
   and that Step 2 will now search `mom/cpp/*.hpp` for its kernel
   declaration. Then proceed to Step 2.

### 2. Locate the kernel declaration
   Run `grep -rln "[[:space:]]<lowercased_$1>[[:space:]]*(" $0/mom/cpp/*.hpp`
   (match the symbol inside `namespace MOM`, case-sensitive on the
   lowercased form generate_amrex_code uses).
   - 0 matches → stop:
     `Error: MOM::<lowercased_$1> not found under <work-directory>/mom/cpp. Run generate_amrex_code first to produce the kernel.`
   - >1 matches → list the candidate files and ask the user which
     `<module>` to use.
   - 1 match → record `<module>` = that header's basename (without
     `.hpp`), and capture the full declaration: parameter count, each
     parameter's type/constness (`const Box&`, `Array4<const Real>
     const&`, `Array4<Real> const&`, `Real`, `bool`, `int`, or an
     opaque pointer like `OceanOBC*`), and their order. Print the
     captured signature and pause for confirmation before continuing —
     the header has no parameter names (Doxygen-only comments), so the
     order/type list is all Step 6 has to work with.

   Do not proceed to Step 3 until the signature is captured and
   confirmed.

### 3. Locate the capture fixture, or ground the mapping in Fortran source
   Search for `<lowercased_$1>.bin` / `<lowercased_$1>.meta` in the
   user's CWD, then `$0/capture` (per "Settle these decisions" #3).
   - **Found** → record `ground_truth=meta_file` and the resolved
     path. Proceed to Step 4.
   - **Not found** → ask the user for a different fixture location. If
     still not found, this run has no fixture to run against, but it
     can still ground a *real* assertion (not a stub) in the Fortran
     shim, which is source-of-truth for the exact capture field names
     regardless of whether a fixture has been recorded yet:
     - If `$2` (mom6-directory) was supplied, run
       `grep -rln "^[[:space:]]*subroutine[[:space:]]\+<lowercased_or_original_$1>\b" $2/src $2/config_src`
       to find the shim (the dispatcher, not the `_fortran`-renamed
       original — it has the `select case (mode)` / `rec%add(...)`
       block).
       - 1 match → open it, locate the `case (TIMH_capture)` arm, and
         record every `rec%add("_<name>", <var>)` call in order,
         together with each `<var>`'s declared Fortran type (from the
         subroutine's dummy-argument list) so it can be mapped to a
         `.meta`-equivalent type: `type(RealArray_t)` → `RealArray_t`,
         `type(Box_t)` → `Box_t`, `real` → `real64`, `logical` →
         `logical`, `integer` → `integer`. Also read the matching
         `bind(C)` interface block (or the TIM-side bridge
         header/`.cpp` under `$0/mom/cpp/turbotmp_<module>_bridge.*`)
         as an independent cross-check on argument order and type —
         if the two disagree, stop and ask rather than picking one.
         Record `ground_truth=fortran_source` plus the file path and
         line range.
       - 0 matches → fall through to the paste ask below.
       - >1 matches → list the candidate files and ask the user which
         to use.
     - If `$2` was omitted, or the grep above found nothing, ask the
       user to paste the shim's `case (TIMH_capture)` block (the
       `rec%add(...)` calls) inline. If pasted, record
       `ground_truth=user_paste`.
     - If neither yields anything, record `ground_truth=none` — Step 7
       will write a `GTEST_SKIP()` stub instead of a real assertion,
       and Step 9 (if run) should expect a SKIPPED result, not a
       failure.

   This is read-only: nothing under `$2` is ever modified, and no
   Fortran source is copied verbatim into C++ — only the field *names*
   and *types* feed Step 5's mapping.

### 4. Parse the field list *(skip if `ground_truth=none`)*
   - If `ground_truth=meta_file`: read `<capture-path>.meta`, each
     line `<name> <type> <byte-offset>`, `<type>` one of `Box_t`,
     `RealArray_t`, `real64`, `logical`, `integer` (see
     `test_mom/common/captured_io.hpp`'s header comment for the exact
     `.bin` layout — do not re-derive it, `CapturedFile` already
     implements the reader).
   - If `ground_truth=fortran_source` or `user_paste`: the field list
     is exactly the `rec%add(...)` calls captured in Step 3, already
     paired with their inferred `.meta`-equivalent types — no `.bin`
     to parse, since none was captured.

   Print the field list either way.

### 5. Map kernel parameters ↔ captured fields *(skip if `ground_truth=none`)*
   For each parameter captured in Step 2, in declaration order, find
   its counterpart in the Step 4 field list:
   - `const Box&` → a `Box_t` entry, commonly named `_bx`, `_bxH`, or
     `_bxC` — use whichever name the field list actually contains,
     never assume one over the other.
   - `Array4<const Real> const&` (pure input) → a `RealArray_t` entry
     named `_<argname>` — accessed via `CapturedFile::fab_device()`.
   - `Array4<Real> const&` (in/out) → **two** `RealArray_t` entries,
     `_<argname>_before` (the pre-call state, via `fab_device()`) and
     `_<argname>_after` (the expected post-call state, via
     `fab_host()`, used only for the assertion — never passed into the
     kernel call).
   - `Real` / `bool` / `int` scalar → a `real64` / `logical` /
     `integer` entry named `_<argname>`, via `CapturedFile::real64()` /
     `logical()` / `integer()`.
   - An opaque pointer parameter (e.g. `OceanOBC*`) has no captured
     counterpart — pass `nullptr` and note it in the eventual output
     ("OBC-inactive only", matching the existing
     `PPM_reconstruction_x`/`_y` tests).

   Print the derived mapping as a table (parameter → captured
   field(s) → accessor) and pause for user confirmation. If the field
   list's entry count, or any candidate field's type, doesn't line up
   cleanly with the captured signature, **stop and ask** — never guess
   a mapping silently, regardless of whether the field list came from
   a `.meta` file or from Fortran source.

### 6. Determine the target test file and CMake wiring
   - If `$0/test_mom/<module>/test_<module>.cpp` exists **and**
     `$0/test_mom/CMakeLists.txt` already has an
     `add_executable(test_mom_<module> ...)` block referencing it →
     this run only appends a new `TEST(...)` case to that file (Step
     7). No CMake changes.
   - Otherwise → this run creates
     `test_mom/<module>/test_<module>.cpp` (containing only the new
     `TEST(...)` case) and appends a new block to
     `test_mom/CMakeLists.txt`, copied structurally from the existing
     `test_mom_continuity_ppm` block: `add_executable(test_mom_<module>
     test_mom_main.cpp common/captured_io.cpp
     common/amrex_assertions.cpp <module>/test_<module>.cpp)`, the
     `AMReX_GPU_BACKEND STREQUAL "CUDA"` source-properties guard over
     the same file list, `target_include_directories(...
     PRIVATE common)`, `target_link_libraries(... PRIVATE <lib>
     GTest::gtest)`, and `gtest_discover_tests(... EXTRA_ARGS
     --data-dir=${TEST_MOM_DATA_DIR} DISCOVERY_MODE PRE_TEST)`.

     Determine `<lib>` by running
     `grep -B5 "<module>\.cpp" $0/mom/cpp/CMakeLists.txt` and reading
     off the enclosing `add_library(<lib> ...)` — do not assume
     `<lib>` equals `<module>` (e.g. `mom_continuity_ppm.cpp` is built
     into the `mom_continuity` library, not a `mom_continuity_ppm`
     one).

### 7. Write the test case
   Test name: `TEST(<PascalCase(function-name)>, MatchesFortranCapture)`
   — title-case each `_`-separated segment of `$1` and concatenate
   (e.g. `PPM_reconstruction_x` → `PpmReconstructionX`), matching the
   existing tests exactly. One `TEST(...)` per kernel; never merge
   kernels into one case, never touch another kernel's existing
   `TEST(...)` in a shared file.

   **If a mapping was derived (Step 3/4/5 ran, `ground_truth` is
   `meta_file`, `fortran_source`, or `user_paste`):**
   ```cpp
   TEST(<Name>, MatchesFortranCapture) {
       test_mom::CapturedFile captured(test_mom::data_dir / "<lowercased_$1>");

       const auto   bx   = captured.box("<box field>");
       const auto   in1  = captured.fab_device("<input field>");
       auto         out1 = captured.fab_device("<inout field>_before");
       const auto   out1_after = captured.fab_host("<inout field>_after");
       const double scalar1 = captured.real64("<scalar field>");
       // ... one binding per Step 5 mapping entry, in declaration order

       MOM::<lowercased_$1>(bx, in1.const_array(), out1.array(), /* ... */);
       amrex::Gpu::synchronize();

       expect_arrays_equal(out1_after, to_host_fab(out1), "<argname>");
       // ... one expect_arrays_equal per in/out array
   }
   ```
   This is the same real-assertion shape whether or not a fixture file
   actually exists on disk yet — `CapturedFile`'s constructor throws
   `std::runtime_error` if the `.bin`/`.meta` pair is missing, which
   GoogleTest reports as a normal (informative) test failure. That is
   the intended behaviour when `ground_truth=fortran_source` or
   `user_paste`: this kernel is expected to be exercised, so a missing
   fixture should read as "capture pending," not as a routinely-passing
   skip. If `ground_truth` is not `meta_file`, add a one-line comment
   above the `TEST(...)` naming the source the mapping was grounded in
   (file path + line range from Step 3, or "user paste"), so a future
   reader can re-verify it once a real fixture lands. Keep this to a
   **citation** — where the mapping came from — not a narrative about
   how it was derived or why. Beyond that one line, every other comment
   in the test states what the kernel/fixture physically represents
   (mirroring the Fortran doc comments), never commentary about the
   porting process itself — see `generate_amrex_code`'s lessons.md §7
   #10, which applies here too.

   Include only the headers already present in
   `mom_continuity_ppm/test_mom_continuity_ppm.cpp`
   (`gtest/gtest.h`, `AMReX_FArrayBox.H`, `AMReX_Gpu.H`,
   `amrex_assertions.hpp`, `captured_io.hpp`, `data_dir.hpp`, and
   `<module>.hpp`), plus `using test_mom::expect_arrays_equal;` /
   `using test_mom::to_host_fab;` if the file is newly created.

   **If `ground_truth=none` (no fixture, no Fortran source, no
   paste):**
   ```cpp
   TEST(<Name>, MatchesFortranCapture) {
       GTEST_SKIP() << "no captured <lowercased_$1>.{bin,meta} fixture yet";
   }
   ```
   matching the existing `PpmLimitCw84` stub verbatim in style. Do not
   fabricate expected values as a substitute for a missing fixture —
   and do not use this branch merely because the fixture is absent if
   a field mapping was in fact derived from Fortran source in Step 3.

### 8. Extend `test_mom/common/` only if needed
   Reuse `CapturedFile` and `expect_arrays_equal`/`to_host_fab` as-is.
   Only add a new accessor (e.g. for a captured type Step 4 found that
   `captured_io.hpp` doesn't yet expose) or a new assertion helper if
   Step 5/7 genuinely needs one that doesn't exist — never duplicate
   an existing one, never change an existing helper's behavior.

### 9. Build and run *(runs only when `--enable_build` is passed; skip otherwise and proceed to Step 11)*
   If `--enable_build` was not supplied, skip this entire step and
   Step 10, and note in the final "Output to the user" report that the
   harness was written but not built. Do not build unprompted — a first build
   compiles AMReX from source (and requires MPI to be discoverable,
   since upstream AMReX defaults `AMReX_MPI=ON`), which can take a
   long time and is a real environment dependency the user may not
   want triggered as a side effect of writing a test file.

   Default build dir: `$0/build/test_mom` (matching the header comment
   in `test_mom/CMakeLists.txt`), unless the user specifies otherwise.

   ```
   cmake -S $0/test_mom -B <build-dir> -DTEST_MOM_DATA_DIR=<capture-dir-from-step-3-or-a-placeholder>
   cmake --build <build-dir>
   ctest --test-dir <build-dir> --output-on-failure
   ```

   `TEST_MOM_DATA_DIR` only needs to be a non-empty, existing path to
   satisfy the project's configure-time check — it does not need to
   contain this kernel's fixture if `ground_truth` was
   `fortran_source`/`user_paste`/`none`; pass the real capture
   directory from Step 3 when one exists (other tests in the target
   need it to pass).

   Confirm the new test name appears (`ctest --test-dir <build-dir>
   -N | grep <Name>`) and either **passes** (`ground_truth=meta_file`
   with a real fixture on disk), is reported **SKIPPED**
   (`ground_truth=none`), or **fails with a `captured_io: cannot open
   ...` error** (`ground_truth=fortran_source`/`user_paste` but no
   fixture file yet — expected until one is recorded; report it as
   such, not as a build problem). Surface any other build or test
   failure verbatim; do not patch kernel or `common/` code to force a
   pass — a genuine mismatch against the Fortran capture is a real bug
   to report, not to work around.

### 10. Confirm no other test regressed *(runs only when `--enable_build` is passed; skip otherwise)*
   Run the full `ctest --test-dir <build-dir> --output-on-failure`
   suite (not just the new test) and confirm every previously-passing
   case still passes. A new CMake block or a new `common/` helper must
   not change any existing test's behavior.

### 11. Commit and push *(gated by the global git-commit preference, overridable per run — see below)*
   **Decide whether to run this step:**
   1. If `--disable_git_commit` was passed, skip this entire step and
      report the files that were modified so the user can commit
      manually — this run's explicit override wins.
   2. Else if `--enable_git_commit` was passed, run this step — this
      run's explicit override wins.
   3. Else, read the global preference: run
      `cat ~/.claude/preferences.json 2>/dev/null` and inspect the
      `git_commit_and_push` key.
      - If it is `"auto"`, run this step.
      - If it is `"manual"`, or the file does not exist, or the key is
        absent, or the JSON fails to parse — skip this entire step and
        report the files that were modified so the user can commit
        manually. Absent or unreadable configuration means "manual";
        never push on an unconfigured machine.

   Branch: `claude_<lowercased_$1>_unit_test` based on `main`. Stage
   all newly created and modified files. Commit message names the
   test added, the kernel it exercises, what `ground_truth` it's based
   on, and the `ctest` result if Step 9 ran (or "not built — pass
   --enable_build to compile and run it" if it didn't). Use a HEREDOC
   so the trailer is preserved verbatim:

   ```
   BRANCH="claude_$(echo "$1" | tr '[:upper:]' '[:lower:]')_unit_test"
   git -C "$0" checkout -B "$BRANCH"
   git -C "$0" add -A
   git -C "$0" commit -m "$(cat <<'EOF'
   <one-line summary of the unit test added for $1>

   <1-3 line body: kernel symbol, ground_truth source (fixture path /
   Fortran source path+lines / user paste), ctest result or "not built">

   Co-authored-by: Claude <noreply@anthropic.com>
   EOF
   )"
   git -C "$0" push -u origin "$BRANCH"
   ```

   If the push is rejected because the branch already exists upstream
   with unrelated history, stop and surface to the user; do not
   force-push.

## Versioning marker

Every C++ file this skill creates or modifies gets a `// SKILLS: 0.3.1`
marker line — the same shared version number `generate_amrex_code` and
the Fortran-side `generate_cpp_bridge` use (`!!SKILLS: 0.3.1`), so
`grep -rn "SKILLS:"` finds every touched file across the Fortran and
C++ trees for a given port. Applies to `test_mom/<module>/test_<module>.cpp`
(Step 7) and to any `test_mom/common/*.{hpp,cpp}` file Step 8 extends —
not to `test_mom/CMakeLists.txt`. Place it as its own line immediately
after the file's top-of-file header comment (every test file in this
tree opens with one, e.g. `// Unit tests for ...`), before the first
`#include`; if the file has no header comment, it is the first line.
If a file already has the marker, update the version number in place
rather than adding a second line. Deliberately grep-able and meant to
be stripped later.

## Hard rules

- Never skip the `// SKILLS: 0.3.1` marker on a file this skill
  touches, and never add a second marker line if one already exists.
- Comments state physical meaning (what the kernel/fixture represents)
  or, per Step 7, a one-line ground-truth citation — never a narrative
  about how or why the port/mapping was done (`generate_amrex_code`'s
  lessons.md §7 #10 applies here too).
- Never modify `mom/cpp/` kernel or bridge sources, or any Fortran
  source under `$2` — read-only in both cases, used only to derive the
  kernel signature (Step 2) or ground a field mapping (Step 3).
- Never fabricate expected values. Assertions come only from a
  captured `_after` array or, when grounded in Fortran source, from
  the field names/types `rec%add(...)` actually records; `GTEST_SKIP()`
  is the only acceptable substitute, and only when `ground_truth=none`.
- Never guess a field-name mapping (Step 5). If the field list's
  entries don't line up cleanly with the kernel's parameter list —
  count, order, or type — stop and ask, regardless of whether that
  list came from a `.meta` file or from Fortran source.
- Never build (Steps 9-10) unless `--enable_build` was passed — no
  exceptions, no asking-then-proceeding-anyway. A missing fixture is
  not a reason to skip writing the real assertion (Step 7), but it is
  never a reason to build without being asked.
- Append only. Never remove, reorder, or rewrite an existing
  `TEST(...)` case in a shared module test file, and never restructure
  an existing `add_executable`/`gtest_discover_tests` block beyond
  what Step 6 adds.
- Mirror the CUDA source-properties guard exactly on every new source
  list — every file compiled into a `test_mom_*` target needs it.
- One `TEST(<PascalCaseKernel>, MatchesFortranCapture)` per kernel; do
  not invent a different naming convention.
- Never patch `mom/cpp/` or `test_mom/common/` code just to make a
  failing assertion pass — report the mismatch instead.

## Output to the user on success

Report:

1. Kernel tested (`MOM::<lowercased_$1>`) and its module.
2. The field mapping used and what it was grounded in — the fixture
   path (`meta_file`), the Fortran source file + line range
   (`fortran_source`), the user's paste (`user_paste`), or "deferred —
   GTEST_SKIP stub written, no fixture and no Fortran source available"
   (`none`).
3. Test file created or extended, and any CMake changes.
4. Build result: either "not built — pass --enable_build to compile
   and run it" (Step 9 skipped), or the new test's `ctest` outcome
   (PASS / SKIPPED / expected capture-pending failure) plus
   confirmation the rest of the suite still passes.
5. Branch pushed to `origin`, or the list of files to commit manually.
