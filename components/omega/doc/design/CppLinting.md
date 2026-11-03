(omega-design-cpp-linting)=
# C++ Linting and Static Analysis

**Table of Contents**
1. [Overview](#1-overview)
2. [Requirements](#2-requirements)
3. [Algorithmic Formulation](#3-algorithmic-formulation)
4. [Design](#4-design)
5. [Verification and Testing](#5-verification-and-testing)

## 1 Overview

Omega's C++ code is currently checked only by `clang-format`, which enforces
formatting but finds no bugs. Since Omega's first `pre-commit` configuration,
the plan has been to also run `clang-tidy`, `cppcheck` and
`include-what-you-use`. Those hooks have been disabled since 2023, waiting for
a CMake build they could use.

The obstacle is that clang-based tools must know how each source file is
compiled: its include paths, preprocessor definitions and language standard.
They read this from a compilation database (`compile_commands.json`), which
CMake can write while configuring a build. Omega's standalone CMake build,
however, always creates a CIME case first, which only works on machines that
CIME supports. GitHub-hosted runners and most laptops are not among them, so
the hooks could never get a compilation database where they run.

This design:

- adds `cppcheck` to the checks run on every commit and in CI;
- adds `clang-tidy` to CI, with an optional hook that developers can run
  locally;
- adds a lint-only CMake project that writes a compilation database in seconds
  without CIME, without building anything and without a Fortran compiler;
- does not adopt `include-what-you-use`, and uses a `clang-tidy` check for
  unused includes instead;
- pins all C++ linting tools to consistent, deliberately updated versions.

A prototype of this design, run on the `develop` branch in September 2026,
found real defects that formatting alone could not. Among them were a test
range check that can never fail and several null-pointer checks that log an
error and then use the pointer anyway.

## 2 Requirements

### 2.1 Requirement: Pinned, consistent tool versions

Developers, `pre-commit` and CI must run the same versions of every linting
tool, and the versions must change only on purpose. Different major releases of
`clang-format` format the same code differently, so an unpinned version leads to
changes that each contributor's tools undo. `clang-format` and `clang-tidy`
must come from the same LLVM release.

### 2.2 Requirement: No dependence on CIME or a supported machine

Every check must run on GitHub-hosted runners and on developers' own machines.
Nothing in the linting workflow may create a CIME case or require a machine
configuration from `cime_config`.

### 2.3 Requirement: No build required

Linting must not require building Omega or its dependencies (Kokkos, SCORPIO,
etc.). Preparing for `clang-tidy` must take seconds, not the length of a build.

### 2.4 Requirement: The same view of the code as the real build

The linters must analyze Omega's C++ sources with the same source lists,
include directories and preprocessor definitions that the real build uses. A
new source file or source directory added to Omega's build must be picked up
without a separate edit for linting.

### 2.5 Requirement: Enabled checks are clean on `develop`

A check is enforced only once `develop` has no findings from it, either because
the findings were fixed or because they were suppressed deliberately. A failure
then always means the pull request introduced a problem, never that it touched
a file with old findings.

### 2.6 Requirement: Fast commit-time checks

Checks that run on every commit must take seconds per changed file and need
nothing beyond Omega's development conda environment. Slower checks run in CI.

### 2.7 Requirement: Explicit, local suppressions

When a finding is a false positive or is accepted on purpose, the suppression
must be visible in review. It must name the specific check and be as close as
possible to the code it applies to.

### 2.8 Desired: Coverage of device and single-precision code paths

Code compiled only for GPU backends (`OMEGA_TARGET_DEVICE`, CUDA, HIP, SYCL) or
only with `OMEGA_SINGLE_PRECISION` would ideally be analyzed too. This design
analyzes the Serial, double-precision configuration only. Other configurations
could be added later as extra compilation databases.

## 3 Algorithmic Formulation

This design involves no numerical algorithms.

## 4 Design

### 4.1 Tool selection

The following tools were evaluated on the `develop` branch in September 2026.
Timings are from a 4-core configuration, comparable to a GitHub-hosted runner,
except where noted.

| Tool | Needs compilation database | Cost | Decision |
|------|---------------------------|------|----------|
| `clang-format` | No | < 1 s per file | Keep |
| `cppcheck` | No | about 1 s per file | Adopt, every commit and CI |
| `clang-tidy` | Yes | about 4 min for all of Omega | Adopt, CI (optional locally) |
| `include-what-you-use` | Yes | about 5 s per file | Do not adopt |
| `cpplint`, `lizard` | No | — | Remove from the environment |

**`cppcheck`** analyzes each file on its own and tolerates missing headers, so it
needs only Omega's source directories as include paths. It does need to be told
about Kokkos macros such as `KOKKOS_LAMBDA`, which is done with a small
`cppcheck` configuration file. Almost all of its bug-oriented findings were
either real problems or were caused by `OMEGA::Error::abort()` not being declared
`[[noreturn]]`. Fixing that declaration removed 45 false positives.

**`clang-tidy`** compiles each file with the clang front end and runs checks on
the resulting syntax tree. It finds a wider range of problems but needs a
compilation database (Section 4.3). Enabled checks are drawn mainly from the
bug-finding (`bugprone-*`), performance (`performance-*`) and MPI (`mpi-*`)
groups, together with the check for unused includes. The style-oriented groups
produced thousands of findings, mostly matters of taste. The clang static
analyzer (`clang-analyzer-*`) roughly tripled the run time and found nothing,
so it is not run on every pull request.

**`include-what-you-use`** reports which headers each file should include. For
Kokkos code it asks for Kokkos's internal headers (`impl/...`, `desul/...`),
which would only be fixed by maintaining hand-written mapping files. Each
release also supports exactly one LLVM version, which constrains upgrades. The
most useful part of what it offers, finding unused includes, is available in
`clang-tidy` as `misc-include-cleaner` with missing-include reports turned off.
That option requires `clang-tidy` 21 or newer.

**`cpplint`** enforces Google's style guide, which conflicts with Omega's
`.clang-format` settings. **`lizard`** reports code complexity, which is a
metric rather than a pass or fail check. Neither is used by any hook, so both
are removed from the development environment. So is `include-what-you-use`.

### 4.2 Tool versions

`clang-format`, `clang-tidy` and the clang builtin headers that `clang-tidy`
needs all come from a single LLVM release, initially LLVM 22. In the conda
development environment they are the `clang-format`, `clang-tools` and `clang`
packages, pinned to the same minor release. The `clang` package matters because
`clang-tidy` cannot parse Omega's code without the builtin headers, which
`clang-tools` alone does not provide.

The `pre-commit` hook for `clang-format` uses a mirror that installs a pinned
`clang-format` into the hook's own environment. The formatter used at commit
time is then independent of whatever `clang-format` is first on the path.

`cppcheck` is pinned to a minor release because its findings change between
releases. Its hook and the optional `clang-tidy` hook check that the version on
the path matches the pinned one, and fail with instructions if it does not.

Updating a tool is a deliberate change: all pins move together, the code is
reformatted in a single commit listed in `.git-blame-ignore-revs`, and any new
findings are fixed or suppressed in the same pull request. Moving from
`clang-format` 18.1.8 to 22 changes about 70 lines in 13 files.

### 4.3 A compilation database without CIME

A separate, lint-only CMake project under `components/omega/lint/` produces the
compilation database. It is not a way to build Omega, and Omega's main CMake
build is unchanged. The project:

- sets the variables that Omega's `src/` and `test/` `CMakeLists.txt` files
  expect, then adds those two directories. Source lists, include directories
  and compile definitions therefore come from Omega's own build files
  (Requirement 2.4).
- configures Kokkos from the EKAT submodule with the Serial backend. Kokkos
  generates configuration headers that Omega's code includes, so it is the one
  dependency that is actually configured.
- replaces every other dependency with a header-only placeholder target:
  spdlog, yaml-cpp, cpptrace, GSW-C, GPTL, Pacer and SCORPIO. Where a dependency
  normally generates a header at configure time (cpptrace's `version.hpp`,
  SCORPIO's `pio_config.h`), the project generates it from the same template.
  Configuring SCORPIO itself would require a Fortran compiler.
- takes MPI, NetCDF, PnetCDF, ParMETIS and METIS headers from the conda
  environment.

Configuring takes about 6 seconds, compiles and links nothing, and works with
any C++ compiler that CMake accepts. A small driver script runs the configure
step. It then filters the database to Omega's own files and removes the
duplicate entries from the single-precision library that the tests also use.

Three alternatives were rejected:

- **A lint mode in Omega's main CMake build.** SCORPIO's CMake requires a
  Fortran compiler, and cpptrace downloads libdwarf while configuring. It would
  also add another mode to `OmegaBuild.cmake`, which is already complex.
- **CIME's `linux-generic` machine on CI.** This keeps CIME in the loop (against
  Requirement 2.2), is slow, and depends on CIME behaving on an unsupported
  platform.
- **A hand-written list of compiler flags.** It would drift from the real build
  whenever a source directory or dependency changes (against Requirement 2.4).

When Omega gains a new external dependency, the lint project needs a matching
placeholder. Forgetting one fails loudly, as a CMake error or a `clang-tidy`
"file not found" error, so the problem cannot go unnoticed.

The project needs these submodules: `externals/ekat` (with its Kokkos
submodule), `externals/scorpio` (for headers only) and the submodules under
`components/omega/external/`. Developers who build Omega already have them.

### 4.4 Where each check runs

| Check | On every commit (`pre-commit`) | In CI |
|-------|------------------------------|-------|
| `clang-format` | Yes, changed files | Yes, through `pre-commit` |
| `cppcheck` | Yes, changed files | Yes, through `pre-commit` |
| `clang-tidy` | Optional hook, run on request | Yes, all of Omega |

The optional `clang-tidy` hook is not part of the default `pre-commit` stage.
Developers run it on request, for example with
`pre-commit run clang-tidy --hook-stage manual --files <files>`. It creates the
compilation database the first time it is needed and reuses it afterwards.

In CI, a new job in the Omega pull request workflow initializes only the
submodules listed in Section 4.3, using shallow clones and HTTPS URLs. It then
creates the conda environment, generates the compilation database and runs
`clang-tidy` on every Omega translation unit. It checks the whole of Omega
rather than only the changed files so that a change to a header is checked
through every file that includes it.

Header files are analyzed through the source files that include them, with
`clang-tidy`'s header filter limited to Omega's `src/` and `test/`
directories. Headers are not analyzed on their own because not all of them are
self-contained. Kokkos and other third-party headers are never reported.

### 4.5 Configuration and suppressions

- `components/omega/.clang-tidy` lists the enabled checks and their options,
  such as treating Kokkos Views as cheap to pass by value. It sets the header
  filter and makes enabled checks errors. `clang-tidy` finds this file
  automatically for every file under `components/omega/`.
- The `cppcheck` configuration and suppression files live in
  `components/omega/lint/`. The configuration defines the Kokkos macros.
- A finding is suppressed inline, preferably with the specific check's name:
  `// NOLINT(check-name)` for `clang-tidy` and
  `// cppcheck-suppress checkId` for `cppcheck`. Each suppression gives a brief
  reason on the same or the preceding line (Requirement 2.7).
- A check that is noisy across the whole code base is disabled in the
  configuration rather than suppressed line by line.

The configuration files are the authoritative list of enabled checks. The
developer's guide ({ref}`omega-dev-linting`) explains how to run the tools and
how to handle findings.

### 4.6 Adoption in stages

Each stage leaves `develop` clean for every check that is enabled (Requirement
2.5):

1. **Tool baseline.** Pin LLVM 22, switch `clang-format` to the pinned hook,
   remove unused tools from the environment, and reformat.
2. **Fixes for issues found during evaluation.** These include
   `[[noreturn]]` on `OMEGA::Error::abort()`, logging macros that behave as
   single statements, compile-time type selection in global reductions, a
   test range check that can never fail, and null-pointer checks that continue
   after logging an error.
3. **`cppcheck`.** Enable the warning-level checks with their fixes, followed
   closely by the performance checks with theirs.
4. **`clang-tidy`.** Add the lint project, `.clang-tidy`, the CI job and the
   optional hook, starting from a set of checks that is already clean. Further
   checks are enabled as their findings are fixed.

## 5 Verification and Testing

### 5.1 Test: Complete compilation database

On a GitHub-hosted runner, the lint project configures without CIME and without
building anything. Every C++ source file in `components/omega/src` and
`components/omega/test` that Omega's build compiles has an entry in the
database. `clang-tidy` reports no compiler errors (`clang-diagnostic-error`)
for any of them.
  - tests requirements 2.2, 2.3, 2.4

### 5.2 Test: Known defects are caught

When each tool is enabled, a temporary change that reintroduces one of the
defects found during evaluation (for example, a range check using `&&` where
`||` is needed) is confirmed to fail the corresponding check. The change is not
committed.
  - tests requirement 2.5

### 5.3 Test: Clean `develop`

After each stage, running all enabled checks on all Omega files on `develop`
produces no findings.
  - tests requirement 2.5

### 5.4 Test: Pinned versions

The CI jobs print the version of each tool. A hook run with a mismatched
`cppcheck` or `clang-tidy` on the path fails with instructions, and
`clang-format` in `pre-commit` reports the pinned version whichever
`clang-format` is first on the path.
  - tests requirement 2.1

### 5.5 Test: Run time

`cppcheck` takes no more than a few seconds per changed file in `pre-commit`.
The CI `clang-tidy` job finishes well within its time limit on a GitHub-hosted
runner.
  - tests requirement 2.6

### 5.6 Test: No change in model behavior

The source fixes in stage 2 pass the Omega CTest suite and the Polaris
`omega_pr` suite. Any change in results is explained by a defect that was
fixed.
  - tests that the fixes are safe
