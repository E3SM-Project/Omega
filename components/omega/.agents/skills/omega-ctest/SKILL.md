---
name: omega-ctest
description: Build Omega and run its CTest unit tests with Polaris's Omega CTest utility. Use only when the developer asks for this helper and CTests are the focus, instead of or before full PR testing.
---

# Omega CTests with Polaris

Polaris's Omega CTest utility is an optional helper for building Omega and
running its CTests. Use it only when the developer wants it, and only when
the CTests are the focus: on their own, or as a quick check before the
full testing in the `omega-pr-testing` skill. Without it, build and test
by hand as described in `doc/devGuide/QuickStart.md` and
`doc/devGuide/CMakeBuild.md`.

## Running it

The utility is `utils/omega/ctest/omega_ctest.py` in E3SM-Project/polaris,
and its `README.md` there lists every option. Read it before running.

1. Work in a Polaris checkout on a machine that both E3SM and Polaris
   support, with the Polaris load script for the compiler sourced. If
   there is no load script for it, ask the developer before running
   `./deploy.py`.
2. Point the utility at the Omega branch to test with `-o <path>`, or at
   an existing build with `-p <path>`. Add `-d` for a debug build.
3. On a login node, the utility writes a job script, and `-s` submits it.
   Ask the developer before submitting.

When the tests finish, the utility prints a summary for the pull request's
`Testing` comment, and writes it to `ctest_output_for_pr.md` in the build
directory. Report it as the `testing-comments` skill describes.
