---
name: omega-pr-testing
description: Run the testing the Omega pull request template asks for (CTests, the Polaris omega_pr suite and coupled CIME tests) with Polaris's Omega PR testing utility. Use only when the developer asks for testing to be run or coordinated this way.
---

# Omega PR testing with Polaris

Polaris's Omega PR testing utility is one way to run the testing the pull
request template asks for. It is not the required way. Developers can also
run the tests by hand, following `doc/devGuide/Testing.md`, or any other
way that covers the template's rows. Use the utility only when the
developer asks for it.

If the developer wants only the CTests, use the `omega-ctest` skill
instead.

## Where the instructions are

The utility is in `utils/omega/pr_testing` of E3SM-Project/polaris, and
its `AGENTS.md` there is the procedure. Read that file in full, in an
up-to-date checkout of Polaris `main`, and follow it. It changes with the
utility, so do not work from a copy or from memory.

In outline, it pins the PR and a baseline, builds both with Polaris on
each machine, runs the CTests, `omega_pr` and coupled tests, and writes a
`Testing` comment per machine. Pushing the test branch, submitting jobs
and posting comments each need the developer's permission first.

## The comments

The utility writes the `Testing` comments. Add your own findings through
its notes option rather than editing its output, and follow the
`testing-comments` skill for those notes.
