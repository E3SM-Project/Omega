<!--
Please add a description of what is accomplished in the PR here at the top.
It becomes the message of the merge commit, which ends up in E3SM's master,
so follow E3SM's commit message rules:
https://e3sm.atlassian.net/wiki/spaces/DOC/pages/16253965/Commit+and+PR+description+template
In short: a few sentences in the imperative ("Add ...", "Fix ..."), wrapped
at about 72 characters, saying what changed and why. Refer to issues and PRs
with the repository (E3SM-Project/Omega#123, not #123).

End the description (above any "---" line, see below) with one of these
tags, and add the matching label:
  [BFB]      answers are bit-for-bit
  [non-BFB]  answers change at roundoff level
  [NCC]      larger than roundoff changes, not believed climate changing
  [CC]       climate changing
Add [NML] as well (and its label) if the PR adds or removes a config option
or changes a default.

Put issue references (such as "Fixes E3SM-Project/Omega#123"), links,
figures and tables below a line containing only "---".
-->

<!--
Below are a few things we ask you or your reviewers to kindly check.
***Remove checks that are not relevant by deleting the line(s) below.***
-->
Checklist
* [ ] The description, above any `---` line, ends with a `[BFB]`, `[non-BFB]`, `[NCC]` or `[CC]` tag (plus `[NML]` if needed), with matching labels. If this tag is inferred without testing, append with ` inferred`.
* [ ] Documentation:
  * [ ] Design document has been generated and [added to the docs](https://github.com/E3SM-Project/Omega/tree/develop/components/omega/doc/design)
  * [ ] [User's Guide](https://github.com/E3SM-Project/Omega/tree/develop/components/omega/doc/userGuide) has been updated
  * [ ] [Developer's Guide](https://github.com/E3SM-Project/Omega/tree/develop/components/omega/doc/devGuide) has been updated
  * [ ] Documentation has been [built locally](https://docs.e3sm.org/Omega/Omega/devGuide/BuildDocs.html) and changes look as expected
* Linting
  * [ ] [Pre-commit](https://docs.e3sm.org/Omega/Omega/devGuide/Linting.html) run locally
* Building
  * [ ] CMake build does not produce any new warnings from changes in this PR
* [ ] Testing[^1]
  <!--
  Remove rows for any machine/compiler combinations you do not plan to run (or ask someone else to run).
  PRs are expected to pass tests on at least one CPU and at least one GPU platform for both standalone
  tests (CTests and omega_pr) and coupled CIME tests.
  -->

  **aurora, intel, mpich**
  - [ ] CTests Pass
  - [ ] Polaris omega_pr Pass
  - [ ] Coupled CIME tests Pass

  **aurora, intelgpu, mpich**
  - [ ] CTests Pass
  - [ ] Polaris omega_pr Pass
  - [ ] Coupled CIME tests Pass

  **chrysalis, intel, openmpi**
  - [ ] CTests Pass
  - [ ] Polaris omega_pr Pass
  - [ ] Coupled CIME tests Pass

  **frontier, craygnu, mpich**
  - [ ] CTests Pass
  - [ ] Polaris omega_pr Pass
  - [ ] Coupled CIME tests Pass

  **frontier, craygnu-mphipcc, mpich**
  - [ ] CTests Pass
  - [ ] Polaris omega_pr Pass
  - [ ] Coupled CIME tests Pass

  **pm-cpu, gnu, mpich**
  - [ ] CTests Pass
  - [ ] Polaris omega_pr Pass
  - [ ] Coupled CIME tests Pass

  **pm-gpu, gnugpu, mpich**
  - [ ] CTests Pass
  - [ ] Polaris omega_pr Pass
  - [ ] Coupled CIME tests Pass

[^1]: Remove rows for any configurations you do not run (or ask someone else to run).
    At a minimum, PRs are expected to pass tests on at least one CPU and at least one GPU platform for both standalone tests (CTests and `omega_pr`) and coupled CIME tests.

* [ ] Provide relevant details in a comment to the PR titled `Testing` with the following:
  * Which machines [CTest unit tests](https://docs.e3sm.org/Omega/Omega/devGuide/QuickStart.html#running-ctests)
        have been run on and indicate that are all passing.
  * The [Polaris omega_pr test suite](https://docs.e3sm.org/Omega/Omega/devGuide/Testing.html)
     has passed, using the Polaris `e3sm_submodules/Omega` baseline
  * Document machine(s), compiler(s), and the build path(s) used for `-p` for both the baseline (Polaris `e3sm_submodules/Omega`) and the PR build
  * Which [coupled CIME tests](https://docs.e3sm.org/Omega/Omega/devGuide/Testing.html#coupled-testing-with-cime)
     (such as `ERS_D` or `PEM_D` tests) have been run and indicate that they are all passing
  * Indicate "All tests passed" or document failing tests
  * Document testing used to verify the changes including any tests that are added/modified/impacted.
* [ ] Performance related PRs: Please include a relevant PACE experiment link documenting performance before and after.
* [ ] New tests:
  * [ ] CTest unit tests for new features have been added per the approved design.
  * [ ] Polaris tests for new features have been added per the approved design (and included in a test suite)
* Stealth Features
  * [ ] If any stealth features are included in the PR, please confirm that they have been documented.

<!--
Please note any issues this fixes using closing keywords: https://help.github.com/articles/closing-issues-using-keywords
-->


