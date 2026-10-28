# Omega CI scripts

Scripts that GitHub Actions runs to build and test Omega in the
[E3SM container images](https://github.com/E3SM-Project/containers), using the
images' CIME machine (`e3sm-container`). Usage is in each script's header.

| File | Purpose | Used by |
|------|---------|---------|
| `get-meshes.sh` | Download the meshes the CTests need (on the runner, outside the container) | [omega-ctest-workflow.yml](../../../.github/workflows/omega-ctest-workflow.yml) |
| `omega-ctest.sh` | Build Omega standalone and run its CTests (in the container) | [omega-ctest-workflow.yml](../../../.github/workflows/omega-ctest-workflow.yml) |
| `omega-cime-test.sh` | Run the CIME system tests given as arguments (in the container) | [omega-cime-test-workflow.yml](../../../.github/workflows/omega-cime-test-workflow.yml) |

Nothing is cached or uploaded as an artifact: meshes and inputdata are
downloaded on every run, and logs of failed runs are printed to the job log.
