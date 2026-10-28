#!/bin/bash
# Run CIME system tests of Omega inside an E3SM container image, using the
# image's CIME machine (e3sm-container). The tests are given as arguments, as
# full CIME test names without machine or compiler. Run as the image's default
# user in a login shell, with the Omega checkout, inputdata and scratch mounted:
#
#   docker run --rm \
#     -v /path/to/Omega:/projects/e3sm/work/Omega \
#     -v /path/to/inputdata:/projects/e3sm/data/inputdata \
#     -v /path/to/scratch:/home/e3sm/e3sm_scratch \
#     ghcr.io/e3sm-project/e3sm-image:gnu13-openmpi4 \
#     bash -l /projects/e3sm/work/Omega/components/omega/ci/omega-cime-test.sh \
#     ERS_D_Ln5.T62_oQU240.GOMEGA-NYF
#
# Missing inputdata is downloaded into the inputdata mount. Test cases end up in
# ~/e3sm_scratch/<test>.omega-ci. All three mounts must be writable by the e3sm
# user (uid 1000).
set -euo pipefail

src_dir=${OMEGA_SRC:-/projects/e3sm/work/Omega}
tests=("$@")
test_id=omega-ci
scratch=$HOME/e3sm_scratch

error() {
    echo "ERROR: $*" >&2
    exit 1
}

if [[ ${#tests[@]} -eq 0 ]]; then
    error "no tests given; usage: $0 TEST [TEST ...]"
fi

echo "== CIME machine"
# The images ship the e3sm-container machine, but it has to be enabled
# shellcheck source=/dev/null
source /projects/e3sm/cime/enable.sh
if [[ "${CIME_MACHINE:-}" != e3sm-container ]]; then
    error "CIME_MACHINE is '${CIME_MACHINE:-}', not e3sm-container"
fi

echo "machine $CIME_MACHINE"
echo "Omega $(git -C "$src_dir" rev-parse HEAD 2>/dev/null || echo unknown)"

echo "== create_test ${tests[*]}"
mkdir -p "$scratch"
cd "$src_dir/cime/scripts"
# create_test won't start a run needing more ranks than its pool, which defaults
# to the core count plus 25%. The machine oversubscribes, so let 8-rank runs
# start on a 4-core runner.
cores=$(nproc)
rc=0
./create_test "${tests[@]}" --wait --test-id "$test_id" \
    --parallel-jobs "$cores" --proc-pool "$((cores > 8 ? cores : 8))" || rc=$?

shopt -s nullglob
for status in "$scratch"/*."$test_id"/TestStatus; do
    echo "-- $status"
    cat "$status"
done

exit "$rc"
