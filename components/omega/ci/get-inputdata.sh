#!/bin/bash
# Download the CIME inputdata a grid and compset need inside an E3SM container
# image, without building. The grid and compset are given as GRID.COMPSET, as in
# a CIME test name. Run as the image's default user in a login shell, with the
# Omega checkout, inputdata and scratch mounted:
#
#   docker run --rm \
#     -v /path/to/Omega:/projects/e3sm/work/Omega \
#     -v /path/to/inputdata:/projects/e3sm/data/inputdata \
#     -v /path/to/scratch:/home/e3sm/e3sm_scratch \
#     ghcr.io/e3sm-project/e3sm-image:gnu13-openmpi4 \
#     bash -l /projects/e3sm/work/Omega/components/omega/ci/get-inputdata.sh \
#     T62_oQU240.GOMEGA-NYF
#
# All three mounts must be writable by the e3sm user (uid 1000).
set -euo pipefail

src_dir=${OMEGA_SRC:-/projects/e3sm/work/Omega}
grid_compset=${1:-}
case_dir=$HOME/e3sm_scratch/get-inputdata

error() {
    echo "ERROR: $*" >&2
    exit 1
}

grid=${grid_compset%%.*}
compset=${grid_compset#*.}
if [[ -z "$grid" || "$compset" == "$grid_compset" ]]; then
    error "usage: $0 GRID.COMPSET"
fi

echo "== CIME machine"
# The images ship the e3sm-container machine, but it has to be enabled
source /projects/e3sm/cime/enable.sh
if [[ "${CIME_MACHINE:-}" != e3sm-container ]]; then
    error "CIME_MACHINE is '${CIME_MACHINE:-}', not e3sm-container"
fi

echo "== create_newcase $grid $compset"
rm -rf "$case_dir"
cd "$src_dir/cime/scripts"
./create_newcase --case "$case_dir" --res "$grid" --compset "$compset" \
    --mach e3sm-container

# preview_namelists writes the Buildconf/*.input_data_list files that
# check_input_data reads
cd "$case_dir"
./case.setup
./preview_namelists
./check_input_data --download

# The workflow keys the inputdata cache on this list
cat Buildconf/*.input_data_list | sort > "$HOME/e3sm_scratch/input_data_list"

rm -rf "$case_dir"
