#!/bin/bash
# Build Omega standalone and run its CTests inside the Omega CI image (ci/Dockerfile), using
# the image's own CIME machine, e3sm-container. Run as the image's default user under a login
# shell, with the Omega checkout (required submodules included) and the meshes mounted:
#
#   docker run --rm -e OMEGA_BUILD_TYPE=Debug \
#     -v /path/to/Omega:/projects/e3sm/work/Omega \
#     -v /path/to/omega-meshes:/projects/e3sm/data/omega-meshes \
#     -v /path/to/build:/home/e3sm/omega-build \
#     omega-ci bash -l /projects/e3sm/work/Omega/components/omega/ci/omega-ctest.sh
#
# The checkout and build mounts must be writable by the e3sm user (uid 1000).
set -euo pipefail

: "${OMEGA_BUILD_TYPE:=Release}"
: "${OMEGA_SRC:=/projects/e3sm/work/Omega}"
: "${OMEGA_BUILD_DIR:=$HOME/omega-build}"
: "${OMEGA_MESH_DIR:=/projects/e3sm/data/omega-meshes}"
: "${CTEST_ARGS:=--timeout 900}"
step() { echo "== $*"; }

step "CIME machine"
# The base images ship the e3sm-container machine without enabling it; opt in
# shellcheck source=/dev/null
. /projects/e3sm/cime/enable.sh
test "${CIME_MACHINE:-}" = e3sm-container ||
    { echo "CIME_MACHINE is '${CIME_MACHINE:-}', not e3sm-container" >&2; exit 1; }
# The image's default compiler: the first one its machine lists
compiler=$(sed -n 's:.*<COMPILERS>\([^,<]*\).*:\1:p' "$HOME/.cime/config_machines.xml")
: "${METIS_ROOT:?not set; is this the Omega CI image (ci/Dockerfile)?}"
: "${PARMETIS_ROOT:?not set; is this the Omega CI image (ci/Dockerfile)?}"
echo "machine $CIME_MACHINE, compiler $compiler, build type $OMEGA_BUILD_TYPE"
echo "Omega $(git -C "$OMEGA_SRC" rev-parse HEAD 2>/dev/null || echo unknown)"

step "configure"
# The build mount may be reused; start from scratch as the dev guide recommends
mkdir -p "$OMEGA_BUILD_DIR"
find "$OMEGA_BUILD_DIR" -mindepth 1 -delete
cmake \
    -DOMEGA_BUILD_TYPE="$OMEGA_BUILD_TYPE" \
    -DOMEGA_CIME_COMPILER="$compiler" \
    -DOMEGA_CIME_MACHINE="$CIME_MACHINE" \
    -DOMEGA_METIS_ROOT="$METIS_ROOT" \
    -DOMEGA_PARMETIS_ROOT="$PARMETIS_ROOT" \
    -DOMEGA_BUILD_TEST=ON \
    -Wno-dev \
    -S "$OMEGA_SRC/components/omega" \
    -B "$OMEGA_BUILD_DIR"
cd "$OMEGA_BUILD_DIR"
# OmegaBuild.cmake ignores create_newcase/case.setup failures; this is where they show
test -x omega_build.sh || { echo "CMake did not write omega_build.sh" >&2; exit 1; }

step "build"
./omega_build.sh

step "meshes"
for mesh in OmegaMesh.nc OmegaSphereMesh.nc OmegaPlanarMesh.nc; do
    test -f "$OMEGA_MESH_DIR/$mesh" ||
        { echo "missing $OMEGA_MESH_DIR/$mesh; run ci/get-meshes.sh" >&2; exit 1; }
    ln -sf "$OMEGA_MESH_DIR/$mesh" "test/$mesh"
done

step "ctest $CTEST_ARGS"
# shellcheck disable=SC2086
./omega_ctest.sh $CTEST_ARGS
