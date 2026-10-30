#!/bin/bash
# Build Omega standalone and run its CTests inside an E3SM container image
# (github.com/E3SM-Project/containers), using the image's CIME machine
# (e3sm-container) and its Metis/ParMETIS. Run as the image's default user in a
# login shell, with the Omega checkout and the meshes mounted:
#
#   docker run --rm -e OMEGA_BUILD_TYPE=Debug \
#     -v /path/to/Omega:/projects/e3sm/work/Omega \
#     -v /path/to/omega-meshes:/projects/e3sm/data/omega-meshes \
#     -v /path/to/build:/home/e3sm/omega-build \
#     ghcr.io/e3sm-project/e3sm-image:gnu13-openmpi4 \
#     bash -l /projects/e3sm/work/Omega/components/omega/ci/omega-ctest.sh
#
# The checkout and build mounts must be writable by the e3sm user (uid 1000).
set -euo pipefail

# Settings, each overridable from the environment
build_type=${OMEGA_BUILD_TYPE:-Release}
src_dir=${OMEGA_SRC:-/projects/e3sm/work/Omega}
build_dir=${OMEGA_BUILD_DIR:-$HOME/omega-build}
mesh_dir=${OMEGA_MESH_DIR:-/projects/e3sm/data/omega-meshes}
read -ra ctest_args <<<"${CTEST_ARGS:---timeout 900}"

error() {
    echo "ERROR: $*" >&2
    exit 1
}

echo "== CIME machine"
# The images ship the e3sm-container machine, but it has to be enabled
# shellcheck source=/dev/null
source /projects/e3sm/cime/enable.sh
if [[ "${CIME_MACHINE:-}" != e3sm-container ]]; then
    error "CIME_MACHINE is '${CIME_MACHINE:-}', not e3sm-container"
fi
if [[ -z "${METIS_ROOT:-}" || -z "${PARMETIS_ROOT:-}" ]]; then
    error "METIS_ROOT or PARMETIS_ROOT not set; was the image built with WITH_PARMETIS=yes?"
fi

# The image's compiler is the first one listed for the machine
compiler=$(sed -n 's:.*<COMPILERS>\([^,<]*\).*:\1:p' "$HOME/.cime/config_machines.xml")

echo "machine $CIME_MACHINE, compiler $compiler, build type $build_type"
echo "Omega $(git -C "$src_dir" rev-parse HEAD 2>/dev/null || echo unknown)"

echo "== configure"
# The build directory may be a mount point, so empty it rather than delete it
mkdir -p "$build_dir"
find "$build_dir" -mindepth 1 -delete

cmake \
    -DOMEGA_BUILD_TYPE="$build_type" \
    -DOMEGA_CIME_COMPILER="$compiler" \
    -DOMEGA_CIME_MACHINE="$CIME_MACHINE" \
    -DOMEGA_METIS_ROOT="$METIS_ROOT" \
    -DOMEGA_PARMETIS_ROOT="$PARMETIS_ROOT" \
    -DOMEGA_BUILD_TEST=ON \
    -Wno-dev \
    -S "$src_dir/components/omega" \
    -B "$build_dir"

cd "$build_dir"

# OmegaBuild.cmake carries on when CIME's create_newcase or case.setup fail, so
# a missing build script is the first sign of it
if [[ ! -x omega_build.sh ]]; then
    error "CMake did not write omega_build.sh; see $build_dir/e3smcase"
fi

echo "== build"
./omega_build.sh

echo "== meshes"
for mesh in OmegaMesh.nc OmegaSphereMesh.nc OmegaPlanarMesh.nc; do
    if [[ ! -f "$mesh_dir/$mesh" ]]; then
        error "missing $mesh_dir/$mesh; run get-meshes.sh"
    fi
    ln -sf "$mesh_dir/$mesh" "test/$mesh"
done

echo "== ctest ${ctest_args[*]}"
./omega_ctest.sh "${ctest_args[@]}"
