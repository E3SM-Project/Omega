#!/bin/bash
# Download the meshes Omega's CTests read into a directory (default: ./omega-meshes),
# skipping any already there. The list matches the Quick Start guide
# (doc/devGuide/QuickStart.md, "Getting test meshes").
#
#   components/omega/ci/get-meshes.sh [mesh_dir]
set -euo pipefail

mesh_dir=${1:-omega-meshes}
base=https://web.lcrc.anl.gov/public/e3sm/polaris/ocean/omega_ctest

meshes=(
    "OmegaMesh.nc ${base}/ocean.QU.240km.omega_vars.260807.nc"
    "OmegaSphereMesh.nc ${base}/cosine_bell_icos480.omega_vars.260911.nc"
    "OmegaPlanarMesh.nc ${base}/PlanarPeriodic48x48.omega_vars.260825.nc"
)

mkdir -p "$mesh_dir"
for entry in "${meshes[@]}"; do
    read -r name url <<< "$entry"
    [ -f "$mesh_dir/$name" ] && continue
    # Into a .part file first, so an interrupted download never looks complete
    curl -fsSL --retry 4 --retry-delay 5 -o "$mesh_dir/$name.part" "$url"
    mv "$mesh_dir/$name.part" "$mesh_dir/$name"
    echo "fetched $name"
done
ls -l "$mesh_dir"
