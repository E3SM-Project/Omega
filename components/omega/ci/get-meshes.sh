#!/bin/bash
# Download the meshes Omega's CTests read into a directory, skipping any already
# there. The list matches doc/devGuide/QuickStart.md ("Getting test meshes").
#
#   components/omega/ci/get-meshes.sh [mesh_dir]   (default: ./omega-meshes)
set -euo pipefail

mesh_dir=${1:-omega-meshes}
base_url=https://web.lcrc.anl.gov/public/e3sm/polaris/ocean/omega_ctest

# Name the tests expect, then the file on the server
meshes=(
    "OmegaMesh.nc        ocean.QU.240km.omega_vars.260807.nc"
    "OmegaSphereMesh.nc  cosine_bell_icos480.omega_vars.260911.nc"
    "OmegaPlanarMesh.nc  PlanarPeriodic48x48.omega_vars.260825.nc"
)

mkdir -p "$mesh_dir"
for entry in "${meshes[@]}"; do
    read -r name file <<< "$entry"
    if [[ -f "$mesh_dir/$name" ]]; then
        echo "have $name"
        continue
    fi

    # Download to a .part file first, so an interrupted download is not mistaken
    # for a complete one on the next run
    curl -fsSL --retry 4 --retry-delay 5 -o "$mesh_dir/$name.part" "$base_url/$file"
    mv "$mesh_dir/$name.part" "$mesh_dir/$name"
    echo "fetched $name"
done
ls -l "$mesh_dir"
