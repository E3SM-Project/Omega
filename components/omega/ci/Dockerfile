# Omega CI image: an E3SM container image (github.com/E3SM-Project/containers) plus the
# static Metis/ParMETIS that Omega's FindParmetis.cmake needs and the base image lacks.
#
# A stopgap to demonstrate Omega CI until metis/parmetis are added to the base images
# upstream (as an optional build arg, like WITH_MOAB); then CI can use the base image
# directly and this file can go away.
#
#   docker build -f components/omega/ci/Dockerfile \
#     --build-arg BASE_IMAGE=ghcr.io/e3sm-project/e3sm-image:gnu13-openmpi4 \
#     -t omega-ci components/omega/ci
ARG BASE_IMAGE=ghcr.io/e3sm-project/e3sm-image:gnu13-openmpi4
FROM ${BASE_IMAGE}

# Same specs Polaris deploys on the supported machines (polaris deploy/spack.yaml.j2), so CI
# exercises the same index/real widths. ~shared because FindParmetis.cmake looks for the
# static libraries only.
ARG METIS_SPEC="metis@5.1.0+int64+real64~shared"
ARG PARMETIS_SPEC="parmetis@4.0.3+int64~shared"

# The base image's spack store and profile.d are root's. Build with the image's own compiler
# (COMPILER_SPEC, from profile.d) and MPI (E3SM_MPI, from image.env), as the base image does
# for its I/O libraries. In spack 1.x "%compiler" binds to the spec it follows, so it is
# given for both packages.
USER root
RUN . /etc/profile && . /etc/e3sm/image.env && \
    spack install --fail-fast \
        "${PARMETIS_SPEC} %${COMPILER_SPEC} ^${METIS_SPEC} %${COMPILER_SPEC} ^${E3SM_MPI}" && \
    { echo "export METIS_ROOT=$(spack location -i "${METIS_SPEC}")"; \
      echo "export PARMETIS_ROOT=$(spack location -i "${PARMETIS_SPEC}")"; } \
        > /etc/profile.d/omega-env.sh && \
    cat /etc/profile.d/omega-env.sh && \
    spack gc -y && \
    spack clean -sdfmp

# Back to the base image's runtime user; its ENTRYPOINT (e3sm-env) is inherited.
USER ${E3SM_USER}
