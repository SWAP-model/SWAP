#!/usr/bin/env bash
# Provision pFUnit for the local unit-test build.
#
# SWAP's unit suite (tests/unit/) builds against pFUnit, discovered by
# tests/unit/meson.build at the path PFUNIT_ROOT points to
# (tests/pFUnit/build/install_gfortran/PFUNIT-4.15 by default — see the
# [activation] env in pixi.toml). pFUnit is not on conda-forge or PyPI, so it
# must be built from source. This script does that idempotently.
#
# pFUnit is pinned to v4.15.0 (the version the install path encodes). The
# tests/pFUnit/ checkout is build-time-provisioned and git-ignored — it is no
# longer a tracked gitlink.
set -euo pipefail

PFUNIT_TAG="v4.15.0"
PFUNIT_REPO="https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git"

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
src_dir="${repo_root}/tests/pFUnit"
build_dir="${src_dir}/build"
prefix="${build_dir}/install_gfortran"
install_dir="${prefix}/PFUNIT-4.15"

if [ -x "${install_dir}/bin/funitproc" ] && [ -f "${install_dir}/include/driver.F90" ]; then
    echo "pFUnit already provisioned at ${install_dir} — nothing to do."
    exit 0
fi

echo "Provisioning pFUnit ${PFUNIT_TAG} -> ${install_dir}"

# Source: clone the pinned tag (with GFE sub-deps) if not already present.
if [ ! -f "${src_dir}/CMakeLists.txt" ]; then
    echo "Cloning pFUnit ${PFUNIT_TAG} (recursive) into ${src_dir} ..."
    rm -rf "${src_dir}"
    git clone --depth 1 --branch "${PFUNIT_TAG}" --recurse-submodules --shallow-submodules \
        "${PFUNIT_REPO}" "${src_dir}"
fi

# Configure + build + install. The suite is serial: no MPI / OpenMP / ESMF.
FC="${FC:-gfortran}" cmake -S "${src_dir}" -B "${build_dir}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${prefix}" \
    -DSKIP_MPI=YES \
    -DSKIP_OPENMP=YES \
    -DSKIP_FHAMCREST=YES

cmake --build "${build_dir}" --parallel
cmake --install "${build_dir}"

if [ ! -x "${install_dir}/bin/funitproc" ]; then
    echo "ERROR: pFUnit install did not produce ${install_dir}/bin/funitproc" >&2
    exit 1
fi
echo "pFUnit ${PFUNIT_TAG} installed at ${install_dir}"
