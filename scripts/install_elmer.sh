#!/usr/bin/env bash
set -euo pipefail

source_dir="${PIXI_PROJECT_ROOT}/.elmer-source"
build_dir="${source_dir}/build"
install_dir="${PIXI_PROJECT_ROOT}/.elmer"

if [[ ! -d "${source_dir}" ]]; then
    git clone --branch release-26.2.1 --depth 1 \
        https://github.com/ElmerCSC/elmerfem.git "${source_dir}"
    git -C "${source_dir}" checkout a19504ac53ec222e3355e182b08f2ff280c2203a
fi

cmake -S "${source_dir}" -B "${build_dir}" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER="${CONDA_PREFIX}/bin/x86_64-conda-linux-gnu-cc" \
    -DCMAKE_CXX_COMPILER="${CONDA_PREFIX}/bin/x86_64-conda-linux-gnu-c++" \
    -DCMAKE_Fortran_COMPILER="${CONDA_PREFIX}/bin/x86_64-conda-linux-gnu-gfortran" \
    -DINCLUDE_ELMER_MODULES="ElasticSolve;ResultOutputSolve" \
    -DCMAKE_INSTALL_PREFIX="${install_dir}" \
    -DWITH_ELMERGUI=FALSE \
    -DWITH_MPI=FALSE
cmake --build "${build_dir}" --parallel "${ELMER_BUILD_JOBS:-8}"
cmake --install "${build_dir}"
