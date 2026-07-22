#!/usr/bin/env bash
set -e

R_SRC_DIR="${CONTAINER_WORKSPACE_FOLDER:-/workspaces/seurat_blog}/commands/r-source-4.5.2"

if [ -d "$R_SRC_DIR" ]; then
    echo "=== Setting up R source build in $R_SRC_DIR ==="
    cd "$R_SRC_DIR"

    # 1. Configure R with debug flags and shlib support
    export FFLAGS="-O0 -g"
    export FCFLAGS="-O0 -g"
    export CFLAGS="-O0 -g"
    export CXXFLAGS="-O0 -g"

    ./configure \
        --enable-R-shlib \
        --without-recommended-packages

    # 2. Set up dummy SVN / documentation files to bypass svnonly checks
    mkdir -p .svn doc
    touch doc/FAQ
    echo "Revision: 86000" > SVN-REVISION
    echo "Revision: 86000" > SVN-REVISION-tmp
    rm -f non-tarball

    # 3. Build R
    make GIT="echo Revision: 86000" SVN="echo"
else
    echo "Directory $R_SRC_DIR not found. Skipping R build."
fi