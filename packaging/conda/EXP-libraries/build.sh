#!/bin/bash

REL_PREFIX=$( realpath --relative-to="${PWD}" "${PREFIX" )
sed -i s#cmake_install_dir=\"\"#cmake_install_dir=\"${REL_PREFIX}\"# setup.py
python -m pip install . -vv --no-deps --no-build-isolation
if errorlevel 1 exit 1
