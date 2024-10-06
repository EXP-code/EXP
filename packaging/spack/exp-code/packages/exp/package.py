# Copyright 2013-2024 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

# ----------------------------------------------------------------------------
# If you submit this package back to Spack as a pull request,
# please first remove this boilerplate and all FIXME comments.
#
# This is a template package file for Spack.  We've put "FIXME"
# next to all the things you'll want to change. Once you've handled
# them, you can save this file and test your package like this:
#
#     spack install libexp
#
# You can edit this file again by typing:
#
#     spack edit libexp
#
# See the Spack documentation for more information on packaging.
# ----------------------------------------------------------------------------

from spack.package import *


class Exp(CMakePackage):
    """FIXME: Put a proper description of your package here."""

    # FIXME: Add a proper url for your package's homepage here.
    homepage = "https://www.example.com"
    url = "https://github.com/EXP-code/EXP/archive/refs/tags/v7.7.99.tar.gz"
    git = "https://github.com/georgiastuart/EXP.git"

    maintainers("georgiastuart")
    license("GPL-3.0-or-later", checked_by="georgiastuart")

    version("develop", branch="build-tools")
    version("7.7.99", sha256="c28394fef7de19ba1f4db771d05a9a373dd18ae3ca1b757a9ea424eee78b460d")

    depends_on("highfive@develop +mpi", type="build") # Needs HighFive 3
    depends_on("eigen@3:", type="build")
    depends_on("boost +math")
    depends_on("mpi")
    depends_on("fftw@3:")
    depends_on("hdf5@1.8.20:1.9,1.10.2: +cxx +hl")
    depends_on("yaml-cpp@0.8:,develop")
    depends_on("zlib")

    depends_on("c", type="build")
    depends_on("cxx", type="build")
    depends_on("fortran", type="build")

    variant("python", default=False, description="Build the pyEXP Python module.")
    depends_on("python@3.8:", when="+python")
    depends_on("py-pybind11", when="+python", type="build")
    
    variant("xdr", default=False, description="Enable RPC/XDR support for Tipsy standard.")
    depends_on("libtirpc", when="+xdr")

    variant("nbody", default=True, description="Build the nbody library and executable.")
    variant("utils", default=True, description="Build utility executables.")

    variant("slurm", default=False, description="Enable SLURM checkpointing support.")
    depends_on("slurm", when="+slurm")

    variant("vtk", default=False, description="Configure VTK.")
    depends_on("vtk", when="+vtk")

    variant("docs", default=False, description="Build the documentation.")
    depends_on("doxygen", when="+docs")

    variant("tests", default=False, description="Build the test suite.")

    def cmake_bool(self, variant):
        return "ON" if self.spec.satisfies(f"+{variant}") else "OFF"
    
    def cmake_args(self):
        args = [
            "-DCMAKE_CXX_STANDARD=17",
            "-DINSTALL_HEADERS=ON",
            "-DBUILD_COMMON_LIBRARIES=ON",
            "-DINSTALL_CMAKE_FIND=ON",
            f"-DENABLE_NBODY={self.cmake_bool('nbody')}",
            f"-DENABLE_PYEXP={self.cmake_bool('python')}",
            f"-DBUILD_UTILS={self.cmake_bool('utils')}",
            f"-DENABLE_SLURM={self.cmake_bool('slurm')}",
            f"-DBUILD_DOCS={self.cmake_bool('docs')}",
            f"-DUSE_SUBMODULES={'ON' if self.spec.satisfies('@:7.7.99') else 'OFF'}",
            "-DENABLE_PNG=OFF", # png++ is only compatible with libpng 1.2.X - defunct
            "-DENABLE_USER=OFF",
            f"-DENABLE_TESTS={self.cmake_bool('tests')}",
            f"-DENABLE_XDR={self.cmake_bool('xdr')}",
            "-DENABLE_DSMC=OFF",
            "-Wno-dev"
        ]
        return args
