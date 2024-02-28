import os
import pathlib

from skbuild import setup

src_root = os.path.relpath(pathlib.Path(__file__).parent.joinpath("../../../").resolve())
# src_root = pathlib.Path(__file__).parent.joinpath("../../../").resolve()
print(pathlib.Path(__file__).parent.joinpath("../../../").resolve())
print(src_root)
version="7.7.28"

setup(
    name="pyEXP",
    version=version,
    description="Nbody EXPansion Code - exputil Library",
    author="",
    license="GPL-3.0",
    packages=["pyEXP"],
    python_requires=">=3.8",
    install_requires=[
      f"exp-libraries=={version}"    ],
    # cmake_source_dir=str(src_root),
    cmake_minimum_required_version="3.21",
    cmake_languages=("C", "CXX", "Fortran"),
    cmake_args=[
      "-DCMAKE_CXX_STANDARD=17",
      "-DBUILD_COMMON_LIBRARIES=OFF",
      "-DENABLE_NBODY=OFF",
      "-DENABLE_PYEXP=ON",
      "-DBUILD_UTILS=OFF",
      "-DUSE_SUBMODULES=OFF",
      "-DENABLE_PNG=OFF",
      "-DENABLE_USER=OFF",
      "-DENABLE_TESTS=OFF",
      "-DENABLE_XDR=ON",
      "-DCMAKE_BUILD_TYPE=RELEASE"
    ]
    # extras_require={
    #   "cuda": []
    # }
)