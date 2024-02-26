import pathlib

from skbuild import setup

src_root = pathlib.Path(__file__).parent.joinpath("../../../").relative_to(pathlib.Path.cwd())
print(src_root.resolve())

setup(
    name="pyEXP",
    version="7.7.27",
    description="Nbody EXPansion Code",
    author="",
    license="GPL-3.0",
    packages=['pyEXP'],
    python_requires=">=3.8",
    package_dir={"": str(src_root)},
    cmake_source_dir=str(src_root),
    cmake_minimum_required_version="3.21",
    cmake_languages=("C", "CXX", "Fortran"),
    cmake_args=[
      "-DCMAKE_CXX_STANDARD=17",
      "-DENABLE_NBODY=ON",
      "-DENABLE_PYEXP=ON",
      "-DUSE_SUBMODULES=OFF",
      "-DENABLE_PNG=ON",
      "-DCMAKE_BUILD_TYPE=RELEASE"
    ]
    # extras_require={
    #   "cuda": []
    # }
)