from skbuild import setup

version="7.7.28"

# Setup for the CopyAndSdist subclass

import os
import pathlib
import shutil

from setuptools.command.sdist import sdist

class CopyAndSdist(sdist):
  file_list = []
  directory_list = []
  def run(self):
    self.copy_src()
    super().run()
    self.clean_src()

  @classmethod
  def clean_src(cls):
    print(f"Cleaning the packaging directory, {os.getcwd()}...")
    for file in cls.file_list:
      try:
        os.remove(file)
        print(f"Removed file {file}")
      except FileNotFoundError:
        print(f"File {file} did not exist")

    for directory in cls.directory_list:
      shutil.rmtree(directory, ignore_errors=True)
      print(f"Removed directory {directory}")

  @classmethod
  def copy_src(cls):
    # Set up the necessary source files
    rootdir = pathlib.Path.cwd().joinpath("../../..").resolve()

    for file in cls.file_list:
      shutil.copy(rootdir.joinpath(file), pathlib.Path.cwd().joinpath(file))

    for directory in cls.directory_list:
      shutil.copytree(rootdir.joinpath(directory), pathlib.Path.cwd().joinpath(directory), dirs_exist_ok=True)


CopyAndSdist.file_list = [
  'AUTHORS',
  'ChangeLog',
  'CMakeLists.txt',
  'config_cmake.h_in',
  'COPYING',
  'INSTALL',
  'INSTALL.md',
  'LICENSE',
  'NEWS',
  'README.md'
]
# Directories
CopyAndSdist.directory_list = [
  'cmake',
  'coefs',
  'exputil',
  'extern/rapidxml',
  'include'
]

# Install specifications

setup(
    name="EXP-libraries",
    version="7.7.28",
    description="Nbody EXPansion Code - Libraries",
    author="",
    license="GPL-3.0",
    packages=["EXP-libexputil", "EXP-libexpcoefs"],
    python_requires=">=3.8",
    package_dir={
      "EXP-libexputil": "exputil",
      "EXP-libexpcoefs": "coefs"
      },
    cmdclass={
      "sdist": CopyAndSdist,
    },
    cmake_install_dir="",
    cmake_minimum_required_version="3.21",
    cmake_languages=("C", "CXX", "Fortran"),
    cmake_args=[
      "-DCMAKE_CXX_STANDARD=17",
      "-DINSTALL_HEADERS=ON",
      "-DINSTALL_CMAKE_FIND=ON",
      "-DENABLE_NBODY=OFF",
      "-DENABLE_PYEXP=OFF",
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

