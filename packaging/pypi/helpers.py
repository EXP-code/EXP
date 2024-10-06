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


