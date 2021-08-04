#
# Copyright (c) 2020, CINECA
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#		  * Redistributions of source code must retain the above copyright notice, this
#				list of conditions and the following disclaimer.
#
#			* Redistributions in binary form must reproduce the above copyright notice,
#				this list of conditions and the following disclaimer in the documentation
#				and/or other materials provided with the distribution.
#
#			* Neither the name of the copyright holder nor the names of its
#				contributors may be used to endorse or promote products derived from
#				this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Author: Daniele Cesarini, CINECA

#[=======================================================================[
FindSLURM
-------

Finds the SLURM library.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported targets, if found:

``SLURM::SLURM``
  The Slurm library

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``SLURM_FOUND``
  True if the system has the SLURM library.
``SLURM_INCLUDE_DIRS``
  Include directories needed to use SLURM.
``SLURM_LIBRARIES``
  Libraries needed to link to SLURM.

Cache Variables
^^^^^^^^^^^^^^^

The following cache variables may also be set:

``SLURM_INCLUDE_DIR``
  The directory containing ``slurm.h`` and ``spank.h``.
``SLURM_LIBRARY``
  The path to the SLURM library.

#]=======================================================================]

find_package(PkgConfig)
pkg_check_modules(PC_SLURM QUIET SLURM)

find_path(
  SLURM_INCLUDE_DIR
  NAMES
    "slurm/slurm.h"
    "slurm/spank.h"
  PATHS
    ENV CPATH
    ENV C_INCLUDE_PATH
    ENV CPLUS_INCLUDE_PATH
    ENV SLURM_PATH
  PATH_SUFFIXES
    "include"
  DOC
    "Path to the SLURM shared library")

find_library(
  SLURM_LIBRARY
  NAMES
    slurm
  PATHS
    ENV LD_LIBRARY_PATH
    ENV LIBRARY_PATH
    ENV PATH
    ENV SLURM_PATH
  PATH_SUFFIXES
    "lib"
    "lib64"
  DOC
    "Path to the SLURM include directory")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SLURM
  FOUND_VAR 
    SLURM_FOUND
  REQUIRED_VARS
    SLURM_LIBRARY
    SLURM_INCLUDE_DIR)

if(SLURM_FOUND)
  set(SLURM_LIBRARIES "${SLURM_LIBRARY}")
  set(SLURM_INCLUDE_DIRS "${SLURM_INCLUDE_DIR}")
  set(SLURM_DEFINITIONS "${PC_SLURM_CFLAGS_OTHER}")

  add_library(SLURM::SLURM UNKNOWN IMPORTED)
  set_target_properties(SLURM::SLURM 
    PROPERTIES
      IMPORTED_LOCATION "${SLURM_LIBRARY}"
      INTERFACE_COMPILE_OPTIONS "${PC_SLURM_CFLAGS_OTHER}"
      INTERFACE_INCLUDE_DIRECTORIES "${SLURM_INCLUDE_DIR}")
endif()
