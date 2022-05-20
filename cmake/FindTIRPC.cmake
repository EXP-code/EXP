# FindTIRPC
# ---------
#
# Find the native TIRPC includes and library.
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This module will set the following variables in your project:
#
# ``TIRPC_INCLUDE_DIRS``
#   where to find rpc.h, etc.
# ``TIRPC_LIBRARIES``
#   the libraries to link against to use TIRPC.
# ``TIRPC_VERSION``
#   the version of TIRPC found.
# ``TIRPC_FOUND``
#   true if the TIRPC headers and libraries were found.
#

find_package(PkgConfig)
pkg_check_modules(PC_TIRPC QUIET TIRPC)

find_path(
  TIRPC_INCLUDE_DIR
  NAMES
  "rpc/types.h"
  "tpc/xdr.h"
  PATHS
  ENV CPATH
  ENV C_INCLUDE_PATH
  ENV CPLUS_INCLUDE_PATH
  ENV RPC_PATH
  PATH_SUFFIXES
  "include/tirpc"
  DOC
  "Path to the Sun RPC include directory"
)

find_library(TIRPC_LIBRARY
  NAMES tirpc
  PATH_SUFFIXES
  "lib"
  "lib64"
  "lib/x86_64-linux-gnu"
  DOC
  "Path to the Sun RPC library"
)

set(TIRPC_VERSION ${PC_TIRPC_VERSION})

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(TIRPC
  FOUND_VAR
  TIRPC_FOUND
  REQUIRED_VARS
  TIRPC_LIBRARY
  TIRPC_INCLUDE_DIR
  )

if(TIRPC_FOUND)
  set(TIRPC_LIBRARIES "${SLURM_LIBRARY}")
  set(TIRPC_INCLUDE_DIRS "${SLURM_INCLUDE_DIR}")
  set(TIRPC_DEFINITIONS "${PC_SLURM_CFLAGS_OTHER}")

  add_library(TIRPC::TIRPC UNKNOWN IMPORTED)
  set_target_properties(TIRPC::TIRPC
    PROPERTIES
      IMPORTED_LOCATION "${TIRPC_LIBRARY}"
      INTERFACE_COMPILE_OPTIONS "${PC_TIRPC_CFLAGS_OTHER}"
      INTERFACE_INCLUDE_DIRECTORIES "${TIRPC_INCLUDE_DIR}")
endif()

