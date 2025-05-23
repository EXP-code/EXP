cmake_minimum_required(VERSION 3.25) # Needed for CUDA, MPI, and CTest features

project(
  EXP
  VERSION "7.8.6"
  HOMEPAGE_URL https://github.com/EXP-code/EXP
  LANGUAGES C CXX Fortran)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler flags.  Not all tested thoroughly...
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # using Clang
  message(STATUS "Clang compiler version is: ${CMAKE_CXX_COMPILER_VERSION}")
  if(CMAKE_SYSTEM_NAME MATCHES "Linux")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 18.0.0)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++")
    endif()
  endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  # using GCC
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  # using intel
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qno-offload")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  # using Visual Studio C++
endif()

# Required compiler features
add_compile_options(-D_REENTRANT)

# Bake in library paths (esp. useful for HPC sites with modules)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Check and enforce that we are a git repository.  Necessary for
# submodules to work correctly.
if(EXISTS "${PROJECT_SOURCE_DIR}/.git")
  message(STATUS "Checking that we are a git repository - good")
else()
  message(STATUS "Checking that we are a git repository - NO")
  message(SEND_ERROR "You need to 'git clone ${CMAKE_PROJECT_HOMEPAGE_URL}'.  Please don't use the zip download.")
endif()

# Build options
option(ENABLE_NBODY "Enable EXP n-body" ON)
option(ENABLE_PYEXP "Enable the Python bindings" ON)
option(ENABLE_PNG "Enable PNG graphics support" FALSE)
option(ENABLE_CUDA "Enable CUDA" FALSE)
option(ENABLE_SLURM "Enable SLURM checkpointing support" FALSE)
option(ENABLE_XDR "Enable RPC/XDR support for Tipsy standard" FALSE)
option(ENABLE_VTK "Configure VTK if available" FALSE)
option(ENABLE_CUDA_SINGLE "Use real*4 instead of real*8 for CUDA" FALSE)
option(ENABLE_USER "Enable basic user modules" ON)
option(ENABLE_SLCHECK "Enable *careful* Sturm-Liouville solutions" TRUE)
option(ENABLE_TESTS "Enable build tests for EXP, pyEXP and helpers" ON)
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)
option(BUILD_DOCS "Build documentation" OFF)

# Set mpirun launcher for CTest

set(EXP_MPI_LAUNCH "mpirun" CACHE STRING "Command to run an MPI application (for unit tests only)")

# Find newest version if multiple versions are available

set(CMAKE_FIND_PACKAGE_SORT_ORDER NATURAL)
set(CMAKE_FIND_PACKAGE_SORT_DIRECTION DEC)

# Package support

find_package(MPI REQUIRED COMPONENTS C CXX)
find_package(OpenMP)
find_package(FFTW REQUIRED)
find_package(HDF5 COMPONENTS C CXX HL REQUIRED)
find_package(TIRPC)	       # Check for alternative Sun rpc support
find_package(Eigen3 3.4...<3.5 REQUIRED)
find_package(PNG)
find_package(Git)

# Check for FE
include(FEENABLE)

# Check for Slurm if user has not set to false
if(ENABLE_SLURM)
  find_package(SLURM)
endif()

# Optional VTK support
if(ENABLE_VTK)
  find_package(VTK)
endif()

# Optional CUDA support
if(ENABLE_CUDA)
  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)
  message(STATUS "Using CUDA architectures ${CMAKE_CUDA_ARCHITECTURES} by default.  Please set the CUDAARCHS environment variable to override this.")
  set(CMAKE_CUDA_FLAGS_RELEASE "-O3 -w"
    CACHE STRING "Global cuda flags for Release build" FORCE)
  set(CMAKE_CUDA_FLAGS_MINSIZEREL "-O3 -w"
    CACHE STRING "Global cuda flags for MinSizeRel build" FORCE)
  set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-O2 -g -w"
    CACHE STRING "Global cuda flags for RelWithDebInfo build" FORCE)
  set(HAVE_LIBCUDA TRUE)
endif()

# For mpi
include_directories(${MPI_C_HEADER_DIR})

# HDF5 version restrictions
if(HDF5_VERSION VERSION_LESS 1.8)
  message(FATAL_ERROR "Your HDF5 version must be >=1.8.20 or >=1.10.2 or >=1.12")
endif()

if(HDF5_VERSION VERSION_GREATER_EQUAL 1.8.0 AND
    HDF5_VERSION VERSION_LESS 1.8.20)
  message(FATAL_ERROR "Your HDF5 version in 1.8 series must be 1.8.20 or greater")
endif()

if(HDF5_VERSION VERSION_GREATER_EQUAL 1.10.0 AND
    HDF5_VERSION VERSION_LESS 1.10.2)
  message(FATAL_ERROR "Your HDF5 version in 1.10 series must be 1.10.2 or greater")
endif()

# For git submodules
include_directories(${PROJECT_SOURCE_DIR}/extern)
# Add OpenMP compile flags
if(OpenMP_FOUND)
  set(HAVE_OMP_H TRUE)
  OPTION (USE_OpenMP "Use OpenMP" ON)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  if(ENABLE_CUDA)
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler='${OpenMP_CXX_FLAGS}'")
  endif()
endif()
# Slurm support
if(SLURM_FOUND)
  set(HAVE_LIBSLURM TRUE)
endif()
# Check for VTK
if(VTK_FOUND)
  set(HAVE_VTK TRUE)
  if (VTK_VERSION VERSION_LESS "8.90.0")
    include(${VTK_USE_FILE}) # I don't like this global setting stuff
  else()
    include_directories("${VTK_PREFIX_PATH}/include/vtk-${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}")
  endif()
endif()
if(HDF5_FOUND)
  set(HAVE_HDF5 TRUE)
endif()
if(FFTW_FOUND)
  set(HAVE_FFTW TRUE)
endif()
if(ENABLE_SLCHECK)
  set(SLEDGE_THROW TRUE)
endif()
if(PNG_FOUND AND ENABLE_PNG)
  set(HAVE_LIBPNGPP TRUE)
endif()
if(ENABLE_CUDA_SINGLE)
  add_compile_definitions(O_SINGLE=1)
endif()

# Checking for RPC support (needed for Tipsy standard format)
if(NOT TIRPC_FOUND)
  include(CheckIncludeFile)
  check_include_file("rpc/types.h" HAVE_RPC_TYPES)
endif()

# Only include RPC support if the installer wants XDR
if(ENABLE_XDR)
  if(TIRPC_FOUND OR HAVE_RPC_TYPES)
    set(HAVE_XDR TRUE CACHE BOOL "We have an XDR implementation")
    message(STATUS "We have an XDR implementation; Tipsy standard files will be supported.")
    if(TIRPC_FOUND)		# Add TIRPC support
      set(HAVE_TIRPC TRUE)
      include_directories(${TIRPC_INCLUDE_DIRS})
    endif()
  else()
    message(STATUS "Could NOT find RPC support; you will not be able to read Tipsy standard files.  Tipsy native (and Bonsai) will still work.")
  endif()
else()
  message(STATUS "RPC/XDR is disabled; you will not be able to read Tipsy standard files.  Tipsy native (and Bonsai) will still work.")
endif()

# For setting the default size of the Cuda real attribute array
set(CUDA_EXP_DATTRIB "4" CACHE STRING
  "Number of real particle attributes for Cuda particle structure")
add_compile_definitions(DATTRIB_CUDA=${CUDA_EXP_DATTRIB})

# Get the current working branch
execute_process(
  COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Git submodule updates
execute_process(
  COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  RESULT_VARIABLE GIT_SUBMOD_RESULT
)

if(NOT GIT_SUBMOD_RESULT EQUAL "0")
   message(FATAL_ERROR "git submodule update --init --recursive failed ${GIT_SUBMOD_RESULT}, please checkout submodules")
else()
   message(STATUS "Submodules updated successfully - good")
endif()

# Get the latest abbreviated commit hash of the working branch
execute_process(
  COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the build date
execute_process(
  COMMAND date +%Y-%m-%d\ %H:%M:%S\ %Z
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE COMPILE_TIME
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

include_directories(${PROJECT_SOURCE_DIR}/extern/yaml-cpp/include)
include_directories(${PROJECT_SOURCE_DIR}/extern/pybind11/include)

# Report to the user
message("Configuring build for ${GIT_BRANCH}/${GIT_COMMIT} at ${COMPILE_TIME}")

add_subdirectory(extern/yaml-cpp)
add_subdirectory(extern/pybind11)

# Set options for the HighFive git submodule in extern
set(HIGHFIVE_EXAMPLES OFF CACHE BOOL "Do not build the examples")
set(HIGHFIVE_BUILD_DOCS OFF CACHE BOOL "Do not build the documentation")
set(HIGHFIVE_USE_BOOST OFF CACHE BOOL "Do not use Boost in HighFive")
set(HIGHFIVE_UNIT_TESTS OFF CACHE BOOL "Turn off internal testing for HighFIve")
set(H5_USE_EIGEN TRUE CACHE BOOL "Eigen3 support in HighFive")

add_subdirectory(extern/HighFive EXCLUDE_FROM_ALL)

# Configure the remaining native subdirectories
add_subdirectory(exputil)
add_subdirectory(expui)
if (ENABLE_NBODY)
  add_subdirectory(src)
endif()
add_subdirectory(utils)
if (ENABLE_PYEXP)
  add_subdirectory(pyEXP)
endif()

add_subdirectory(extern/user-modules)

# Build the tests; set ENABLE_TEST=OFF to disable
if(ENABLE_TESTS)
  include(CTest)
  add_subdirectory(tests)
endif()

# try to find pybind11 and build wrapper python module
find_package(Python3 COMPONENTS Interpreter Development)
message(STATUS "python3 include dirs: ${Python3_INCLUDE_DIRS}")


# Force installation of the yaml-cpp libraries
install(TARGETS yaml-cpp DESTINATION lib)

# Check for doxygen is the user wants web docs
if (BUILD_DOCS)
  find_package(Doxygen REQUIRED)
endif()

if(DOXYGEN_FOUND)
  set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/doc/exp.cfg)

  message("Making documenation using Doxygen.")

  add_custom_target(doc ALL
    COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile_in}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM)

  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc/html
    DESTINATION share/EXP/doc)
endif()

if (NOT DOXYGEN_FOUND AND BUILD_DOCS)
  message("Doxygen is needed to build the documentation. This is optional.")
endif()

# Add sanitizer build types
set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE}
  CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel Tsan Asan Lsan Msan Ubsan"
  FORCE)

# ThreadSanitizer
set(CMAKE_C_FLAGS_TSAN
  "-fsanitize=thread -g -O1"
  CACHE STRING "Flags used by the C compiler during ThreadSanitizer builds."
  FORCE)

set(CMAKE_CXX_FLAGS_TSAN
  "-fsanitize=thread -g -O1"
  CACHE STRING "Flags used by the C++ compiler during ThreadSanitizer builds."
  FORCE)

# AddressSanitize
set(CMAKE_C_FLAGS_ASAN
  "-fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer -g -O1"
  CACHE STRING "Flags used by the C compiler during AddressSanitizer builds."
  FORCE)
set(CMAKE_CXX_FLAGS_ASAN
  "-fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer -g -O1"
  CACHE STRING "Flags used by the C++ compiler during AddressSanitizer builds."
  FORCE)

# LeakSanitizer
set(CMAKE_C_FLAGS_LSAN
  "-fsanitize=leak -fno-omit-frame-pointer -g -O1"
  CACHE STRING "Flags used by the C compiler during LeakSanitizer builds."
  FORCE)
set(CMAKE_CXX_FLAGS_LSAN
  "-fsanitize=leak -fno-omit-frame-pointer -g -O1"
  CACHE STRING "Flags used by the C++ compiler during LeakSanitizer builds."
  FORCE)

# MemorySanitizer
set(CMAKE_C_FLAGS_MSAN
  "-fsanitize=memory -fno-optimize-sibling-calls -fsanitize-memory-track-origins=2 -fno-omit-frame-pointer -g -O2"
  CACHE STRING "Flags used by the C compiler during MemorySanitizer builds."
  FORCE)
set(CMAKE_CXX_FLAGS_MSAN
  "-fsanitize=memory -fno-optimize-sibling-calls -fsanitize-memory-track-origins=2 -fno-omit-frame-pointer -g -O2"
  CACHE STRING "Flags used by the C++ compiler during MemorySanitizer builds."
  FORCE)

# UndefinedBehaviour
set(CMAKE_C_FLAGS_UBSAN
  "-fsanitize=undefined"
  CACHE STRING "Flags used by the C compiler during UndefinedBehaviourSanitizer builds."
  FORCE)
set(CMAKE_CXX_FLAGS_UBSAN
  "-fsanitize=undefined"
  CACHE STRING "Flags used by the C++ compiler during UndefinedBehaviourSanitizer builds."
  FORCE)

# Make the config_exp.h file in the build directory, add the build dir
# to the CMake include path.  This allows for multiple configurations
# from the same source.
configure_file(${CMAKE_SOURCE_DIR}/config_cmake.h_in ${CMAKE_BINARY_DIR}/config_exp.h)
include_directories(${PROJECT_BINARY_DIR})

