cmake_minimum_required(VERSION 3.10)
project(
  EXP
  VERSION 7.697
  HOMEPAGE_URL https://mdweinberg@bitbucket.org/mdweinberg/exp.git
  LANGUAGES C CXX Fortran)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)

# compiler flags
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using intel
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qno-offload")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
endif()

add_compile_options(-D_REENTRANT)

option(DSMC_ENABLED "Enable DSCM module" TRUE)
option(ENABLE_CUDA "Enable CUDA" TRUE)
option(USER_ENABLED "Enable basic user modules" TRUE)
option(USER_ALL "Enable all user modules" FALSE)
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

# required compiler features
find_package(MPI REQUIRED)
find_package(Boost REQUIRED COMPONENTS filesystem program_options random serialization)
find_package(OpenMP)
find_package(VTK REQUIRED)
find_package(HDF5 COMPONENTS C CXX HL)
find_package(SLURM)
find_package(Eigen3 REQUIRED)
if(ENABLE_CUDA)
  find_package(CUDA REQUIRED)
  enable_language(CUDA)
endif()
find_package(PNG)
include(${VTK_USE_FILE}) # I don't like this global setting stuff, but it's not really
# documented otherwise

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
set(DEP_LIB ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

if(OpenMP_FOUND)
  set(HAVE_OMP_H TRUE)
endif()
if(SLURM_FOUND)
  set(HAVE_LIBSLURM TRUE)
endif()
if(VTK_FOUND)
  set(HAVE_VTK TRUE)
endif()
if(CUDA_FOUND)
  set(HAVE_LIBCUDA TRUE)
endif()
if(HDF5_FOUND)
  set(HAVE_HDF5 TRUE)
endif()
if(PNG_FOUND)
  include_directories(${PROJECT_SOURCE_DIR}/extern/png++)
  set(HAVE_LIBPNGPP TRUE)
endif()

# Get the current working branch
execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest abbreviated commit hash of the working branch
execute_process(
  COMMAND git rev-parse HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the build date
execute_process(
  COMMAND date -u +%Y-%m-%d\ %H:%M:%S\ UTC
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE COMPILE_TIME
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

include_directories(${PROJECT_SOURCE_DIR}/extern/yaml-cpp/include)

message("Configuring build for ${GIT_BRANCH}/${GIT_COMMIT} at ${COMPILE_TIME}")

configure_file(config_cmake.h.in config.h)

add_subdirectory(extern/yaml-cpp)
add_subdirectory(exputil)
add_subdirectory(src)
add_subdirectory(utils)

find_package(Doxygen)

option(BUILD_DOCUMENTATION "Create and install the HTML based API
documentation (requires Doxygen)" ${DOXYGEN_FOUND})

if(BUILD_DOCUMENTATION)
  if(NOT DOXYGEN_FOUND)
    message(FATAL_ERROR "Doxygen is needed to build the documentation.")
  endif()

  set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/doc/exp.cfg)

  configure_file(${doxyfile_in} @ONLY)

  message("Doxygen build started.")

  add_custom_target(doc
    COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile_in}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM)

  #    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION     share/doc)
endif()
