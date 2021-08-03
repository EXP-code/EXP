cmake_minimum_required(VERSION 3.10)
project(
  EXP
  VERSION 1.0.0
  LANGUAGES C CXX CUDA)

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

set(MPI_CXX_SKIP_MPICXX
  true
  CACHE BOOL
  "If true, the MPI-2 C++ bindings are disabled using definitions.")

# required compiler features
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
find_package(VTK REQUIRED)
include(${VTK_USE_FILE}) # I don't like this global setting stuff, but it's not really
# documented otherwise

# FIXME: Add variable CUDA_SRC (and detect cuda, obviously)

find_package(MKL)
if(MKL_FOUND)
  set(DEP_LIB ${MKL_SDL_LIBRARY})
else()
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  set(DEP_LIB ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES})
endif()

add_subdirectory(extern/yaml-cpp)

add_subdirectory(exputil)
