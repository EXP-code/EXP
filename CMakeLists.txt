cmake_minimum_required(VERSION 3.10)
project(
  EXP
  VERSION 7.699
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

option(ENABLE_CUDA "Enable CUDA" FALSE)
option(ENABLE_DSMC "Enable DSCM module" FALSE)
option(ENABLE_USER "Enable basic user modules" ON)
option(ENABLE_USER_ALL "Enable all user modules" FALSE)
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)
option(BUILD_DOCS "Build documentation" OFF)

# required compiler features
find_package(MPI REQUIRED)
find_package(Boost REQUIRED COMPONENTS filesystem program_options random serialization)
find_package(OpenMP)
find_package(VTK)
find_package(HDF5 COMPONENTS C CXX HL)
find_package(SLURM)
find_package(Eigen3 REQUIRED)
if(ENABLE_CUDA)
  find_package(CUDA REQUIRED)
  enable_language(CUDA)
  find_cuda_helper_libs(nvToolsExt)
  cuda_select_nvcc_arch_flags(ARCH_FLAGS 5.2 6.0 6.1 7.0 7.2 7.5)
  list(APPEND CUDA_NVCC_FLAGS ${ARCH_FLAGS})
  message("CUDA arch flags are: ${ARCH_FLAGS}")
endif()
find_package(PNG)

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
  include(${VTK_USE_FILE}) # I don't like this global setting stuff,
			   # but it's not really documented otherwise
endif()
if(CUDA_FOUND)
  set(HAVE_LIBCUDA TRUE)
endif()
if(HDF5_FOUND)
  set(HAVE_HDF5 TRUE)
endif()
if(PNG_FOUND)
  include_directories(${PROJECT_SOURCE_DIR}/extern)
  set(HAVE_LIBPNGPP TRUE)
endif()
if(ENABLE_DSMC)
  set(DSMC_ENABLED 1)
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

configure_file(${CMAKE_SOURCE_DIR}/config_cmake.h.in ${CMAKE_SOURCE_DIR}/config.h)

add_subdirectory(extern/yaml-cpp)
add_subdirectory(exputil)
add_subdirectory(src)
add_subdirectory(utils)

# Force installation of the yaml-cpp libraries
install(TARGETS yaml-cpp DESTINATION lib)

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
