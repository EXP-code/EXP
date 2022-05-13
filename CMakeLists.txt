cmake_minimum_required(VERSION 3.12)
project(
  EXP
  VERSION 7.705
  HOMEPAGE_URL https://mdweinberg@bitbucket.org/mdweinberg/exp.git
  LANGUAGES C CXX Fortran)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_CXX_STANDARD 17)
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

# Required compiler features
add_compile_options(-D_REENTRANT)

# Check that we are a git repository
if(EXISTS "${PROJECT_SOURCE_DIR}/.git")
  message(STATUS "Checking that we are a git repository - good")
else()
  message(STATUS "Checking that we are a git repository - NO")
  message(SEND_ERROR "You need to 'git clone ${CMAKE_PROJECT_HOMEPAGE_URL}'.  Please don't use the zip download.")
endif()

# Build options
option(ENABLE_NBODY "Enable EXP n-body" ON)
option(ENABLE_CUDA "Enable CUDA" FALSE)
option(ENABLE_CUDA_SINGLE "Use real*4 instead of real*8 for CUDA" FALSE)
option(ENABLE_DSMC "Enable DSCM module" FALSE)
option(ENABLE_USER "Enable basic user modules" ON)
option(ENABLE_USER_ALL "Enable all user modules" FALSE)
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)
option(BUILD_DOCS "Build documentation" OFF)

# Package support
find_package(MPI REQUIRED)
find_package(Boost COMPONENTS serialization)
find_package(OpenMP)
find_package(FFTW)
find_package(VTK)
find_package(HDF5 COMPONENTS C CXX HL REQUIRED)
find_package(SLURM)
find_package(Eigen3 REQUIRED)
if(ENABLE_CUDA)
  find_package(CUDA REQUIRED)
  enable_language(CUDA)
  find_cuda_helper_libs(nvToolsExt)
  cuda_select_nvcc_arch_flags(ARCH_FLAGS 5.2 6.0 6.1 7.0 7.2 7.5)
  # list(APPEND CUDA_NVCC_FLAGS ${ARCH_FLAGS})
  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} ${ARCH_FLAGS}" CACHE STRING "Global compile flags for nvcc")
  message("CUDA arch flags are: ${ARCH_FLAGS}")
  set(CMAKE_CUDA_FLAGS_RELEASE "-O3 -w"
    CACHE STRING "Global cuda flags for Release build" FORCE)
  set(CMAKE_CUDA_FLAGS_MINSIZEREL "-O3 -w"
    CACHE STRING "Global cuda flags for MinSizeRel build" FORCE)
  set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-O2 -g -w"
    CACHE STRING "Global cuda flags for RelWithDebInfo build" FORCE)
endif()

find_package(PNG)

# Check for Sun rpc support
include(CheckIncludeFile)
include(FindPkgConfig)
pkg_check_modules(TIRPC libtirpc)
check_include_file("rpc/types.h" HAVE_RPC_TYPES)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
set(DEP_LIB ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# For PNG and rapidxml
include_directories(${PROJECT_SOURCE_DIR}/extern)

if(OpenMP_FOUND)
  set(HAVE_OMP_H TRUE)
endif()
if(SLURM_FOUND)
  set(HAVE_LIBSLURM TRUE)
endif()
if(VTK_FOUND)
  set(HAVE_VTK TRUE)
  if (VTK_VERSION VERSION_LESS "8.90.0")
    include(${VTK_USE_FILE}) # I don't like this global setting stuff
  else()
    include_directories("${VTK_PREFIX_PATH}/include/vtk-${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}")
  endif()
endif()
if(CUDA_FOUND)
  set(HAVE_LIBCUDA TRUE)
endif()
if(HDF5_FOUND)
  set(HAVE_HDF5 TRUE)
endif()
if(FFTW_FOUND)
  set(HAVE_FFTW TRUE)
endif()
if(PNG_FOUND)
  set(HAVE_LIBPNGPP TRUE)
endif()
if(ENABLE_DSMC)
  if(NOT BOOST_FOUND)
    message(SEND_ERROR "You need Boost to compile DSMC")
  endif()
  set(DSMC_ENABLED 1)
endif()
if(ENABLE_CUDA_SINGLE)
  add_compile_definitions(O_SINGLE=1)
endif()
if(TIRPC_FOUND OR HAVE_RPC_TYPES)
  set(HAVE_XDR TRUE CACHE BOOL "We have an XDR implementation")
endif()

# For setting the default size of the Cuda real attribute array
set(CUDA_EXP_DATTRIB "0" CACHE STRING
  "Number of real particle attributes for Cuda particle structure")
add_compile_definitions(DATTRIB_CUDA=${CUDA_EXP_DATTRIB})

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
include_directories(${PROJECT_SOURCE_DIR}/extern/pybind11/include)


message("Configuring build for ${GIT_BRANCH}/${GIT_COMMIT} at ${COMPILE_TIME}")

configure_file(${CMAKE_SOURCE_DIR}/config_cmake.h.in ${CMAKE_SOURCE_DIR}/config.h)

add_subdirectory(extern/yaml-cpp)
add_subdirectory(extern/pybind11)

add_subdirectory(exputil)
if (ENABLE_NBODY)
  add_subdirectory(src)
endif()
add_subdirectory(utils)

# try to find pybind11 and build wrapper python module
find_package(Python3 COMPONENTS Interpreter Development)
message(STATUS "python3 include dirs: ${Python3_INCLUDE_DIRS}")


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
