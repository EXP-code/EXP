set(bin_PROGRAMS nativetoh5 h5compare viewcoefs h5power makecoefs testread)

set(common_LINKLIB OpenMP::OpenMP_CXX MPI::MPI_CXX yaml-cpp exputil
  ${VTK_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})

if(PNG_FOUND)
  list(APPEND common_LINKLIB PNG::PNG)
endif()

set(common_INCLUDE $<INSTALL_INTERFACE:include>
  # $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  # $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/> ${CMAKE_BINARY_DIR}
  ${DEP_INC}
  # ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..
  # ${CMAKE_BINARY_DIR}/config
  ${HDF5_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} ${FFTW_INCLUDE_DIRS})

set(common_INCLUDE_DIRS_PRIVATE
    ${CMAKE_BINARY_DIR}/config
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${HighFive_SOURCE_DIR}/include
    )

if(ENABLE_CUDA)
  list(APPEND common_INCLUDE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} ${CUDAToolkit_INCLUDE_DIRS})
  list(APPEND common_LINKLIB CUDA::cudart CUDA::nvToolsExt)
endif()

if(SLURM_FOUND)
  list(APPEND common_LINKLIB ${SLURM_LIBRARY})
endif()

if(ENABLE_XDR AND TIRPC_FOUND)
  list(APPEND common_INCLUDE ${TIRPC_INCLUDE_DIRS})
  list(APPEND common_LINKLIB ${TIRPC_LIBRARIES})
endif()

# Make the expcoef shared library
#
set(expcoef_SOURCES BasisFactory.cc BiorthBasis.cc FieldBasis.cc
  CoefContainer.cc CoefStruct.cc FieldGenerator.cc expMSSA.cc
  Coefficients.cc KMeans.cc Centering.cc ParticleIterator.cc
  Koopman.cc)
add_library(expcoef ${expcoef_SOURCES})

# TODO: Replace this globbing with a list of actually public headers
file(GLOB expcoef_PUBLIC_HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.[hH]")
set_target_properties (expcoef PROPERTIES OUTPUT_NAME expcoef)
set_target_properties (expcoef PROPERTIES EXPORT_NAME EXPcoefs)
set_target_properties (expcoef PROPERTIES PUBLIC_HEADER "${expcoef_PUBLIC_HEADERS}")

target_include_directories(expcoef PUBLIC ${common_INCLUDE})
target_include_directories(expcoef PRIVATE ${common_INCLUDE_DIRS_PRIVATE})
target_link_libraries(expcoef PUBLIC ${common_LINKLIB})

install (TARGETS expcoef
  EXPORT EXPcoefsTargets
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include/expcoefs)

# Set up cmake configs
set(CMAKEPACKAGE_INSTALL_DIR share/EXPcoefs/cmake)
install (EXPORT EXPcoefsTargets
  FILE EXPcoefsTargets.cmake
  NAMESPACE EXPcoefs::
  DESTINATION ${CMAKEPACKAGE_INSTALL_DIR})

include(CMakePackageConfigHelpers)
configure_package_config_file(EXPcoefsConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/EXPcoefsConfig.cmake
  INSTALL_DESTINATION ${CMAKEPACKAGE_INSTALL_DIR}
)

write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/EXPcoefsVersion.cmake
                                 VERSION ${CMAKE_PROJECT_VERSION}
                                 COMPATIBILITY SameMinorVersion )

install (FILES ${CMAKE_CURRENT_BINARY_DIR}/EXPcoefsConfig.cmake
               ${CMAKE_CURRENT_BINARY_DIR}/EXPcoefsVersion.cmake
        DESTINATION ${CMAKEPACKAGE_INSTALL_DIR} )

# Configure and build the test routines
#
if (ENABLE_TESTS)
  add_executable(nativetoh5 coefstoh5.cc)
  add_executable(h5compare  h5compare.cc)
  add_executable(viewcoefs  viewcoefs.cc)
  add_executable(h5power    h5power.cc)
  add_executable(makecoefs  makecoefs.cc)
  add_executable(testread   testread.cc)

  foreach(program ${bin_PROGRAMS})
    target_link_libraries(${program} expcoef exputil ${common_LINKLIB})
    target_include_directories(${program} PUBLIC ${common_INCLUDE})
    target_compile_options(${program} PUBLIC ${OpenMP_CXX_FLAGS})
  endforeach()

  install(TARGETS ${bin_PROGRAMS} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()