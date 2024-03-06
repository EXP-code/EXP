find_package(EXPutil)
find_package(EXPcoefs)
find_package(EXPcommon)

set(CUDA_SRC)
if (ENABLE_CUDA)
  list(APPEND CUDA_SRC cudaPolarBasis.cu cudaSphericalBasis.cu
    cudaCylinder.cu cudaEmpCylSL.cu cudaComponent.cu NVTX.cc
    cudaHOT.cu cudaIncpos.cu cudaIncvel.cu cudaMultistep.cu
    cudaOrient.cu cudaBiorthCyl.cu cudaCube.cu)
endif()

set(exp_SOURCES Basis.cc Bessel.cc CBrock.cc Component.cc
  CBrockDisk.cc Cube.cc Cylinder.cc ExternalCollection.cc
  ExternalForce.cc Hernquist.cc Orient.cc PotAccel.cc ScatterMFP.cc
  PeriodicBC.cc SphericalBasis.cc AxisymmetricBasis.cc Sphere.cc
  TwoDCoefs.cc TwoCenter.cc EJcom.cc global.cc begin.cc ddplgndr.cc
  Direct.cc Shells.cc NoForce.cc end.cc OutputContainer.cc OutPS.cc
  OutPSQ.cc OutPSN.cc OutPSP.cc OutPSR.cc OutCHKPT.cc OutCHKPTQ.cc
  Output.cc externalShock.cc CylEXP.cc generateRelaxation.cc
  HaloBulge.cc incpos.cc incvel.cc ComponentContainer.cc OutAscii.cc
  OutMulti.cc OutRelaxation.cc OrbTrace.cc OutDiag.cc OutLog.cc
  OutVel.cc OutCoef.cc multistep.cc parse.cc Slab.cc SlabSL.cc step.cc
  tidalField.cc ultra.cc ultrasphere.cc MPL.cc OutFrac.cc OutCalbr.cc
  ParticleFerry.cc pCell.cc chkSlurm.c chkTimer.cc pHOT.cc
  GravKernel.cc ${CUDA_SRC} CenterFile.cc PolarBasis.cc FlatDisk.cc
  signals.cc)

set(common_INCLUDE_DIRS $<INSTALL_INTERFACE:include>
  ${HDF5_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR}
  ${EXPutil_INCLUDE_DIRS} ${EXPcoefs_INCLUDE_DIRS} ${EXPcommon_INCLUDE_DIRS})

set(common_INCLUDE_DIRS_PRIVATE
    ${CMAKE_BINARY_DIR}/config
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/extern/>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

if (NOT (EXPutil_FOUND AND EXPcoefs_FOUND))
  list(APPEND common_INCLUDE_DIRS_PRIVATE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
endif()

# set(common_INCLUDE_DIRS
#   $<INSTALL_INTERFACE:include>
#   $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
#   $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/pyEXP/>
#   $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/extern/DSMC/src>
#   $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/extern/yaml-cpp/include>
#   ${CMAKE_BINARY_DIR} ${DEP_INC}
#   ${CMAKE_CURRENT_SOURCE_DIR}
#   ${EIGEN3_INCLUDE_DIR})

set(DSMC_LIBS)
if (ENABLE_DSMC)
  add_subdirectory(${PROJECT_SOURCE_DIR}/extern/DSMC/src DSMC)
  list(APPEND DSMC_LIBS expdsmc)
endif()

set(common_LINKLIB ${DSMC_LIBS}
  $<IF:$<BOOL:${EXPutil_LIBRARIES}>,${EXPutil_LIBRARIES},exputil>
  $<IF:$<BOOL:${EXPcoefs_LIBRARIES}>,${EXPcoefs_LIBRARIES},expcoef>
  ${OpenMP_CXX_LIBRARIES} ${MPI_CXX_LIBRARIES}
  ${YAML_CPP_LIBRARIES} ${VTK_LIBRARIES})

if(PNG_FOUND)
  list(APPEND common_LINKLIB ${PNG_LIBRARIES})
endif()

if(ENABLE_CUDA)
  list(APPEND common_INCLUDE_DIRS ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} ${CUDAToolkit_INCLUDE_DIRS})
  list(APPEND common_LINKLIB CUDA::cudart CUDA::nvToolsExt)
endif()

if(SLURM_FOUND)
  list(APPEND common_LINKLIB ${SLURM_LIBRARY})
endif()

if(ENABLE_XDR AND TIRPC_FOUND)
  list(APPEND common_LINKLIB ${TIRPC_LIBRARY})
endif()

add_library(EXPlib ${exp_SOURCES})
set_target_properties(EXPlib PROPERTIES OUTPUT_NAME EXPlib)
set_target_properties (EXPlib PROPERTIES EXPORT_NAME EXPnbody)
target_include_directories(EXPlib PUBLIC ${common_INCLUDE_DIRS})
target_include_directories(EXPlib PRIVATE ${common_INCLUDE_DIRS_PRIVATE})
target_link_libraries(EXPlib PUBLIC ${common_LINKLIB})

if (ENABLE_USER)
  add_subdirectory(user)
endif()

add_executable(exp expand.cc)
target_include_directories(exp PUBLIC ${common_INCLUDE_DIRS})
target_include_directories(exp PRIVATE ${common_INCLUDE_DIRS_PRIVATE})
target_link_libraries(exp PUBLIC ${common_LINKLIB} EXPlib)

install(TARGETS EXPlib exp
  EXPORT EXPnbodyTargets
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include/exputil)

if (INSTALL_CMAKE_FIND)
  # Set up cmake configs
  set(CMAKEPACKAGE_INSTALL_DIR share/EXPnbody/cmake)
  install (EXPORT EXPnbodyTargets
    FILE EXPnbodyTargets.cmake
    NAMESPACE EXP::
    DESTINATION ${CMAKEPACKAGE_INSTALL_DIR})

  include(CMakePackageConfigHelpers)
  configure_package_config_file(EXPnbodyConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/EXPnbodyConfig.cmake
    INSTALL_DESTINATION ${CMAKEPACKAGE_INSTALL_DIR}
    PATH_VARS CMAKEPACKAGE_INSTALL_DIR CMAKE_INSTALL_PREFIX
  )
  write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/EXPnbodyVersion.cmake
                                  VERSION ${CMAKE_PROJECT_VERSION}
                                  COMPATIBILITY SameMinorVersion )
  install (FILES ${CMAKE_CURRENT_BINARY_DIR}/EXPnbodyConfig.cmake
                ${CMAKE_CURRENT_BINARY_DIR}/EXPnbodyVersion.cmake
          DESTINATION ${CMAKEPACKAGE_INSTALL_DIR} )
endif()