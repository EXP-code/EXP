find_package(EXPutil)
find_package(EXPui)

if(EXPutil_FOUND)
  set(EXPUTIL EXP::EXPutil)
  set(EXPUTIL_INCLUDE ${EXPutil_INCLUDE_DIRS})
else()
  set(EXPUTIL exputil) #internal name
  set(EXPUTIL_INCLUDE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/exputil/>)
endif()

if(EXPui_FOUND)
  set(EXPUI EXP::EXPui)
  set(EXPUI_INCLUDE ${EXPui_INCLUDE_DIRS})
else()
  set(EXPUI expui) #internal name
  set(EXPUI_INCLUDE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/expui/>)
endif()

set(CUDA_SRC)
if (ENABLE_CUDA)
  list(APPEND CUDA_SRC cudaPolarBasis.cu cudaSphericalBasis.cu
    cudaCylinder.cu cudaEmpCylSL.cu cudaComponent.cu NVTX.cc
    cudaIncpos.cu cudaIncvel.cu cudaMultistep.cu cudaOrient.cu
    cudaBiorthCyl.cu cudaCube.cu cudaSlabSL.cu)
endif()

set(exp_SOURCES Basis.cc Bessel.cc Component.cc
  Cube.cc Cylinder.cc ExternalCollection.cc
  ExternalForce.cc Orient.cc PotAccel.cc ScatterMFP.cc
  PeriodicBC.cc SphericalBasis.cc AxisymmetricBasis.cc Sphere.cc
  TwoDCoefs.cc TwoCenter.cc EJcom.cc global.cc begin.cc ddplgndr.cc
  Direct.cc Shells.cc NoForce.cc end.cc OutputContainer.cc OutPS.cc
  OutPSQ.cc OutPSN.cc OutPSP.cc OutPSR.cc OutCHKPT.cc OutCHKPTQ.cc
  Output.cc externalShock.cc CylEXP.cc generateRelaxation.cc
  HaloBulge.cc incpos.cc incvel.cc ComponentContainer.cc OutAscii.cc
  OutMulti.cc OutRelaxation.cc OrbTrace.cc OutDiag.cc OutLog.cc
  OutVel.cc OutCoef.cc multistep.cc parse.cc SlabSL.cc step.cc
  tidalField.cc ultra.cc ultrasphere.cc MPL.cc OutFrac.cc OutCalbr.cc
  ParticleFerry.cc chkSlurm.c chkTimer.cc GravKernel.cc ${CUDA_SRC}
  CenterFile.cc PolarBasis.cc FlatDisk.cc signals.cc)

# set(exp_HEADERS AxisymmetricBasis.H
# Basis.H
# Bessel.H
# Component.H
# ComponentContainer.H
# Cube.H
# CylEXP.H
# Cylinder.H
# Direct.H
# EJcom.H
# EnumBitset.H
# ExternalCollection.H
# ExternalForce.H
# FlatDisk.H
# GravKernel.H
# HaloBulge.H
# MixtureBasis.H
# NVTX.H
# NoForce.H
# OrbTrace.H
# Orient.H
# OutAscii.H
# OutCHKPT.H
# OutCHKPTQ.H
# OutCalbr.H
# OutCoef.H
# OutDiag.H
# OutFrac.H
# OutLog.H
# OutMulti.H
# OutPS.H
# OutPSN.H
# OutPSP.H
# OutPSQ.H
# OutPSR.H
# OutRelaxation.H
# OutVel.H
# Output.H
# OutputContainer.H
# ParticleFerry.H
# PeriodicBC.H
# PolarBasis.H
# PotAccel.H
# ScatterMFP.H
# Shells.H
# SlabSL.H
# Sphere.H
# SphericalBasis.H
# StatsMPI.H
# TwoCenter.H
# TwoDCoefs.H
# chkTimer.H
# coef.H
# expand.H
# externalShock.H
# generateRelaxation.H
# global.H
# global_key_set.H
# tidalField.H)

set(common_INCLUDE 
  $<INSTALL_INTERFACE:include/exputil>
  ${EXPUTIL_INCLUDE} 
  ${EXPUI_INCLUDE}
  ${HDF5_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})

set(common_INCLUDE_PRIVATE
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/config/>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/extern/>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

# if (NOT (EXPutil_FOUND AND EXPui_FOUND))
#   list(APPEND common_INCLUDE_PRIVATE

#     $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
# endif()

# set(common_INCLUDE
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

set(common_LINKLIBS
  ${EXPUTIL} ${EXPUI}
  ${DSMC_LIBS}
  $<$<BOOL:${HighFive_FOUND}>:HighFive>
  ${OpenMP_CXX_LIBRARIES} ${MPI_CXX_LIBRARIES}
  ${YAML_CPP_LIBRARIES} ${VTK_LIBRARIES})

if(PNG_FOUND)
  list(APPEND common_LINKLIBS PNG::PNG)
  list(APPEND common_INCLUDE ${PNG_INCLUDE_DIRS})
  add_definitions(${PNG_DEFINITIONS})
endif()

if(ENABLE_CUDA)
  list(APPEND common_INCLUDE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} ${CUDAToolkit_INCLUDE_DIRS})
  list(APPEND common_LINKLIBS CUDA::cudart CUDA::nvToolsExt)
endif()

if(SLURM_FOUND)
  list(APPEND common_LINKLIBS ${SLURM_LIBRARY})
endif()

if(ENABLE_XDR AND TIRPC_FOUND)
  list(APPEND common_LINKLIBS ${TIRPC_LIBRARY})
endif()


add_library(EXPlib ${exp_SOURCES})

set_target_properties(EXPlib PROPERTIES OUTPUT_NAME EXPlib)
set_target_properties (EXPlib PROPERTIES EXPORT_NAME EXPnbody)
target_include_directories(EXPlib PUBLIC ${common_INCLUDE})
target_include_directories(EXPlib PRIVATE ${common_INCLUDE_PRIVATE})
target_link_libraries(EXPlib PUBLIC ${common_LINKLIBS})

if (ENABLE_USER)
  add_subdirectory(user)
endif()

add_executable(exp expand.cc expand.H)
target_include_directories(exp PRIVATE ${common_INCLUDE_PRIVATE})
target_include_directories(exp PUBLIC ${common_INCLUDE})
target_link_libraries(exp PUBLIC EXPlib ${common_LINKLIBS})
target_compile_options(exp PUBLIC ${OpenMP_CXX_FLAGS})

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