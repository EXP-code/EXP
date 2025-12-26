set(exp_SOURCES Basis.cc Bessel.cc Component.cc Cube.cc Cylinder.cc
  ExternalCollection.cc CBDisk.cc ExternalForce.cc Orient.cc
  PotAccel.cc ScatterMFP.cc PeriodicBC.cc SphericalBasis.cc
  AxisymmetricBasis.cc Sphere.cc TwoDCoefs.cc TwoCenter.cc EJcom.cc
  global.cc begin.cc ddplgndr.cc Direct.cc Shells.cc NoForce.cc end.cc
  OutputContainer.cc OutPS.cc OutPSQ.cc OutPSN.cc OutPSP.cc OutPSR.cc
  OutCHKPT.cc OutCHKPTQ.cc Output.cc externalShock.cc CylEXP.cc
  generateRelaxation.cc HaloBulge.cc incpos.cc incvel.cc
  ComponentContainer.cc OutAscii.cc OutHDF5.cc OutMulti.cc
  OutRelaxation.cc OrbTrace.cc OutDiag.cc OutLog.cc OutVel.cc
  OutCoef.cc OutSample.cc multistep.cc parse.cc SlabSL.cc step.cc
  tidalField.cc ultra.cc ultrasphere.cc MPL.cc OutFrac.cc OutCalbr.cc
  ParticleFerry.cc chkSlurm.c chkTimer.cc GravKernel.cc CenterFile.cc
  PolarBasis.cc FlatDisk.cc signals.cc)

if (ENABLE_CUDA)
  list(APPEND exp_SOURCES cudaPolarBasis.cu cudaSphericalBasis.cu
    cudaCylinder.cu cudaEmpCylSL.cu cudaComponent.cu NVTX.cc
    cudaIncpos.cu cudaIncvel.cu cudaMultistep.cu cudaOrient.cu
    cudaBiorthCyl.cu cudaCube.cu cudaSlabSL.cu cudaPeriodicBC.cu)
endif()

set(common_INCLUDE_DIRS 
  $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/pyEXP/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/extern/yaml-cpp/include>
  ${CMAKE_BINARY_DIR} ${DEP_INC}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${EIGEN3_INCLUDE_DIR})

set(common_LINKLIB exputil expui OpenMP::OpenMP_CXX MPI::MPI_CXX
  yaml-cpp ${VTK_LIBRARIES})

if(PNG_FOUND)
  list(APPEND common_LINKLIB PNG::PNG)
endif()

if(ENABLE_CUDA)
  list(APPEND common_LINKLIB CUDA::toolkit CUDA::cudart)
  if (CUDAToolkit_VERSION VERSION_GREATER_EQUAL 12)
    list(APPEND common_LINKLIB CUDA::nvtx3)
    set_source_files_properties(${exp_SOURCES} PROPERTIES LANGUAGE CUDA)
  else ()
    list(APPEND common_LINKLIB CUDA::nvToolsExt)
  endif ()
endif()

if(SLURM_FOUND)
  list(APPEND common_LINKLIB ${SLURM_LIBRARY})
endif()

if(ENABLE_XDR AND TIRPC_FOUND)
  list(APPEND common_LINKLIB ${TIRPC_LIBRARY})
endif()

add_library(EXPlib ${exp_SOURCES})
set_target_properties(EXPlib PROPERTIES OUTPUT_NAME EXPlib)
target_include_directories(EXPlib PUBLIC ${common_INCLUDE_DIRS})
target_link_libraries(EXPlib PUBLIC ${common_LINKLIB})

if (ENABLE_USER)
  add_subdirectory(user)
endif()

add_executable(exp expand.cc)
target_include_directories(exp PUBLIC ${common_INCLUDE_DIRS})
target_link_libraries(exp PUBLIC ${common_LINKLIB} EXPlib)

if (ENABLE_CUDA)
  set_target_properties(exp PROPERTIES LINKER_LANGUAGE CUDA)
endif ()

install(TARGETS EXPlib DESTINATION lib)
install(TARGETS exp DESTINATION bin)
