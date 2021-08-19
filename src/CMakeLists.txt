message(${CMAKE_CURRENT_SOURCE_DIR})
set(CUDA_SRC cudaParticle.cu cudaSphericalBasis.cu cudaCylinder.cu
  cudaEmpCylSL.cu cudaComponent.cu NVTX.cc cudaHOT.cu
  cudaIncpos.cu cudaIncvel.cu cudaMultistep.cu cudaOrient.cu)

set(exp_SOURCES Basis.cc Bessel.cc CBrock.cc Component.cc
  CBrockDisk.cc Cube.cc Cylinder.cc EmpCylSL.cc
  ExternalCollection.cc ExternalForce.cc Hernquist.cc Orient.cc
  Particle.cc PotAccel.cc ScatterMFP.cc SphericalBasis.cc
  AxisymmetricBasis.cc Sphere.cc SphericalCoefs.cc
  CylindricalCoefs.cc TwoCenter.cc EJcom.cc global.cc begin.cc
  ddplgndr.cc Direct.cc Shells.cc NoForce.cc end.cc
  OutputContainer.cc OutPS.cc OutPSQ.cc OutPSN.cc OutPSP.cc
  OutCHKPT.cc OutCHKPTQ.cc Output.cc header.cc externalShock.cc
  generateRelaxation.cc HaloBulge.cc incpos.cc incvel.cc
  ComponentContainer.cc OutAscii.cc OutMulti.cc OutRelaxation.cc
  OrbTrace.cc OutDiag.cc OutLog.cc OutCoef.cc multistep.cc
  parse.cc Slab.cc SlabSL.cc step.cc sbessj.cc sbessz.cc
  tidalField.cc ultra.cc ultrasphere.cc MPL.cc zbrent.cc
  OutFrac.cc OutCalbr.cc ParticleFerry.cc pCell.cc chkSlurm.c
  chkTimer.cc pHOT.cc VtkPCA.cc GravKernel.cc ${CUDA_SRC})

# # shared lib
# add_library(exp SHARED ${exp_SOURCES})
# set_target_properties(exp PROPERTIES OUTPUT_NAME exp)
# target_include_directories(
#   exp
#   PUBLIC $<INSTALL_INTERFACE:include>
#   $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
#   $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/>
#   ${CMAKE_BINARY_DIR} ${DEP_INC}
#   ${CMAKE_CURRENT_SOURCE_DIR}
#   ${VTK_INCLUDE})
# target_link_libraries(exp PUBLIC ${DEP_LIB} OpenMP::OpenMP_CXX MPI::MPI_CXX Boost::random yaml-cpp ${CUDA_LIBRARIES})


add_library(exp_static STATIC ${exp_SOURCES})
set_target_properties(exp_static PROPERTIES OUTPUT_NAME exp)
target_include_directories(
  exp_static
  PUBLIC $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/>
  ${CMAKE_BINARY_DIR} ${DEP_INC}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${VTK_INCLUDE})
target_link_libraries(exp_static PUBLIC ${DEP_LIB} OpenMP::OpenMP_CXX MPI::MPI_CXX Boost::random yaml-cpp ${CUDA_LIBRARIES})

add_subdirectory(DSMC)
