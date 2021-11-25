
set(CUDA_SRC)
if (ENABLE_CUDA)
  list(APPEND CUDA_SRC 
    cudaParticle.cu cudaSphericalBasis.cu cudaCylinder.cu
    cudaEmpCylSL.cu cudaComponent.cu NVTX.cc cudaHOT.cu
    cudaIncpos.cu cudaIncvel.cu cudaMultistep.cu cudaOrient.cu)
endif()

set(exp_SOURCES Basis.cc Bessel.cc CBrock.cc Component.cc
  CBrockDisk.cc Cube.cc Cylinder.cc EmpCylSL.cc ExternalCollection.cc
  ExternalForce.cc Hernquist.cc Orient.cc Particle.cc PotAccel.cc
  ScatterMFP.cc SphericalBasis.cc AxisymmetricBasis.cc Sphere.cc
  SphericalCoefs.cc CylindricalCoefs.cc TwoCenter.cc EJcom.cc
  global.cc begin.cc ddplgndr.cc Direct.cc Shells.cc NoForce.cc end.cc
  OutputContainer.cc OutPS.cc OutPSQ.cc OutPSN.cc OutPSP.cc
  OutCHKPT.cc OutCHKPTQ.cc Output.cc header.cc externalShock.cc
  generateRelaxation.cc HaloBulge.cc incpos.cc incvel.cc
  ComponentContainer.cc OutAscii.cc OutMulti.cc OutRelaxation.cc
  OrbTrace.cc OutDiag.cc OutLog.cc OutCoef.cc multistep.cc parse.cc
  Slab.cc SlabSL.cc step.cc sbessj.cc sbessz.cc tidalField.cc ultra.cc
  ultrasphere.cc MPL.cc zbrent.cc OutFrac.cc OutCalbr.cc
  ParticleFerry.cc pCell.cc chkSlurm.c chkTimer.cc pHOT.cc
  GravKernel.cc ${CUDA_SRC})

if(HAVE_VTK)
  list(APPEND exp_SOURCES VtkPCA.cc)
endif()

set(common_INCLUDE_DIRS 
  $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/DSMC/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/extern/yaml-cpp/include>
  ${CMAKE_BINARY_DIR} ${DEP_INC}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${EIGEN3_INCLUDE_DIR})

set(DSMC_LIBS)
if (ENABLE_DSMC)
  add_subdirectory(DSMC)
  list(APPEND DSMC_LIBS expdsmc)
endif()

set(common_LINKLIB ${DSMC_LIBS} exputil OpenMP::OpenMP_CXX
  MPI::MPI_CXX yaml-cpp ${VTK_LIBRARIES})

if(CUDA_FOUND)
  list(APPEND common_INCLUDE_DIRS ${CUDA_TOOLKIT_INCLUDE})
  list(APPEND common_LINKLIB ${CUDA_LIBRARIES} ${CUDA_nvToolsExt_LIBRARY})
endif()

if(SLURM_FOUND)
  list(APPEND common_LINKLIB ${SLURM_LIBRARY})
endif()

add_library(EXP ${exp_SOURCES})
set_target_properties(EXP PROPERTIES OUTPUT_NAME EXP)
target_include_directories(EXP PUBLIC ${common_INCLUDE_DIRS})
target_link_libraries(EXP PUBLIC ${common_LINKLIB})

if (ENABLE_USER OR ENABLE_USER_ALL)
  add_subdirectory(user)
endif()

add_executable(exp expand.cc)
target_include_directories(exp PUBLIC ${common_INCLUDE_DIRS})
target_link_libraries(exp PUBLIC ${common_LINKLIB} EXP)

install(TARGETS EXP DESTINATION lib)
install(TARGETS exp DESTINATION bin)
