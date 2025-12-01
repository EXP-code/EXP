set(bin_PROGRAMS nativetoh5 h5compare viewcoefs h5power makecoefs testread)

set(common_LINKLIB OpenMP::OpenMP_CXX MPI::MPI_CXX yaml-cpp exputil
  ${VTK_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})

if(PNG_FOUND)
  list(APPEND common_LINKLIB PNG::PNG)
endif()

set(common_INCLUDE $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/>
  ${CMAKE_BINARY_DIR} ${DEP_INC} ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/.. ${HighFive_SOURCE_DIR}/include
  ${HDF5_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} ${FFTW_INCLUDE_DIRS})


if(ENABLE_CUDA)
  list(APPEND common_LINKLIB CUDA::toolkit CUDA::cudart)
  if (CUDAToolkit_VERSION VERSION_GREATER_EQUAL 12)
    list(APPEND common_LINKLIB CUDA::nvtx3)
  else ()
    list(APPEND common_LINKLIB CUDA::nvToolsExt)
  endif ()
endif()

if(SLURM_FOUND)
  list(APPEND common_LINKLIB ${SLURM_LIBRARY})
endif()

if(ENABLE_XDR AND TIRPC_FOUND)
  list(APPEND common_INCLUDE ${TIRPC_INCLUDE_DIRS})
  list(APPEND common_LINKLIB ${TIRPC_LIBRARIES})
endif()

# Make the expui shared library
#

set(expui_SOURCES BasisFactory.cc BiorthBasis.cc FieldBasis.cc
  Covariance.cc CoefContainer.cc CoefStruct.cc FieldGenerator.cc
  expMSSA.cc Coefficients.cc KMeans.cc Centering.cc
  ParticleIterator.cc Koopman.cc BiorthBess.cc SvdSignChoice.cc
  KDdensity.cc)
add_library(expui ${expui_SOURCES})
set_target_properties(expui PROPERTIES OUTPUT_NAME expui)
target_include_directories(expui PUBLIC ${common_INCLUDE})
target_link_libraries(expui PUBLIC ${common_LINKLIB})

# Add public headers from local directory
target_sources(expui PUBLIC FILE_SET HEADERS
      BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}"
      FILES
	${CMAKE_CURRENT_SOURCE_DIR}/Coefficients.H
	${CMAKE_CURRENT_SOURCE_DIR}/CoefStruct.H
	${CMAKE_CURRENT_SOURCE_DIR}/BasisFactory.H
	${CMAKE_CURRENT_SOURCE_DIR}/BiorthBasis.H
	${CMAKE_CURRENT_SOURCE_DIR}/BiorthBess.H
	${CMAKE_CURRENT_SOURCE_DIR}/FieldGenerator.H
)		     

# Make the exp config header part of the public install
target_sources(expui PUBLIC FILE_SET HEADERS
      BASE_DIRS "${CMAKE_BINARY_DIR}"
      FILES
	${CMAKE_BINARY_DIR}/config_exp.h
)		     

# Add public headers from the main include directory
target_sources(expui PUBLIC FILE_SET HEADERS
      BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/../include"
      FILES
	${CMAKE_CURRENT_SOURCE_DIR}/../include/DiskDensityFunc.H
	${CMAKE_CURRENT_SOURCE_DIR}/../include/ParticleReader.H
	${CMAKE_CURRENT_SOURCE_DIR}/../include/header.H
	${CMAKE_CURRENT_SOURCE_DIR}/../include/Particle.H
	${CMAKE_CURRENT_SOURCE_DIR}/../include/gadget.H
	${CMAKE_CURRENT_SOURCE_DIR}/../include/tipsy.H
	${CMAKE_CURRENT_SOURCE_DIR}/../include/BiorthCube.H
	${CMAKE_CURRENT_SOURCE_DIR}/../include/BiorthCyl.H
	${CMAKE_CURRENT_SOURCE_DIR}/../include/sltableMP2.H
	${CMAKE_CURRENT_SOURCE_DIR}/../include/SLGridMP2.H
	${CMAKE_CURRENT_SOURCE_DIR}/../include/YamlCheck.H
	${CMAKE_CURRENT_SOURCE_DIR}/../include/EmpCylSL.H
	${CMAKE_CURRENT_SOURCE_DIR}/../include/localmpi.H
	${CMAKE_CURRENT_SOURCE_DIR}/../include/exputils.H
	${CMAKE_CURRENT_SOURCE_DIR}/../include/massmodel.H
	${CMAKE_CURRENT_SOURCE_DIR}/../include/interp.H
	${CMAKE_CURRENT_SOURCE_DIR}/../include/orbit.H
	${CMAKE_CURRENT_SOURCE_DIR}/../include/StringTok.H
	${CMAKE_CURRENT_SOURCE_DIR}/../include/OrthoFunction.H
	${CMAKE_CURRENT_SOURCE_DIR}/../include/gaussQ.H
	${CMAKE_CURRENT_SOURCE_DIR}/../include/PseudoAccel.H
	${CMAKE_CURRENT_SOURCE_DIR}/../include/QuadLS.H
	${CMAKE_CURRENT_SOURCE_DIR}/../include/MonotCubicInterpolator.H
	${CMAKE_CURRENT_SOURCE_DIR}/../include/QPDistF.H
	${CMAKE_CURRENT_SOURCE_DIR}/../include/DiskWithHalo.H
	${CMAKE_CURRENT_SOURCE_DIR}/../include/EmpCyl2d.H
	${CMAKE_CURRENT_SOURCE_DIR}/../include/EXPmath.H
	${CMAKE_CURRENT_SOURCE_DIR}/../include/libvars.H
	${CMAKE_CURRENT_SOURCE_DIR}/../include/Timer.H
	${CMAKE_CURRENT_SOURCE_DIR}/../include/coef.H
)

install(TARGETS expui FILE_SET HEADERS DESTINATION include/EXP)
install(TARGETS expui DESTINATION lib)

# Configure and build the test routines
#
if(ENABLE_UTILS)
  add_executable(nativetoh5 coefstoh5.cc)
  add_executable(h5compare  h5compare.cc)
  add_executable(viewcoefs  viewcoefs.cc)
  add_executable(h5power    h5power.cc)
  add_executable(makecoefs  makecoefs.cc)
  add_executable(testread   testread.cc)

  foreach(program ${bin_PROGRAMS})
    target_link_libraries(${program} expui exputil ${common_LINKLIB})
    target_include_directories(${program} PUBLIC ${common_INCLUDE})
    # temporary Apple clang fix
    if(NOT (APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"))
      target_compile_options(${program} PUBLIC ${OpenMP_CXX_FLAGS})
    endif()    
    install(TARGETS ${program} DESTINATION bin)
  endforeach()

endif()
