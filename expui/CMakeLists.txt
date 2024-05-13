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

# Make the expui shared library
#
set(expui_SOURCES BasisFactory.cc BiorthBasis.cc FieldBasis.cc
  CoefContainer.cc CoefStruct.cc FieldGenerator.cc expMSSA.cc
  Coefficients.cc KMeans.cc Centering.cc ParticleIterator.cc
  Koopman.cc BiorthBess.cc)
add_library(expui ${expui_SOURCES})
set_target_properties(expui PROPERTIES OUTPUT_NAME expui)
target_include_directories(expui PUBLIC ${common_INCLUDE})
target_link_libraries(expui PUBLIC ${common_LINKLIB})

install(TARGETS expui DESTINATION lib)

# Configure and build the test routines
#
add_executable(nativetoh5 coefstoh5.cc)
add_executable(h5compare  h5compare.cc)
add_executable(viewcoefs  viewcoefs.cc)
add_executable(h5power    h5power.cc)
add_executable(makecoefs  makecoefs.cc)
add_executable(testread   testread.cc)

foreach(program ${bin_PROGRAMS})
  target_link_libraries(${program} expui exputil ${common_LINKLIB})
  target_include_directories(${program} PUBLIC ${common_INCLUDE})
  target_compile_options(${program} PUBLIC ${OpenMP_CXX_FLAGS})
  install(TARGETS ${program} DESTINATION bin)
endforeach()
