set(bin_PROGRAMS exp_haloN disk_noise halo_noise expmssa)

set(common_LINKLIB yaml-cpp exputil ${VTK_LIBRARIES}
  ${EIGEN3_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES}
  ${FFTW_DOUBLE_LIB} ${FFTW3_LIBRARIES} ${PNG_LIBRARIES})

if(PNG_FOUND)
  list(APPEND common_LINKLIB PNG::PNG)
endif()

if(ENABLE_CUDA)
  list(APPEND common_LINKLIB CUDA::cudart CUDA::nvToolsExt)
endif()

if(SLURM_FOUND)
  list(APPEND common_LINKLIB ${SLURM_LIBRARY})
endif()

if(ENABLE_XDR AND TIRPC_FOUND)
  list(APPEND common_LINKLIB ${TIRPC_LIBRARIES})
endif()

set(common_INCLUDE $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/>
  ${CMAKE_BINARY_DIR} ${DEP_INC} ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/.. ${HDF5_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR}
  ${YAML_CPP_INCLUDE_DIR} ${FFTW_INCLUDE_DIRS} ${FFTW3_INCLUDE_DIRS} ${PNG_INCLUDE_DIRS})

add_executable(exp_haloN     exp_haloN.cc Coefs.cc)
add_executable(disk_noise    exp_disk_noise.cc Coefs.cc)
add_executable(halo_noise    exp_halo_noise.cc Coefs.cc)
add_executable(expmssa       expmssa.cc Coefs.cc CoefDB.cc KMeans.cc)

foreach(program ${bin_PROGRAMS})
  target_link_libraries(${program} ${common_LINKLIB} OpenMP::OpenMP_CXX MPI::MPI_CXX)
  target_include_directories(${program} BEFORE PUBLIC ${common_INCLUDE})
  install(TARGETS ${program} DESTINATION bin)
endforeach()

