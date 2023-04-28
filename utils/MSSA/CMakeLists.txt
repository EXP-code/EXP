set(bin_PROGRAMS exp_haloN disk_noise halo_noise expmssa)

set(common_LINKLIB
  yaml-cpp PNG::PNG exputil ${VTK_LIBRARIES} ${EIGEN3_LIBRARIES}
  ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} ${FFTW_DOUBLE_LIB})

set(common_INCLUDE $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/> ${CMAKE_BINARY_DIR}
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/extern/yaml-cpp/include>
  ${DEP_INC} ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/.. ${EIGEN3_INCLUDE_DIR} ${FFTW_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS})

if(CUDA_FOUND)
  list(APPEND common_LINKLIB ${CUDA_LIBRARIES} ${CUDA_nvToolsExt_LIBRARY})
endif()

if(SLURM_FOUND)
  list(APPEND common_LINKLIB ${SLURM_LIBRARY})
endif()

if(ENABLE_XDR AND TIRPC_FOUND)
  list(APPEND common_LINKLIB ${TIRPC_LIBRARIES})
endif()

add_executable(exp_haloN     exp_haloN.cc Coefs.cc)
add_executable(disk_noise    exp_disk_noise.cc Coefs.cc)
add_executable(halo_noise    exp_halo_noise.cc Coefs.cc)
add_executable(expmssa       expmssa.cc Coefs.cc CoefDB.cc KMeans.cc)

foreach(program ${bin_PROGRAMS})
  target_link_libraries(${program} ${common_LINKLIB} OpenMP::OpenMP_CXX MPI::MPI_CXX)
  target_include_directories(${program} PUBLIC ${common_INCLUDE})
  install(TARGETS ${program} DESTINATION bin)
endforeach()

