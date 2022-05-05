set(bin_PROGRAMS exp_haloN disk_noise halo_noise expmssa)

set(common_LINKLIB
yaml-cpp PNG::PNG ${VTK_LIBRARIES} ${EIGEN3_LIBRARIES}
${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} ${FFTW_DOUBLE_LIB})

set(common_INCLUDE $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/> ${CMAKE_BINARY_DIR}
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/extern/yaml-cpp/include>
  ${DEP_INC} ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/.. ${EIGEN3_INCLUDE_DIR} ${FFTW_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS})

add_executable(exp_haloN     exp_haloN.cc TransformFFT.cc Coefs.cc)
add_executable(disk_noise    exp_disk_noise.cc Coefs.cc)
add_executable(halo_noise    exp_halo_noise.cc Coefs.cc)
add_executable(expmssa       expmssa.cc TransformFFT.cc Coefs.cc CoefDB.cc KMeans.cc)

foreach(program ${bin_PROGRAMS})
  target_link_libraries(${program} ${common_LINKLIB} OpenMP::OpenMP_CXX)
  target_include_directories(${program} PUBLIC ${common_INCLUDE})
  install(TARGETS ${program} DESTINATION bin)
endforeach()

