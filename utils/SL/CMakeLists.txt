
set(bin_PROGRAMS slcheck slshift orthochk diskpot qtest eoftest
  oftest slabchk)
		
set(common_LINKLIB OpenMP::OpenMP_CXX MPI::MPI_CXX yaml-cpp exputil
  ${VTK_LIBRARIES})

if(ENABLE_CUDA)
  list(APPEND common_LINKLIB CUDA::cudart)
  if (CUDAToolkit_VERSION VERSION_GREATER_EQUAL 12)
    list(APPEND common_LINKLIB CUDA::nvtx3)
  else ()
    list(APPEND common_LINKLIB  CUDA::nvToolsExt)
  endif ()
endif()

if(ENABLE_XDR AND TIRPC_FOUND)
  list(APPEND common_LINKLIB ${TIRPC_LIBRARIES})
endif()

set(common_INCLUDE
  $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/DSMC/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/>
  ${CMAKE_BINARY_DIR} ${DEP_INC}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..)

add_executable(slcheck slcheck.cc)
add_executable(slabchk slabchk.cc Model1d.cc)
add_executable(slshift slshift.cc SLSphere.cc)
add_executable(orthochk orthochk.cc)
add_executable(diskpot diskpot.cc CylindricalDisk.cc SLSphere.cc)
add_executable(qtest qtest.cc)
add_executable(eoftest EOF2d.cc)
add_executable(oftest oftest.cc)

foreach(program ${bin_PROGRAMS})
  target_link_libraries(${program} ${common_LINKLIB})
  target_include_directories(${program} PUBLIC ${common_INCLUDE})
  install(TARGETS ${program} DESTINATION bin)
endforeach()
