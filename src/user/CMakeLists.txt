
set(USER_MODULES user bar logpot disk halo mndisk mwgala)

set (common_INCLUDE_DIRS
  $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/>
  ${CMAKE_BINARY_DIR} ${DEP_INC}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..)

set (common_LINKLIB ${DEP_LIB} OpenMP::OpenMP_CXX MPI::MPI_CXX
  exputil EXPlib yaml-cpp)

set(user_SRC      UserTest.cc)
set(bar_SRC       UserBar.cc)
set(logpot_SRC    UserLogPot.cc)
set(disk_SRC      UserDisk.cc)
set(halo_SRC      UserHalo.cc)
set(mndisk_SRC    UserMNdisk.cc)
set(mwgala_SRC    UserMW.cc)

if(ENABLE_CUDA)
  list(APPEND common_LINKLIB CUDA::toolkit CUDA::cudart)
  if (CUDAToolkit_VERSION VERSION_GREATER_EQUAL 12)
    list(APPEND common_LINKLIB CUDA::nvtx3)
  else ()
    list(APPEND common_LINKLIB  CUDA::nvToolsExt)
  endif ()
  # set_source_files_properties(UserBar.cc UserDisk.cc UserHalo.cc 
  #   UserLogPot.cc UserMNdisk.cc UserMW.cc UserTest.cc UserTestCuda.cc 
  #   PROPERTIES LANGUAGE CUDA)
  set(cudatest_SRC UserTestCuda.cc cudaUserTest.cu)
  list(APPEND bar_SRC cudaUserBar.cu)
endif()

foreach(mlib ${USER_MODULES})
  add_library(${mlib} ${${mlib}_SRC})
  set_target_properties(${mlib} PROPERTIES OUTPUT_NAME ${mlib})
  target_include_directories(${mlib} PUBLIC ${common_INCLUDE_DIRS})
  target_link_libraries(${mlib} PUBLIC ${common_LINKLIB})
  install(TARGETS ${mlib} DESTINATION lib/user)
endforeach()

