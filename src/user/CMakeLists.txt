
set(USER_MODULES user bar logpot disk halo mndisk mwgala)

set (common_INCLUDE_DIRS
  $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/>
  ${CMAKE_BINARY_DIR} ${DEP_INC}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..)

set (common_LINKLIBS ${DEP_LIB} OpenMP::OpenMP_CXX MPI::MPI_CXX
  exputil EXPlib yaml-cpp)

if(ENABLE_CUDA)
  list(APPEND common_INCLUDE_DIRS ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} ${CUDAToolkit_INCLUDE_DIRS})
  list(APPEND common_LINKLIBS CUDA::cudart CUDA::nvToolsExt)
endif()

set(user_SRC      UserTest.cc)
set(bar_SRC       UserBar.cc)
set(logpot_SRC    UserLogPot.cc)
set(disk_SRC      UserDisk.cc)
set(halo_SRC      UserHalo.cc)
set(mndisk_SRC    UserMNdisk.cc)
set(mwgala_SRC    UserMW.cc)

if(ENABLE_CUDA)
  set(cudatest_SRC UserTestCuda.cc cudaUserTest.cu)
endif()

foreach(mlib ${USER_MODULES})
  add_library(${mlib} ${${mlib}_SRC})
  set_target_properties(${mlib} PROPERTIES OUTPUT_NAME ${mlib})
  target_include_directories(${mlib} PUBLIC ${common_INCLUDE_DIRS})
  target_link_libraries(${mlib} PUBLIC ${common_LINKLIBS})
  install(TARGETS ${mlib} DESTINATION lib/user)
endforeach()

