
set(USER_MODULES user bar logpot disk halo periodic reflect trigger
  mndisk tidalrad satfix mwgala)

set (common_INCLUDE_DIRS
  $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/>
  ${CMAKE_BINARY_DIR} ${DEP_INC}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..)

set (common_LINKLIBS ${DEP_LIB} OpenMP::OpenMP_CXX MPI::MPI_CXX
  exputil EXPlib yaml-cpp)

if(CUDA_FOUND)
  list(APPEND common_INCLUDE_DIRS ${CUDA_TOOLKIT_INCLUDE})
  list(APPEND common_LINKLIB ${CUDA_LIBRARIES} ${CUDA_nvToolsExt_LIBRARY})
endif()

set(user_SRC      UserTest.cc)
set(bar_SRC       UserBar.cc)
set(satfix_SRC    SatFix.cc)
set(satfixorb_SRC SatFixOrb.cc)
set(mwgala_SRC    UserMW.cc)
set(logpot_SRC    UserLogPot.cc)
set(halo_SRC      UserHalo.cc)
set(disk_SRC      UserDisk.cc)
set(mndisk_SRC    UserMNdisk.cc)
set(tidalrad_SRC  UserTidalRad.cc)
set(periodic_SRC  UserPeriodic.cc)
set(reflect_SRC   UserReflect.cc)
set(trigger_SRC   UserBarTrigger.cc)

if (ENABLE_CUDA)
  set(cudatest_Src UserTestCuda.cc cudaUserTest.cu)
endif()

foreach(mlib ${USER_MODULES})
  add_library(${mlib} ${${mlib}_SRC})
  set_target_properties(${mlib} PROPERTIES OUTPUT_NAME ${mlib})
  target_include_directories(${mlib} PUBLIC ${common_INCLUDE_DIRS})
  target_link_libraries(${mlib} PUBLIC ${common_LINKLIBS})
  install(TARGETS ${mlib} DESTINATION lib/user)
endforeach()

