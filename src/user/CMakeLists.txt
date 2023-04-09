
set(USER_MODULES user bar ebar logpot disk ebarp halo ebars ebarn
  periodic reflect trigger mndisk tidalrad massevo addmass sat
  satfix mwgala agnnoise)

if (ENABLE_USER_ALL)
  list(APPEND USER_MODULES diffuse torque diffrot respot respotn
    respotorb rptest wake atmos satwake userslab shear snheat
    profile)
endif()


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

set (common_SRC ResPot.cc ResPotOrb.cc BarForcing.cc Perturbation.cc
  EllipForce.cc SatelliteOrbit.cc CircularOrbit.cc
  TwoBodyDiffuse.cc TimeSeriesCoefs.cc LinearOrbit.cc
  UnboundOrbit.cc Synchronize.cc
  )

set (swake_SRC RespMat.cc cratint.cc rombeveccmp.cc)

add_library(expuser_common ${common_SRC})
target_include_directories(expuser_common PUBLIC ${common_INCLUDE_DIRS})
target_link_libraries(expuser_common PUBLIC ${common_LINKLIBS})
install(TARGETS expuser_common DESTINATION lib)

if(ENABLE_USER_ALL)
  add_library(expuser_swake ${swake_SRC})
  target_include_directories(expuser_swake PUBLIC ${common_INCLUDE_DIRS})
  target_link_libraries(expuser_swake PUBLIC ${common_LINKLIBS})
  install(TARGETS expuser_swake DESTINATION lib)
endif()

set(user_SRC      UserTest.cc)
set(bar_SRC       UserBar.cc)
set(ebar_SRC      UserEBar.cc)
set(ebarn_SRC     UserEBarN.cc)
set(ebarp_SRC     UserEBarP.cc)
set(ebars_SRC     UserEBarS.cc)
set(sat_SRC       UserSat.cc)
set(diffuse_SRC   UserDiffuse.cc)
set(satfix_SRC    SatFix.cc)
set(mwgala_SRC    UserMW.cc)
set(satfixorb_SRC SatFixOrb.cc)
set(diffrot_SRC   UserDiffRot.cc)
set(logpot_SRC    UserLogPot.cc)
set(disk_SRC      UserDisk.cc)
set(disk_SRC      UserDisk.cc)
set(mndisk_SRC    UserMNdisk.cc)
set(tidalrad_SRC  UserTidalRad.cc)
set(massevo_SRC   UserMassEvo.cc)
set(addmass_SRC   UserAddMass.cc)
set(halo_SRC      UserHalo.cc)
set(torque_SRC    UserTorque.cc)
set(respot_SRC    UserResPot.cc)
set(respotn_SRC   UserResPotN.cc)
set(respotorb_SRC UserResPotOrb.cc)
set(rptest_SRC    UserRPtest.cc)
set(periodic_SRC  UserPeriodic.cc)
set(snheat_SRC    UserSNheat.cc)
set(atmos_SRC     UserAtmos.cc)
set(reflect_SRC   UserReflect.cc)
set(satwake_SRC   UserSatWake.cc)
set(userslab_SRC  UserSlabHalo.cc)
set(trigger_SRC   UserBarTrigger.cc)
set(shear_SRC     UserShear.cc)
set(profile_SRC   UserProfile.cc)
set(wake_SRC      UserWake.cc)
set(agnnoise_SRC  UserAgnNoise.cc)

if (ENABLE_CUDA)
  list(APPEND sat_SRC cudaUserSat.cu)
  set(cudatest_Src UserTestCuda.cc cudaUserTest.cu)
endif()

foreach(mlib ${USER_MODULES})
  add_library(${mlib} ${${mlib}_SRC})
  set_target_properties(${mlib} PROPERTIES OUTPUT_NAME ${mlib})
  target_include_directories(${mlib} PUBLIC ${common_INCLUDE_DIRS})
  target_link_libraries(${mlib} PUBLIC expuser_common ${common_LINKLIBS})
  install(TARGETS ${mlib} DESTINATION lib/user)
endforeach()

if(ENABLE_USER_ALL)
  target_link_libraries(satwake PUBLIC expuser_common expuser_swake
    ${common_LINKLIBS})
endif()
