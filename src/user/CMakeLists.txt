message(${CMAKE_CURRENT_SOURCE_DIR})

set(USER_MODULES user bar ebar logpot disk ebarp halo ebars ebarn
  periodic atmos reflect trigger mndisk tidalrad massevo addmass )

if (USER_ALL)
  list(APPEND USER_MODULES diffuse torque sat satfix diffrot respot
    respotn respotorb rptest wake atmos satwake pst userslab shear
    snheat profile)
endif()


set (PUBLIC_INCLUDE_LIST
  $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/>
  ${CMAKE_BINARY_DIR} ${DEP_INC}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${VTK_INCLUDE})

set (PUBLIC_LINK_LIST ${DEP_LIB} OpenMP::OpenMP_CXX MPI::MPI_CXX
Boost::random Boost::serialization Boost::filesystem exputil
yaml-cpp ${CUDA_LIBRARIES})

set (common_SRC ResPot.cc ResPotOrb.cc BarForcing.cc Perturbation.cc
  EllipForce.cc SatelliteOrbit.cc CircularOrbit.cc
  TwoBodyDiffuse.cc TimeSeriesCoefs.cc LinearOrbit.cc
  UnboundOrbit.cc Synchronize.cc
  )

set (swake_SRC respmat6.cc cratint.cc rombeveccmp.cc)

add_library(common ${common_SRC})
target_include_directories(common PUBLIC ${PUBLIC_INCLUDE_LIST})
target_link_libraries(common PUBLIC ${PUBLIC_LINK_LIST})

if (USER_ALL)
  add_library(swake ${common_SRC})
  target_include_directories(swake PUBLIC ${PUBLIC_INCLUDE_LIST})
  target_link_libraries(swake PUBLIC ${PUBLIC_LINK_LIST})
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
set(slab_SRC      UserSlabHalo.cc)
set(trigger_SRC   UserBarTrigger.cc)
set(shear_SRC     UserShear.cc)
set(profile_SRC   UserProfile.cc)

if (ENABLE_CUDA)
  list(APPEND sat_SRC cudaUserSat.cu)
  set(cudatest_Src UserTestCuda.cc cudaUserTest.cu)
endif()

foreach(mlib ${USER_MODULES})
  add_library(${mlib} ${${mlib}_SRC})
  set_target_properties(${mlib} PROPERTIES OUTPUT_NAME ${mlib})
  target_include_directories(${mlib} PUBLIC ${PUBLIC_INCLUDE_LIST})
  target_link_libraries(${mlib} PUBLIC common ${PUBLIC_LINK_LIST})
  install(TARGETS ${mlib} DESTINATION lib/user)
endforeach()

if (USER_ALL)
  target_link_libraries(satwake PUBLIC common swake ${PUBLIC_LINK_LIST})
endif()

