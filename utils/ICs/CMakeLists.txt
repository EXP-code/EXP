
set(bin_PROGRAMS shrinkics gensph gendisk gendisk2d gsphere pstmod
  empinfo empdump eofbasis eofcomp testcoefs testcoefs2 testdeval
  forcetest hdf52accel forcetest2 cylcache modelfit test2d addsphmod
  cubeics zangics slabics)

set(common_LINKLIB OpenMP::OpenMP_CXX MPI::MPI_CXX yaml-cpp exputil
  ${VTK_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES}
  ${HDF5_CXX_LIBRARIES})

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
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/extern/DSMC/src>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/utils/SL/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/>
  ${CMAKE_BINARY_DIR} ${DEP_INC}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${HDF5_INCLUDE_DIRS})

add_executable(shrinkics shrinkics.cc)

add_executable(gensph gensph.cc SphericalSL.cc EllipForce.cc
  legendre.cc SParticle.cc ddplgndr.cc sinecosine.cc)
target_link_libraries(gensph ${FFTW_DOUBLE_LIB})
target_include_directories(gensph PUBLIC ${FFTW_INCLUDE_DIRS})

add_executable(gendisk initial.cc DiskHalo.cc SphericalSL.cc
  legendre.cc SParticle.cc ddplgndr.cc sinecosine.cc AddDisk.cc)

add_executable(gendisk2d initial2d.cc Disk2d.cc Disk2dHalo.cc
  SphericalSL.cc legendre.cc SParticle.cc ddplgndr.cc sinecosine.cc
  AddDisk.cc)

add_executable(gsphere bonnerebert.cc)

add_executable(pstmod pst_model.cc ../SL/CylindricalDisk.cc
  EllipsoidForce.cc ../SL/SLSphere.cc)

add_executable(empinfo empinfo.cc)

add_executable(empdump empdump.cc)

add_executable(eofbasis eof_basis.cc)

add_executable(eofcomp eof_compare.cc)

add_executable(testcoefs check_coefs.cc SphericalSL.cc legendre.cc
  sinecosine.cc)

add_executable(testcoefs2 check_coefs2.cc DiskEval.cc SphericalSL.cc
  legendre.cc sinecosine.cc)

add_executable(testdeval test_deval.cc DiskEval.cc)

add_executable(forcetest force_test.cc DiskEval.cc)

add_executable(forcetest2 force_test2.cc DiskEval.cc)

add_executable(cylcache cylcache.cc)

add_executable(hdf52accel hdf52accel.cc)

add_executable(modelfit ModelFit.cc)

add_executable(test2d test2d.cc Disk2d.cc)

add_executable(addsphmod addsphmod.cc AddSpheres.cc)

add_executable(cubeics cubeICs.cc)

add_executable(zangics ZangICs.cc)

add_executable(slabics genslab.cc massmodel1d.cc)

foreach(program ${bin_PROGRAMS})
  target_link_libraries(${program} ${common_LINKLIB})
  target_include_directories(${program} PUBLIC ${common_INCLUDE})
  install(TARGETS ${program} DESTINATION bin)
endforeach()
