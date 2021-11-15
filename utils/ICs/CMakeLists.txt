
set(bin_PROGRAMS shrinkics gensph gendisk gsphere pstmod empinfo empdump
  eofbasis makeIon phIon phConv eofcomp testcoefs testcoefs2 testdeval
  forcetest hdf52accel forcetest2)

set(common_LINKLIB OpenMP::OpenMP_CXX MPI::MPI_CXX yaml-cpp EXP
  exputil ${VTK_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})

if(CUDA_FOUND)
  list(APPEND common_LINKLIB ${CUDA_LIBRARIES} ${CUDA_nvToolsExt_LIBRARY})
endif()

set(common_INCLUDE
  $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/DSMC/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/utils/SL/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/>
  ${CMAKE_BINARY_DIR} ${DEP_INC}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${HDF5_INCLUDE_DIRS})

add_executable(shrinkics shrinkics.cc)

add_executable(gensph gensph.cc SphericalSL.cc EllipForce.cc
  legendre.cc SParticle.cc ddplgndr.cc sinecosine.cc)
target_link_libraries(gensph fftw3)

add_executable(gendisk initial.cc DiskHalo.cc SphericalSL.cc
  legendre.cc SParticle.cc ddplgndr.cc sinecosine.cc AddDisk.cc)

add_executable(gsphere bonnerebert.cc)

add_executable(pstmod pst_model.cc ../SL/CylindricalDisk.cc
  EllipsoidForce.cc ../SL/SLSphere.cc)

add_executable(empinfo empinfo.cc)

add_executable(empdump empdump.cc)

add_executable(eofbasis eof_basis.cc)

add_executable(eofcomp eof_compare.cc)

add_executable(makeIon makeIonIC.cc
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/DSMC/atomic_constants.cc>)

add_executable(phIon ph_ion.cc)

add_executable(phConv ph_conv.cc)

add_executable(testcoefs check_coefs.cc SphericalSL.cc legendre.cc
  sinecosine.cc)

add_executable(testcoefs2 check_coefs2.cc DiskEval.cc SphericalSL.cc
  legendre.cc sinecosine.cc)

add_executable(testdeval test_deval.cc DiskEval.cc)

add_executable(forcetest force_test.cc DiskEval.cc)

add_executable(forcetest2 force_test2.cc DiskEval.cc)

add_executable(hdf52accel hdf52accel.cc)

foreach(program ${bin_PROGRAMS})
  target_link_libraries(${program} ${common_LINKLIB})
  target_include_directories(${program} PUBLIC ${common_INCLUDE})
  install(TARGETS ${program} DESTINATION bin)
endforeach()
