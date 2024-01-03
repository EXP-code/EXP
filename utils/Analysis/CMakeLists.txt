find_package (Python3 COMPONENTS Interpreter Development)

set(bin_PROGRAMS diskcoef halocoef diskcoef diskprof haloprof slabprof
  psphisto gas2d makegrid pspbox eofpeek diskfreqs mssaprofD mssaprofH
  pcatest readcoefsD readcoefsH sphprof pspxvalH pspxvalH2 pspxvalD
  pspxvalD2 pspKLH pspKLD diskeof yaml_diff kdtest)

set(common_LINKLIB OpenMP::OpenMP_CXX MPI::MPI_CXX yaml-cpp exputil
  ${VTK_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})

if(PNG_FOUND)
  list(APPEND common_LINKLIB PNG::PNG)
endif()

if(CUDA_FOUND)
  list(APPEND common_LINKLIB ${CUDA_CUDART_LIBRARY} ${CUDA_nvToolsExt_LIBRARY})
endif()

if(SLURM_FOUND)
  list(APPEND common_LINKLIB ${SLURM_LIBRARY})
endif()

if(ENABLE_XDR AND TIRPC_FOUND)
  list(APPEND common_LINKLIB ${TIRPC_LIBRARIES})
endif()

set(common_INCLUDE $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/>
  ${CMAKE_BINARY_DIR} ${DEP_INC} ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/.. ${HDF5_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})

add_executable(diskcoef    diskprof_coef.cc)
add_executable(halocoef    haloprof_coef.cc SphSL.cc)
add_executable(diskprof    diskprof.cc)	    
add_executable(haloprof    haloprof.cc SphSL.cc)
add_executable(slabprof    slabprof.cc)
add_executable(psphisto    psphisto.cc)
add_executable(gas2d       gas2dcyl.cc)
add_executable(makegrid    makegrid.cc)
add_executable(pspbox      pspbox.cc)
add_executable(eofpeek     eofpeek.cc)
add_executable(diskfreqs   diskfreqs.cc SphSL.cc)
add_executable(mssaprofD   mssaprof_disk.cc)
add_executable(mssaprofH   mssaprof_halo.cc SphSL.cc)
add_executable(pcatest     pcatest.cc)
add_executable(readcoefsD  readcoefsD.cc)
add_executable(readcoefsH  readcoefsH.cc)
add_executable(sphprof     sphprof.cc SphSL.cc)
add_executable(pspxvalH    cross_validation_sph.cc SphSL.cc)
add_executable(pspxvalH2   cross_validation_sph2.cc SphSL.cc)
add_executable(pspxvalD    cross_validation_cyl.cc dtorh1.f rout.f)
add_executable(pspxvalD2   cross_validation_cyl2.cc)
add_executable(pspKLH      KL_sph.cc SphSL.cc)
add_executable(pspKLD      KL_cyl.cc)
add_executable(diskeof     diskeof.cc)
add_executable(yaml_diff   yaml_diff.cc)
add_executable(kdtest      KDcyltest.cc)


foreach(program ${bin_PROGRAMS})
  target_link_libraries(${program} ${common_LINKLIB})
  target_include_directories(${program} PUBLIC ${common_INCLUDE})
  install(TARGETS ${program} DESTINATION bin)
endforeach()
