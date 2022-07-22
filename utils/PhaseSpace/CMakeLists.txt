
set(bin_PROGRAMS pspinfo pspstat psp2ascii ascii2psp psp2histo
  psp2histo1d snap2histo1d psp2histoE psp2histoH pspreal pspstatS
  psp2range psp2histoT psp2histoTC pspstatT spl2psp ascii2psp2
  psp2rings psp2bess psp2lagu psp2interp diffpsp diffpspN psp2hdf5
  pspmono snap2ascii psporbv)
		
if(HAVE_VTK)
  list(APPEND bin_PROGRAMS psp2eHisto psp2vtu)
endif()

if(HAVE_XDR)
  list(APPEND bin_PROGRAMS tipstd2psp)
endif()

set(common_LINKLIB OpenMP::OpenMP_CXX MPI::MPI_CXX yaml-cpp exputil
  ${VTK_LIBRARIES} ${HDF5_LIBRARIES})

if(CUDA_FOUND)
  list(APPEND common_LINKLIB ${CUDA_LIBRARIES} ${CUDA_nvToolsExt_LIBRARY})
endif()

if(TIRPC_FOUND)
  list(APPEND common_LINKLIB ${TIRPC_LIBRARIES})
endif()

set(common_INCLUDE
  $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/DSMC/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/>
  ${CMAKE_BINARY_DIR} ${DEP_INC}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${HDF5_INCLUDE_DIRS})

add_executable(pspinfo pspinfo.cc PSP.cc)
add_executable(pspstat pspstat.cc PSP.cc)
add_executable(pspstatS pspstatS.cc PSP.cc)
add_executable(pspstatT pspstatT.cc PSP.cc)
add_executable(psp2ascii psp2ascii.cc PSP.cc)
add_executable(snap2ascii snap2ascii.cc)
add_executable(ascii2psp ascii2psp.cc PSP.cc)
add_executable(ascii2psp2 ascii2psp2.cc PSP.cc)
add_executable(pspreal pspreal.cc PSP.cc)
add_executable(psp2histo psp2histo.cc PSP.cc)
add_executable(snap2histo1d snap2histo1d.cc)
add_executable(psp2histo1d psp2histo1d.cc PSP.cc)
add_executable(psp2histoE psp2histoE.cc PSP.cc)
add_executable(psp2histoH psp2histoH.cc PSP.cc)
add_executable(psp2histoT psp2histoT.cc PSP.cc)
add_executable(psp2histoTC psp2histoTC.cc PSP.cc)
add_executable(psp2range psp2range.cc PSP.cc)
add_executable(spl2psp spl2psp.cc PSP.cc)
add_executable(psp2rings psp2rings.cc PSP.cc)
add_executable(psp2bess psp2bess.cc PSP.cc Bess.cc)
add_executable(psp2lagu psp2lagu.cc PSP.cc)
add_executable(psp2interp psp2interp.cc PSP.cc)
add_executable(diffpsp diffpsp.cc MakeModel.cc PSP.cc)
add_executable(pspmono pspmono.cc MakeModel.cc PSP.cc)
add_executable(diffpspN diffpspN.cc PSP.cc)
add_executable(psp2hdf5 psp2hdf5.cc PSP.cc)
add_executable(psporbv psporbv.cc PSP.cc)


if(HAVE_VTK)
  add_executable(psp2eHisto psp2eHisto.cc PSP.cc)
  add_executable(psp2vtu psp2vtu.cc PSP.cc)
endif()

if(HAVE_XDR)
  add_executable(tipstd2psp tipstd2psp.cc PSP.cc tipsyxdr.c)
endif()

foreach(program ${bin_PROGRAMS})
  target_link_libraries(${program} ${common_LINKLIB})
  target_include_directories(${program} PUBLIC ${common_INCLUDE})
  install(TARGETS ${program} DESTINATION bin)
endforeach()

# Reset tipstd2psp targets for TIRPC
if(HAVE_TIRPC)
  target_link_libraries(tipstd2psp ${common_LINKLIB} ${TIRPC_LDFLAGS})
  target_include_directories(tipstd2psp PUBLIC ${common_INCLUDE} ${TIRPC_INCLUDE_DIR})
endif()

