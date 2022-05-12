
set(bin_PROGRAMS readSL testBarrier testConfig)

set(common_LINKLIB OpenMP::OpenMP_CXX MPI::MPI_CXX exputil yaml-cpp
  ${VTK_LIBRARIES})

set(common_INCLUDE
  $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/>
  ${CMAKE_BINARY_DIR} ${DEP_INC}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..)

if (ENABLE_CUDA)
  list(APPEND common_LINKLIB ${CUDA_LIBRARIES} ${CUDA_nvToolsExt_LIBRARY} )
endif()

add_executable(readSL         read_slsph.cc)
add_executable(testBarrier    test_barrier.cc)
add_executable(testConfig     test_config.cc)

foreach(program ${bin_PROGRAMS})
  target_link_libraries(${program} ${common_LINKLIB})
  target_include_directories(${program} PUBLIC ${common_INCLUDE})
  # Don't install these tests
  # install(TARGETS ${program} DESTINATION bin)
endforeach()
