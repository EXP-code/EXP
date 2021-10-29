
set(bin_PROGRAMS testConfig simpleConfig readSL testBarrier)

set(common_LINKLIB OpenMP::OpenMP_CXX MPI::MPI_CXX
  Boost::program_options yaml-cpp exputil EXP ${VTK_LIBRARIES})

set(common_INCLUDE
  $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/>
  ${CMAKE_BINARY_DIR} ${DEP_INC}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${VTK_INCLUDE}
  )

if (ENABLE_CUDA)
  list(APPEND common_LINKLIB ${CUDA_LIBRARIES} ${CUDA_nvToolsExt_LIBRARY} )
endif()

add_executable(testConfig     test_config.cc)
add_executable(simpleConfig   simple_config.cc)
add_executable(readSL         read_slsph.cc)
add_executable(testBarrier    test_barrier.cc)

foreach(program ${bin_PROGRAMS})
  target_link_libraries(${program} ${common_LINKLIB})
  target_include_directories(${program} PUBLIC ${common_INCLUDE})
  # install(TARGETS ${program} DESTINATION bin)
endforeach()