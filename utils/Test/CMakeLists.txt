
set(bin_PROGRAMS testBarrier expyaml orthoTest)

set(common_LINKLIB OpenMP::OpenMP_CXX MPI::MPI_CXX expui exputil
  yaml-cpp ${VTK_LIBRARIES})

if(PNG_FOUND)
  list(APPEND common_LINKLIB PNG::PNG)
endif()

set(common_INCLUDE
  $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/>
  ${CMAKE_BINARY_DIR} ${DEP_INC}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..)

if(ENABLE_CUDA)
  list(APPEND common_LINKLIB CUDA::toolkit CUDA::cudart)
  if (CUDAToolkit_VERSION VERSION_GREATER_EQUAL 12)
    list(APPEND common_LINKLIB CUDA::nvtx3)
  else ()
    list(APPEND common_LINKLIB  CUDA::nvToolsExt)
  endif ()
endif()

if(ENABLE_XDR AND TIRPC_FOUND)
  list(APPEND common_LINKLIB ${TIRPC_LIBRARIES})
endif()

add_executable(testBarrier    test_barrier.cc)
add_executable(expyaml        test_config.cc)
add_executable(orthoTest      orthoTest.cc Biorth2Ortho.cc)

foreach(program ${bin_PROGRAMS})
  target_link_libraries(${program} ${common_LINKLIB})
  target_include_directories(${program} PUBLIC ${common_INCLUDE})
  # Don't install these tests
  # install(TARGETS ${program} DESTINATION bin)
endforeach()

install(TARGETS expyaml DESTINATION bin)
