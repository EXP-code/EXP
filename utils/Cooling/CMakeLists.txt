message(${CMAKE_CURRENT_SOURCE_DIR})

set(bin_PROGRAMS hctest)

add_executable(hctest hctest.cc)

target_link_libraries(hctest OpenMP::OpenMP_CXX MPI::MPI_CXX
  Boost::program_options expdsmc yaml-cpp exputil EXP
  ${CUDA_LIBRARIES} ${VTK_LIBRARIES})

target_include_directories(hctest
  PUBLIC $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/DSMC/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/>
  ${CMAKE_BINARY_DIR} ${DEP_INC}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${VTK_INCLUDE}
  )

