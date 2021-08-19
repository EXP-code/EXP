
set(DSMC_SRC TreeDSMC.cc Collide.cc CollideLTE.cc HeatCool.cc
  CollideIon.cc Ion.cc interactSelect.cc Elastic.cc TopBase.cc
  NTC.cc Quantile.cc QuantileBag.cc phfit2.f Badnell.cc
  atomic_constants.cc spline.cc cudaCollideIon.cu)

set(TESTDSMC_SRC Ion.cc interactSelect.cc Elastic.cc TopBase.cc
  Quantile.cc QuantileBag.cc phfit2.f atomic_constants.cc
  Badnell.cc spline.cc cudaIon.cu)

add_library(expdsmc STATIC ${DSMC_SRC})
set_target_properties(expdsmc PROPERTIES OUTPUT_NAME expdsmc)
target_include_directories(
  expdsmc
  PUBLIC $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/>
  ${CMAKE_BINARY_DIR} ${DEP_INC}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${VTK_INCLUDE})
target_link_libraries(expdsmc PUBLIC ${DEP_LIB} OpenMP::OpenMP_CXX MPI::MPI_CXX Boost::random Boost::serialization Boost::filesystem exputil_static yaml-cpp ${CUDA_LIBRARIES})

add_library(testdsmc STATIC ${TESTDSMC_SRC})
set_target_properties(testdsmc PROPERTIES OUTPUT_NAME testdsmc)
target_include_directories(
  testdsmc
  PUBLIC $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/>
  ${CMAKE_BINARY_DIR} ${DEP_INC}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${VTK_INCLUDE})
target_link_libraries(testdsmc PUBLIC ${DEP_LIB} OpenMP::OpenMP_CXX MPI::MPI_CXX Boost::random Boost::serialization Boost::filesystem exputil_static yaml-cpp ${CUDA_LIBRARIES})

set(bin_PROGRAMS testCross testIonRecomb genIonRecomb testQuantile testExciteRate testRatio makeRecombCache)
if(ENABLE_CUDA)
  list(APPEND bin_PROGRAMS testCrossCuda)
endif()

foreach(program ${bin_PROGRAMS})
  add_executable(${program} ${program}.cc)

  target_link_libraries(${program} ${DEP_LIB} OpenMP::OpenMP_CXX MPI::MPI_CXX Boost::program_options testdsmc yaml-cpp exputil_static ${CUDA_LIBRARIES} ${VTK_LIBRARIES})

  target_include_directories(${program}
    PUBLIC $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/>
    ${CMAKE_BINARY_DIR} ${DEP_INC}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${VTK_INCLUDE}
    )

endforeach()
