message(${CMAKE_CURRENT_SOURCE_DIR})

set(DSMC_SRC TreeDSMC.cc Collide.cc CollideLTE.cc HeatCool.cc
  CollideIon.cc Ion.cc interactSelect.cc Elastic.cc TopBase.cc
  NTC.cc Quantile.cc QuantileBag.cc phfit2.f Badnell.cc
  atomic_constants.cc spline.cc cudaCollideIon.cu)

set(TESTDSMC_SRC Ion.cc interactSelect.cc Elastic.cc TopBase.cc
  Quantile.cc QuantileBag.cc phfit2.f atomic_constants.cc
  Badnell.cc spline.cc cudaIon.cu)

set(COMMON_INCLUDE_DIRS
  $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/>
  ${CMAKE_BINARY_DIR} ${DEP_INC}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${EIGEN3_INCLUDE_DIR}
  ${CUDA_TOOLKIT_INCLUDE}
  ${VTK_INCLUDE})

set(COMMON_LINK_LIBS ${DEP_LIB} OpenMP::OpenMP_CXX MPI::MPI_CXX
  Boost::random Boost::serialization Boost::filesystem
  Boost::program_options exputil yaml-cpp ${CUDA_LIBRARIES}
  ${VTK_LIBRARIES})

add_library(expdsmc ${DSMC_SRC})
set_target_properties(expdsmc PROPERTIES OUTPUT_NAME expdsmc)
target_include_directories(expdsmc PUBLIC ${COMMON_INCLUDE_DIRS})
target_link_libraries(expdsmc PUBLIC ${COMMON_LINK_LIBS})

add_library(testdsmc ${TESTDSMC_SRC})
set_target_properties(testdsmc PROPERTIES OUTPUT_NAME testdsmc)
target_include_directories(testdsmc PUBLIC ${COMMON_INCLUDE_DIRS})
target_link_libraries(testdsmc PUBLIC ${COMMON_LINK_LIBS})

set(bin_PROGRAMS testCross testIonRecomb genIonRecomb testQuantile
  testExciteRate testRatio makeRecombCache)

if(ENABLE_CUDA)
  list(APPEND bin_PROGRAMS testCrossCuda)
endif()

foreach(program ${bin_PROGRAMS})
  add_executable(${program} ${program}.cc)
  target_link_libraries(${program} testdsmc ${COMMON_LINK_LIBS})
  target_include_directories(${program} PUBLIC ${COMMON_INCLUDE_DIRS})
endforeach()

install(TARGETS expdsmc DESTINATION lib)

