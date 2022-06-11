
set(DSMC_SRC TreeDSMC.cc Collide.cc CollideLTE.cc HeatCool.cc
  CollideIon.cc Ion.cc interactSelect.cc Elastic.cc TopBase.cc
  NTC.cc Quantile.cc QuantileBag.cc phfit2.f Badnell.cc
  atomic_constants.cc spline.cc cudaCollideIon.cu)

set(TESTDSMC_SRC Ion.cc interactSelect.cc Elastic.cc HeatCool.cc
  TopBase.cc Quantile.cc QuantileBag.cc phfit2.f atomic_constants.cc
  Badnell.cc spline.cc cudaIon.cu)

set(common_INCLUDE_DIRS
  $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/>
  ${CMAKE_BINARY_DIR} ${DEP_INC}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${EIGEN3_INCLUDE_DIR})

set(common_LINKLIB ${DEP_LIB} OpenMP::OpenMP_CXX MPI::MPI_CXX
  Boost::serialization exputil yaml-cpp ${VTK_LIBRARIES})

if(CUDA_FOUND)
  list(APPEND common_INCLUDE_DIRS ${CUDA_TOOLKIT_INCLUDE})
  list(APPEND common_LINKLIB ${CUDA_LIBRARIES} ${CUDA_nvToolsExt_LIBRARY})
endif()

if(SLURM_FOUND)
  list(APPEND common_LINKLIB ${SLURM_LIBRARY})
endif()

if(TIRPC_FOUND)
  list(APPEND common_LINKLIB ${TIRPC_LIBRARY})
endif()

add_library(expdsmc ${DSMC_SRC})
set_target_properties(expdsmc PROPERTIES OUTPUT_NAME expdsmc)
target_include_directories(expdsmc PUBLIC ${common_INCLUDE_DIRS})
target_link_libraries(expdsmc PUBLIC ${common_LINKLIB})
install(TARGETS expdsmc DESTINATION lib)

add_library(testdsmc ${TESTDSMC_SRC})
set_target_properties(testdsmc PROPERTIES OUTPUT_NAME testdsmc)
target_include_directories(testdsmc PUBLIC ${common_INCLUDE_DIRS})
target_link_libraries(testdsmc PUBLIC ${common_LINKLIB})
install(TARGETS testdsmc DESTINATION lib)

set(bin_PROGRAMS testCross testIonRecomb genIonRecomb testQuantile
  testExciteRate testRatio makeRecombCache)

if(ENABLE_CUDA)
  list(APPEND bin_PROGRAMS testCrossCuda)
endif()

foreach(program ${bin_PROGRAMS})
  add_executable(${program} ${program}.cc)
  target_link_libraries(${program} testdsmc ${common_LINKLIB})
  target_include_directories(${program} PUBLIC ${common_INCLUDE_DIRS})
  install(TARGETS ${program} DESTINATION bin)
endforeach()


