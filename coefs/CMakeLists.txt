find_package (Python3 COMPONENTS Interpreter Development)

set(bin_PROGRAMS nativetoh5 h5compare viewcoefs h5power makecoefs)

set(common_LINKLIB OpenMP::OpenMP_CXX MPI::MPI_CXX PNG::PNG yaml-cpp
  exputil ${VTK_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})

set(common_INCLUDE $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/>
  ${CMAKE_BINARY_DIR} ${DEP_INC} ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${HighFive_SOURCE_DIR}/include ${HDF5_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})

if(CUDA_FOUND)
  list(APPEND common_INCLUDE ${CUDA_TOOLKIT_INCLUDE})
  list(APPEND common_LINKLIB ${CUDA_LIBRARIES} ${CUDA_nvToolsExt_LIBRARY})
endif()

if(SLURM_FOUND)
  list(APPEND common_LINKLIB ${SLURM_LIBRARY})
endif()

if(TIRPC_FOUND)
  list(APPEND common_INCLUDE ${TIRPC_INCLUDE_DIRS})
  list(APPEND common_LINKLIB ${TIRPC_LIBRARIES})
endif()

# shared library

set(expcoef_SOURCES BasisFactory.cc CoefFactory.cc CoefStruct.cc
  FieldGenerator.cc)
add_library(expcoef ${expcoef_SOURCES})
set_target_properties(expcoef PROPERTIES OUTPUT_NAME expcoef)
target_include_directories(expcoef PUBLIC ${common_INCLUDE})
target_link_libraries(expcoef PUBLIC ${common_LINKLIBS})

install(TARGETS expcoef DESTINATION lib)

# Test routines

add_executable(nativetoh5 coefstoh5.cc)
add_executable(h5compare  h5compare.cc)
add_executable(viewcoefs  viewcoefs.cc)
add_executable(h5power    h5power.cc)
add_executable(makecoefs  makecoefs.cc)

foreach(program ${bin_PROGRAMS})
  target_link_libraries(${program} expcoef exputil ${common_LINKLIB})
  target_include_directories(${program} PUBLIC ${common_INCLUDE})
  install(TARGETS ${program} DESTINATION bin)
endforeach()

pybind11_add_module(pyEXP PyWrappers.cc CoefWrappers.cc
  BasisWrappers.cc FieldWrappers.cc ParticleReaderWrappers.cc)
target_link_libraries(pyEXP PUBLIC expcoef exputil ${common_LINKLIB})
target_include_directories(pyEXP PUBLIC ${common_INCLUDE})

install(TARGETS pyEXP
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages")
