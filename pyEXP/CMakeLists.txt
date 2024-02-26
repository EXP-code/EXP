find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 REQUIRED)

# Default build environment
#
set(common_LINKLIB OpenMP::OpenMP_CXX MPI::MPI_CXX yaml-cpp
  exputil ${VTK_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})

if(PNG_FOUND)
  list(APPEND common_LINKLIB PNG::PNG)
endif()

set(common_INCLUDE $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/>
  ${CMAKE_BINARY_DIR} ${DEP_INC} ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/.. ${HighFive_SOURCE_DIR}/include
  ${HDF5_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})

if(ENABLE_CUDA)
  list(APPEND common_LINKLIB CUDA::cudart CUDA::nvToolsExt)
endif()

if(SLURM_FOUND)
  list(APPEND common_LINKLIB ${SLURM_LIBRARY})
endif()

if(ENABLE_XDR AND TIRPC_FOUND)
  list(APPEND common_LINKLIB ${TIRPC_LIBRARIES})
endif()

# Build and configure the installation of the Python bindings
#
pybind11_add_module(pyEXP PyWrappers.cc CoefWrappers.cc
  UtilWrappers.cc BasisWrappers.cc FieldWrappers.cc
  ParticleReaderWrappers.cc MSSAWrappers.cc EDMDWrappers.cc)
target_link_libraries(pyEXP PUBLIC expcoef exputil ${common_LINKLIB})
target_include_directories(pyEXP PUBLIC ${common_INCLUDE})
target_compile_options(pyEXP PUBLIC ${OpenMP_CXX_FLAGS})
get_target_property(cxxflags pyEXP COMPILE_OPTIONS)
message("The project has set the following flags: ${cxxflags}")

# install(TARGETS pyEXP
#   # RUNTIME DESTINATION bin
#   # LIBRARY DESTINATION "lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages"
#   # ARCHIVE DESTINATION "lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages"
# )

