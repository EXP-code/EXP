find_package(Python3 COMPONENTS Development REQUIRED)
find_package(pybind11 REQUIRED)
find_package(EXPutil)
find_package(EXPcoefs)
find_package(EXPcommon)
message(STATUS "EXPutil_FOUND=${EXPutil_FOUND} ${EXPutil_INCLUDE_DIRS} ${EXPutil_INCLUDE_DIR}")

# Default build environment
#
# set(common_LINKLIB OpenMP::OpenMP_CXX MPI::MPI_CXX yaml-cpp
#   EXP::EXPutil EXP::EXPcoefs EXP::EXPcommon HighFive ${VTK_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
set(common_LINKLIB ${OpenMP_CXX_LIBRARIES} ${MPI_CXX_LIBRARIES}
  # ${YAML_CPP_LIBRARIES} # Not directly used in pyexp
  ${VTK_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES}
  $<IF:$<BOOL:${EXPutil_LIBRARIES}>,${EXPutil_LIBRARIES},exputil>
  $<IF:$<BOOL:${EXPcoefs_LIBRARIES}>,${EXPcoefs_LIBRARIES},expcoef>
  $<$<BOOL:${HighFive_FOUND}>:HighFive>)
  # $<$<BOOL:${PNG_FOUND}>:${PNG_LIBRARIES}>)


# if(PNG_FOUND)
#   list(APPEND common_LINKLIB PNG::PNG)
# endif()

set(common_INCLUDE $<INSTALL_INTERFACE:include>
  ${DEP_INC} ${CMAKE_CURRENT_SOURCE_DIR}
  ${HDF5_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR}
  ${EXPutil_INCLUDE_DIRS} ${EXPcoefs_INCLUDE_DIRS} ${EXPcommon_INCLUDE_DIRS})

set(common_INCLUDE_DIRS_PRIVATE
    ${CMAKE_BINARY_DIR}/config)

if (NOT (EXPutil_FOUND AND EXPcoefs_FOUND))
  list(APPEND common_INCLUDE_DIRS_PRIVATE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
endif()

if(ENABLE_CUDA)
  list(APPEND common_LINKLIB CUDA::cudart CUDA::nvToolsExt)
endif()

if(SLURM_FOUND)
  list(APPEND common_LINKLIB ${SLURM_LIBRARY})
endif()

# Not used in pyEXP
# if(ENABLE_XDR AND TIRPC_FOUND)
#   list(APPEND common_LINKLIB ${TIRPC_LIBRARIES})
# endif()

# Build and configure the installation of the Python bindings
#
pybind11_add_module(_pyEXP PyWrappers.cc CoefWrappers.cc
  UtilWrappers.cc BasisWrappers.cc FieldWrappers.cc
  ParticleReaderWrappers.cc MSSAWrappers.cc EDMDWrappers.cc)
target_link_libraries(_pyEXP PUBLIC ${common_LINKLIB})
target_include_directories(_pyEXP PUBLIC ${common_INCLUDE})
target_include_directories(_pyEXP PRIVATE ${common_INCLUDE_DIRS_PRIVATE})
target_compile_options(_pyEXP PUBLIC ${OpenMP_CXX_FLAGS})
get_target_property(cxxflags _pyEXP COMPILE_OPTIONS)
message("The project has set the following flags: ${cxxflags}")

install(TARGETS _pyEXP LIBRARY DESTINATION pyEXP)
# install(TARGETS _pyEXP
#   # RUNTIME DESTINATION bin
#   # LIBRARY DESTINATION "lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages"
#   # ARCHIVE DESTINATION "lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages"
# )

