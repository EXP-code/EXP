find_package(Python3 REQUIRED)
find_package(PythonExtensions QUIET)

# Copy Python init file to build directory
configure_file(__init__.py __init__.py COPYONLY)

# Default build environment
set(common_LINKLIBS ${OpenMP_CXX_LIBRARIES} ${MPI_CXX_LIBRARIES}
  ${VTK_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES}
  $<$<BOOL:${HighFive_FOUND}>:HighFive>)

set(common_INCLUDE $<INSTALL_INTERFACE:include>
  ${DEP_INC} ${CMAKE_CURRENT_SOURCE_DIR}
  ${HDF5_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR}
  ${EXPutil_INCLUDE_DIRS} ${EXPui_INCLUDE_DIRS} ${EXPcommon_INCLUDE_DIRS})

set(common_INCLUDE_PRIVATE
    ${CMAKE_BINARY_DIR}/config)

find_package(EXPutil QUIET)
if(EXPutil_FOUND)
  list(APPEND common_LINKLIBS EXP::EXPutil)
  list(APPEND common_INCLUDE ${EXPutil_INCLUDE_DIRS})
else()
  list(APPEND common_LINKLIBS exputil) #internal name
  list(APPEND common_INCLUDE_PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>)
endif()

find_package(EXPui QUIET)
if(EXPui_FOUND)
  list(APPEND common_LINKLIBS EXP::EXPui)
  list(APPEND common_INCLUDE ${EXPui_INCLUDE_DIRS})
else()
  list(APPEND common_LINKLIBS expui) #internal name
  list(APPEND common_INCLUDE_PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/expui/>)
endif()

if(USE_SUBMODULES)
  list(APPEND common_INCLUDE_PRIVATE 
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/extern/yaml-cpp/include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/extern/HighFive/include>)
  list(APPEND common_LINKLIBS yaml-cpp)
endif()

if(ENABLE_CUDA)
  list(APPEND common_LINKLIBS CUDA::cudart CUDA::nvToolsExt)
endif()

if(SLURM_FOUND)
  list(APPEND common_LINKLIBS ${SLURM_LIBRARY})
endif()

# Build and configure the installation of the Python bindings
pybind11_add_module(_pyEXP PyWrappers.cc CoefWrappers.cc
  UtilWrappers.cc BasisWrappers.cc FieldWrappers.cc
  ParticleReaderWrappers.cc MSSAWrappers.cc EDMDWrappers.cc)
target_link_libraries(_pyEXP PUBLIC ${common_LINKLIBS} )
target_include_directories(_pyEXP PUBLIC ${common_INCLUDE})
target_include_directories(_pyEXP PRIVATE ${common_INCLUDE_PRIVATE})
target_compile_options(_pyEXP PUBLIC ${OpenMP_CXX_FLAGS})

if(PYTHON_RELATIVE_SITE_PACKAGES_DIR)
  install(TARGETS _pyEXP DESTINATION ${PYTHON_RELATIVE_SITE_PACKAGES_DIR}/pyEXP)
  install(FILES "__init__.py" DESTINATION ${PYTHON_RELATIVE_SITE_PACKAGES_DIR}/pyEXP)
else()
  install(TARGETS _pyEXP
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION "lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages/pyEXP"
    ARCHIVE DESTINATION "lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages/pyEXP")
  install(FILES "__init__.py" DESTINATION  lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages/pyEXP)
endif()

