# Make the expui shared library
set(expui_SOURCES BasisFactory.cc BiorthBasis.cc BiorthBess.cc FieldBasis.cc
  CoefContainer.cc CoefStruct.cc FieldGenerator.cc expMSSA.cc
  Coefficients.cc KMeans.cc Centering.cc ParticleIterator.cc
  Koopman.cc)

add_library(expui ${expui_SOURCES})
add_dependencies(expui exputil)

set(common_LINKLIBS ${DEP_LIB} ${OpenMP_CXX_LIBRARIES} 
  ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}
  ${MPI_CXX_LIBRARIES}
  ${YAML_CPP_LIBRARIES} ${EXPUTIL}
  ${VTK_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES}
  HighFive)

set(common_INCLUDE 
  ${DEP_INC} ${HDF5_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR} ${FFTW_INCLUDE_DIRS}
  $<INSTALL_INTERFACE:include/exputil>)

set(common_INCLUDE_PRIVATE
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/config/>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

#   if(HighFive_FOUND)
#   list(APPEND common_LINKLIBS HighFive)
# else()
#   list(APPEND common_INCLUDE_PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/extern/HighFive/include>)
# endif()

if(USE_SUBMODULES)
  list(APPEND common_INCLUDE_PRIVATE 
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/extern/yaml-cpp/include>)
  list(APPEND common_LINKLIBS yaml-cpp)
endif()

find_package(EXPutil QUIET)

if(EXPutil_FOUND)
  list(PREPEND common_LINKLIBS EXP::EXPutil)
  list(PREPEND common_INCLUDE ${EXPutil_INCLUDE_DIRS})
else()
  list(PREPEND common_LINKLIBS exputil) #internal name
  list(PREPEND common_INCLUDE_PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>)
endif()

if(PNG_FOUND AND ENABLE_PNG)
  list(APPEND common_LINKLIBS PNG::PNG)
  list(APPEND common_INCLUDE ${PNG_INCLUDE_DIRS})
  add_definitions(${PNG_DEFINITIONS})
endif()

if(ENABLE_CUDA)
  list(APPEND common_INCLUDE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} ${CUDAToolkit_INCLUDE_DIRS})
  list(APPEND common_LINKLIBS CUDA::cudart CUDA::nvToolsExt)
endif()

if(SLURM_FOUND)
  list(APPEND common_LINKLIBS ${SLURM_LIBRARY})
endif()

# TODO: Replace this globbing with a list of actually public headers
file(GLOB expui_PUBLIC_HEADERS "*.[hH]")
set_target_properties (expui PROPERTIES 
  OUTPUT_NAME expui
  EXPORT_NAME EXPui
  PUBLIC_HEADER "${expui_PUBLIC_HEADERS}")

target_include_directories(expui PUBLIC ${common_INCLUDE})
target_include_directories(expui PRIVATE ${common_INCLUDE_PRIVATE})
target_link_libraries(expui PUBLIC ${common_LINKLIBS})

# Configure and build the test routines
if (ENABLE_TESTS)
  set(bin_PROGRAMS nativetoh5 h5compare viewcoefs h5power makecoefs testread)

  add_executable(nativetoh5 coefstoh5.cc)
  add_executable(h5compare  h5compare.cc)
  add_executable(viewcoefs  viewcoefs.cc)
  add_executable(h5power    h5power.cc)
  add_executable(makecoefs  makecoefs.cc)
  add_executable(testread   testread.cc)

  foreach(program ${bin_PROGRAMS})
    target_link_libraries(${program} PUBLIC expui ${common_LINKLIBS})
    target_include_directories(${program} PUBLIC ${common_INCLUDE})
    target_include_directories(${program} PRIVATE ${common_INCLUDE_PRIVATE})
    target_compile_options(${program} PUBLIC ${OpenMP_CXX_FLAGS})
  endforeach()
endif()

install (TARGETS expui
  EXPORT EXPuiTargets
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  PUBLIC_HEADER DESTINATION include/expui)

if (INSTALL_CMAKE_FIND)
  # Set up cmake configs
  set(CMAKEPACKAGE_INSTALL_DIR share/EXPui/cmake)
  install (EXPORT EXPuiTargets
    FILE EXPuiTargets.cmake
    NAMESPACE EXP::
    DESTINATION ${CMAKEPACKAGE_INSTALL_DIR})

  include(CMakePackageConfigHelpers)
  configure_package_config_file(EXPuiConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/EXPuiConfig.cmake
    INSTALL_DESTINATION ${CMAKEPACKAGE_INSTALL_DIR}
    PATH_VARS CMAKEPACKAGE_INSTALL_DIR CMAKE_INSTALL_PREFIX
  )

  write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/EXPuiVersion.cmake
                                  VERSION ${CMAKE_PROJECT_VERSION}
                                  COMPATIBILITY SameMinorVersion )

  install (FILES ${CMAKE_CURRENT_BINARY_DIR}/EXPuiConfig.cmake
                ${CMAKE_CURRENT_BINARY_DIR}/EXPuiVersion.cmake
          DESTINATION ${CMAKEPACKAGE_INSTALL_DIR} )
endif()
