set(ODE_SRC rk4.cc bs.cc odesolve.cc)
set(SYMP_SRC sia4.cc)
set(ROOT_SRC rtbis.cc rtsafe.cc brent.cc zbrent.cc mnbrak.cc zbrak.cc zbrac.cc)
set(QUAD_SRC qadapt.cc gauleg.cc qadapt2d.cc gint2.cc rombe2d.cc
  Trapsum.cc Splsum.cc)
set(UTIL_SRC nrutil.cc elemfunc.cc euler.cc euler_slater.cc # Hankel.cc
  rotmatrix.cc wordSplit.cc FileUtils.cc BarrierWrapper.cc stack.cc
  localmpi.cc TableGrid.cc writePVD.cc libvars.cc TransformFFT.cc QDHT.cc
  YamlCheck.cc parseVersionString.cc EXPmath.cc laguerre_polynomial.cpp
  YamlConfig.cc orthoTest.cc OrthoFunction.cc)

if(HAVE_VTK)
  list(APPEND UTIL_SRC VtkGrid.cc VtkPCA.cc)
endif()

set(OPTIMIZATION_SRC SimAnn.cc)
set(POLY_SRC zroots.cc laguer.cc poly.cc cpoly.cc horner.cc syndiv.cc)
set(PHASE_SRC phase.cc ensemble.cc io_ensemble.cc move_ensemble.cc
  diag_ensemble.cc map.cc pc.cc models.cc prolate.cc perfect.cc
  rotcurv.cc needle.cc hubble.cc quad.cc)
set(SPECFUNC_SRC gammln.cc bessel.cc OrthoPoly.cc CauchyPV.cc) # modbessel.cc
set(INTERP_SRC Spline.cc SplintE.cc Vodd2.cc Vlocate.cc levsurf.cc Interp1d.cc Cheby1d.cc MonotCubicInterpolator.cc)
set(MASSMODEL_SRC massmodel.cc massmodel_dist.cc embedded.cc isothermal.cc realize_model.cc GenPoly.cc mestel.cc
  toomre.cc exponential.cc)
set(ORBIT_SRC orbit.cc orbit_trans.cc FindOrb.cc)

set(BIORTH_SRC biorth_wake.cc biorth.cc biorth2d.cc biorth_grid.cc
  sbessz.cc ultra.cc bessz.cc sphereSL.cc biorth1d.cc Coefs.cc
  biorth_wake_orientation.cc SLGridMP2.cc scalarprod.cc EmpCylSL.cc
  EmpCyl2d.cc BiorthCyl.cc BiorthCube.cc)

set(GAUSS_SRC gaussQ.cc GaussCore.c Hermite.c Jacobi.c Laguerre.c)
set(QPDISTF_SRC QPDistF.cc qld.c)
set(SLEDGE_SRC sledge.f)
set(PARTICLE_SRC Particle.cc ParticleReader.cc header.cc)
set(CUDA_SRC cudaParticle.cu cudaSLGridMP2.cu)

set(exputil_SOURCES ${ODE_SRC} ${ROOT_SRC} ${QUAD_SRC}
  ${RANDOM_SRC} ${UTIL_SRC} ${SPECFUNC_SRC}
  ${PHASE_SRC}  ${SYMP_SRC} ${INTERP_SRC} ${MASSMODEL_SRC}
  ${ORBIT_SRC} ${BIORTH_SRC} ${POLY_SRC} ${GAUSS_SRC}
  ${QPDISTF_SRC} ${BESSEL_SRC} ${OPTIMIZATION_SRC}
  ${SLEDGE_SRC} ${PARTICLE_SRC} ${CUDA_SRC})

set(common_INCLUDE 
  $<INSTALL_INTERFACE:include/exputil>
  ${YAML_CPP_INCLUDE_DIR} ${DEP_INC} ${EIGEN3_INCLUDE_DIR}
  ${HDF5_INCLUDE_DIRS} ${FFTW_INCLUDE_DIRS})

set(common_INCLUDE_PRIVATE
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/config/>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/extern/>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/>)

set(common_LINKLIBS
  ${DEP_LIB} ${OpenMP_CXX_LIBRARIES} ${MPI_CXX_LIBRARIES}
  ${YAML_CPP_LIBRARIES} ${VTK_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES}
  ${FFTW_DOUBLE_LIB} HighFive)

# if(HighFive_FOUND)
#   list(APPEND common_LINKLIBS HighFive)
# else()
#   list(APPEND common_INCLUDE_PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/extern/HighFive/include>)
# endif()

if(USE_SUBMODULES)
  list(APPEND common_INCLUDE_PRIVATE 
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/extern/yaml-cpp/include>)
  list(APPEND common_LINKLIBS yaml-cpp)
endif()

if(ENABLE_CUDA)
  list(APPEND common_INCLUDE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} ${CUDAToolkit_INCLUDE_DIRS})
  list(APPEND common_LINKLIBS CUDA::cudart CUDA::nvToolsExt)
endif()

if(ENABLE_XDR AND TIRPC_FOUND)
  list(APPEND common_INCLUDE ${TIRPC_INCLUDE_DIRS})
  list(APPEND common_LINKLIBS ${TIRPC_LIBRARIES})
endif()

if(PNG_FOUND AND ENABLE_PNG)
  list(APPEND common_LINKLIBS PNG::PNG)
  list(APPEND common_INCLUDE ${PNG_INCLUDE_DIRS})
  add_definitions(${PNG_DEFINITIONS})
endif()

# TODO: Replace this globbing with a list of actually public headers
file(GLOB exputil_PUBLIC_HEADERS "../include/*.[hH]")

# shared lib
add_library(exputil ${exputil_SOURCES} ${exputil_PUBLIC_HEADERS})

set_target_properties(exputil PROPERTIES 
  OUTPUT_NAME exputil
  EXPORT_NAME EXPutil
  PUBLIC_HEADER "${exputil_PUBLIC_HEADERS}")

target_include_directories(exputil PUBLIC ${common_INCLUDE})
target_include_directories(exputil PRIVATE ${common_INCLUDE_PRIVATE})
target_link_libraries(exputil PUBLIC ${common_LINKLIBS})

install(TARGETS exputil
  EXPORT EXPutilTargets
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include/exputil)

install (DIRECTORY $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/extern/png++>
  DESTINATION include/exputil)

# install (DIRECTORY $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/extern/HighFive/include/highfive>
#   DESTINATION include/exputil)

install (DIRECTORY $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/extern/rapidxml>
  DESTINATION include/exputil)

if (INSTALL_CMAKE_FIND)
  # Set up cmake configs
  set(CMAKEPACKAGE_INSTALL_DIR share/EXPutil/cmake)
  install (EXPORT EXPutilTargets
    FILE EXPutilTargets.cmake
    NAMESPACE EXP::
    DESTINATION ${CMAKEPACKAGE_INSTALL_DIR})

  include(CMakePackageConfigHelpers)
  configure_package_config_file(EXPutilConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/EXPutilConfig.cmake
    INSTALL_DESTINATION ${CMAKEPACKAGE_INSTALL_DIR}
    PATH_VARS CMAKEPACKAGE_INSTALL_DIR CMAKE_INSTALL_PREFIX
    )
  write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/EXPutilVersion.cmake
                                  VERSION ${CMAKE_PROJECT_VERSION}
                                  COMPATIBILITY SameMinorVersion )
  install (FILES ${CMAKE_CURRENT_BINARY_DIR}/EXPutilConfig.cmake
                ${CMAKE_CURRENT_BINARY_DIR}/EXPutilVersion.cmake
          DESTINATION ${CMAKEPACKAGE_INSTALL_DIR})
endif()
