
set(ODE_SRC rk4.cc bs.cc odesolve.cc)
set(SYMP_SRC sia4.cc)
set(ROOT_SRC rtbis.cc rtsafe.cc brent.cc zbrent.cc mnbrak.cc zbrak.cc zbrac.cc)
set(QUAD_SRC qadapt.cc gauleg.cc qadapt2d.cc gint2.cc rombe2d.cc Trapsum.cc)
set(UTIL_SRC nrutil.cc elemfunc.cc euler.cc euler_slater.cc
  rotmatrix.cc wordSplit.cc FileUtils.cc BarrierWrapper.cc
  stack.cc VtkGrid.cc localmpi.cc)
set(OPTIMIZATION_SRC simann2.cc)
set(POLY_SRC zroots.cc laguer.cc poly.cc cpoly.cc horner.cc syndiv.cc)
set(PHASE_SRC phase.cc ensemble.cc io_ensemble.cc move_ensemble.cc
  diag_ensemble.cc map.cc pc.cc models.cc prolate.cc perfect.cc
  rotcurv.cc needle.cc hubble.cc quad.cc)
set(SPECFUNC_SRC gammln.cc bessel.cc OrthoPoly.cc CauchyPV.cc) # modbessel.cc
set(INTERP_SRC Spline.cc SplintE.cc Vodd2.cc Vlocate.cc levsurf.cc ChebFit.cc Interp1d.cc Cheby1d.cc)
set(MASSMODEL_SRC massmodel.cc massmodel_dist.cc embedded.cc isothermal.cc realize_model.cc
  toomre.cc exponential.cc)
set(ORBIT_SRC orbit.cc orbit_trans.cc FindOrb.cc)
set(BIORTH_SRC biorth_wake.cc biorth.cc biorth2d.cc biorth_grid.cc
  sbessz.cc ultra.cc jnsph.cc bessz.cc sphereSL.cc biorth1d.cc
  biorth_wake_orientation.cc SLGridMP2.cc scalarprod.cc)
set(GAUSS_SRC gaussQ.cc GaussCore.c Hermite.c Jacobi.c Laguerre.c)
set(QPDISTF_SRC QPDistF.cc qld.c)
set(SLEDGE_SRC sledge.f)
set(CUDA_SRC cudaSLGridMP2.cu)

set(exputil_SOURCES ${ODE_SRC} ${ROOT_SRC} ${QUAD_SRC}
  ${RANDOM_SRC} ${UTIL_SRC} ${SPECFUNC_SRC}
  ${PHASE_SRC}  ${SYMP_SRC} ${INTERP_SRC} ${MASSMODEL_SRC}
  ${ORBIT_SRC} ${BIORTH_SRC} ${POLY_SRC} ${GAUSS_SRC}
  ${QPDISTF_SRC} ${BESSEL_SRC} ${OPTIMIZATION_SRC}
  ${SLEDGE_SRC} ${CUDA_SRC})

# shared lib
add_library(exputil SHARED ${exputil_SOURCES})
set_target_properties(exputil PROPERTIES OUTPUT_NAME exputil)
target_include_directories(
  exputil
  PUBLIC $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/>
  ${CMAKE_BINARY_DIR} ${DEP_INC}
  ${VTK_INCLUDE})
target_link_libraries(exputil PUBLIC ${DEP_LIB} OpenMP::OpenMP_CXX MPI::MPI_CXX Boost::random ${CUDA_LIBRARIES})

add_library(exputil_static STATIC ${exputil_SOURCES})
set_target_properties(exputil_static PROPERTIES OUTPUT_NAME exputil)
target_include_directories(
  exputil_static
  PUBLIC $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/>
  ${CMAKE_BINARY_DIR} ${DEP_INC}
  ${VTK_INCLUDE})
target_link_libraries(exputil_static PUBLIC ${DEP_LIB} OpenMP::OpenMP_CXX MPI::MPI_CXX Boost::random ${CUDA_LIBRARIES})
