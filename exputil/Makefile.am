# Process this file with automake to produce Makefile.in

CC = mpicc
CXX = mpiCC

LIBS = -lpthread -lmpi -lm

if ENABLE_DEBUG
ADD_CFLAGS2   = -DDEBUG
ADD_CXXFLAGS2 = -DDEBUG
endif

if ENABLE_DMALLOC
ADD_CFLAGS3   = -DUSE_DMALLOC -g
ADD_CXXFLAGS3 = -DUSE_DMALLOC -g
endif

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/include $(CUDA_CFLAGS)
AM_CFLAGS   = -fPIC -D_REENTRANT $(INCLUDES) $(ADD_CFLAGS2) $(ADD_CFLAGS3)
AM_CXXFLAGS = -fPIC -D_REENTRANT $(INCLUDES) $(ADD_CXXFLAGS2) $(ADD_CXXFLAGS3) @exp_vtk_include@ $(CUDA_FLAGS) @exp_eigen3_include@

if ENABLE_CUDA
  CUDA_SRC       = cudaSLGridMP2.cu
  ifeq ($(GENCODE_FLAGS),)
    $(foreach sm,@sm_vers@,$(eval GENCODE_FLAGS += -gencode arch=compute_$(sm),code=sm_$(sm)))
  endif
  CUDA_LIBFLAGS  = -lcudart
  NVCC_CXXFLAGS  = -I/usr/include/mpi -I. $(GENCODE_FLAGS) $(NVCC_EXTRA) -Xcompiler $(AM_CPPFLAGS)
endif

lib_LTLIBRARIES = libexputil.la

ODE_SRC = rk4.cc bs.cc odesolve.cc
SYMP_SRC = sia4.cc
ROOT_SRC = rtbis.cc rtsafe.cc brent.cc zbrent.cc mnbrak.cc zbrak.cc zbrac.cc
QUAD_SRC = qadapt.cc gauleg.cc qadapt2d.cc gint2.cc rombe2d.cc Trapsum.cc
UTIL_SRC = nrutil.cc elemfunc.cc euler.cc euler_slater.cc		\
	rotmatrix.cc wordSplit.cc FileUtils.cc BarrierWrapper.cc	\
	stack.cc VtkGrid.cc localmpi.cc
OPTIMIZATION_SRC = simann2.cc
POLY_SRC = zroots.cc laguer.cc poly.cc cpoly.cc horner.cc syndiv.cc
PHASE_SRC = phase.cc ensemble.cc io_ensemble.cc move_ensemble.cc \
	diag_ensemble.cc map.cc pc.cc models.cc prolate.cc perfect.cc \
	rotcurv.cc needle.cc hubble.cc quad.cc
SPECFUNC_SRC = gammln.cc bessel.cc OrthoPoly.cc CauchyPV.cc # modbessel.cc 
INTERP_SRC = Spline.cc SplintE.cc Vodd2.cc Vlocate.cc levsurf.cc \
	ChebFit.cc Interp1d.cc Cheby1d.cc
MASSMODEL_SRC = massmodel.cc massmodel_dist.cc embedded.cc isothermal.cc \
        realize_model.cc toomre.cc exponential.cc
ORBIT_SRC = orbit.cc orbit_trans.cc FindOrb.cc
BIORTH_SRC = biorth_wake.cc biorth.cc biorth2d.cc biorth_grid.cc \
	sbessz.cc ultra.cc jnsph.cc bessz.cc sphereSL.cc biorth1d.cc \
	biorth_wake_orientation.cc SLGridMP2.cc scalarprod.cc
GAUSS_SRC = gaussQ.cc GaussCore.c Hermite.c Jacobi.c Laguerre.c
QPDISTF_SRC = QPDistF.cc qld.c
SLEDGE_SRC = sledge.f

#	all source files

libexputil_la_SOURCES = $(ODE_SRC) $(ROOT_SRC) $(QUAD_SRC) \
	$(RANDOM_SRC) $(UTIL_SRC) $(SPECFUNC_SRC) \
	$(PHASE_SRC)  $(SYMP_SRC) $(INTERP_SRC) $(MASSMODEL_SRC) \
	$(ORBIT_SRC) $(BIORTH_SRC) $(POLY_SRC) $(GAUSS_SRC) \
	$(QPDISTF_SRC) $(BESSEL_SRC) $(OPTIMIZATION_SRC) \
	$(SLEDGE_SRC) $(CUDA_SRC)

libexputil_la_LIBADD = $(LAPACK_LIBS) $(BLAS_LIBS) $(CUDA_LDFLAGS)	\
$(CUDA_LIBFLAGS) @exp_vtk_ldlib@ @exp_vtk_lib@ $(LIBS) $(FLIBS)

# the rule that can recognize *.cu files and compile the *.cu files to
# *.lo and *.o files

.cu.lo:
	$(top_builddir)/cudalt.py $@ $(NVCC) -c $(NVCC_CXXFLAGS) $<

.cu.o :
	$(NVCC) $(CUDA_CFLAGS) $(NVCC_CXXFLAGS) -c -o $@ $<


# special case for f77 for now

sledge.o:	sledge.f
		f77 -O2 -g -o sledge.o -c sledge.f

## libexputil/Makefile.am ends here

