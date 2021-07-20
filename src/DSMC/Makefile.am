## Process this file with automake to produce Makefile.in

CC = mpicc
CXX = mpiCC

LIBS = -lpthread -lmpi -lm
EXTRA_LDFLAGS = -L$(top_srcdir)/exputil/.libs

if ENABLE_DEBUG
ADD_CFLAGS2   = -DDEBUG
ADD_CXXFLAGS2 = -DDEBUG
endif

if ENABLE_DMALLOC
ADD_CFLAGS3   = -DUSE_DMALLOC -g
ADD_CXXFLAGS3 = -DUSE_DMALLOC -g
endif

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src		\
-I$(top_srcdir)/src/DSMC -I$(top_srcdir)/include $(CUDA_CFLAGS)
AM_CFLAGS   = -fPIC -D_REENTRANT $(INCLUDES) $(ADD_CFLAGS2) $(ADD_CFLAGS3)
AM_CXXFLAGS = -fPIC -D_REENTRANT $(INCLUDES) $(ADD_CXXFLAGS2) $(ADD_CXXFLAGS3) @exp_vtk_include@ @exp_eigen3_include@ $(CUDA_FLAGS)

bin_PROGRAMS = testCrossSection testIonRecomb genIonRecomb	\
testQuantile testExciteRate testRatio makeRecombCache

if ENABLE_CUDA
  ADD_CUDA0  = cudaCollideIon.cu
  ADD_CUDA1  = cudaIon.cu
  ifeq ($(GENCODE_FLAGS),)
    $(foreach sm,@sm_vers@,$(eval GENCODE_FLAGS += -gencode arch=compute_$(sm),code=sm_$(sm)))
  endif
  CUDA_LIBFLAGS  = -lcudart -lnvToolsExt
  NVCC_CXXFLAGS  = -I/usr/include/mpi -I. $(GENCODE_FLAGS) $(NVCC_EXTRA) @exp_eigen3_include@ -Xcompiler $(AM_CPPFLAGS)
  bin_PROGRAMS += testCrossCuda
endif

noinst_LIBRARIES = libexpdsmc.a libtestdsmc.a

DSMC_SRC = TreeDSMC.cc Collide.cc CollideLTE.cc HeatCool.cc		\
	CollideIon.cc Ion.cc interactSelect.cc Elastic.cc TopBase.cc	\
	NTC.cc Quantile.cc QuantileBag.cc phfit2.f Badnell.cc		\
	atomic_constants.cc spline.cc $(ADD_CUDA0)

TESTDSMC_SRC = Ion.cc interactSelect.cc Elastic.cc TopBase.cc	\
	Quantile.cc QuantileBag.cc phfit2.f atomic_constants.cc	\
	Badnell.cc spline.cc $(ADD_CUDA1)

EXPFLAGS = -L$(top_srcdir)/exputil/.libs -L$(top_srcdir)/src/DSMC $(BOOST_LDFLAGS)

libexpdsmc_a_SOURCES = $(DSMC_SRC) 
libexpdsmc_a_LDFLAGS = $(EXPFLAGS)
# libexpdsmc_a_LIBADD =  $(ADD_CUDA_LINK0) # $(CUDA_LDFLAGS) $(CUDA_LIBFLAGS) $(LIBS)

libtestdsmc_a_SOURCES = $(TESTDSMC_SRC)
libtestdsmc_a_LDFLAGS = $(EXPFLAGS)
# libtestdsmc_a_LIBADD =   $(ADD_CUDA_LINK1) # $(CUDA_LDFLAGS) $(CUDA_LIBFLAGS) $(LIBS)

testCrossSection_SOURCES = testCross.cc
testCrossSection_DEPENDENCIES = libtestdsmc.a
testCrossSection_LDFLAGS = $(EXTRA_LDFLAGS) $(BOOST_LDFLAGS) $(BUILD_FLAGS) @exp_vtk_lib@ $(CUDA_LDFLAGS)
testCrossSection_LDADD   = $(BOOST_FILESYSTEM_LIB) $(BOOST_SERIALIZATION_LIB) $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_SYSTEM_LIB) libtestdsmc.a -ltestdsmc -lexputil @exp_vtk_ldlib@ $(LAPACK_LIBS) $(BLAS_LIBS) $(CUDA_LIBFLAGS)

testIonRecomb_SOURCES = testIonRecomb.cc
testIonRecomb_DEPENDENCIES = libtestdsmc.a
testIonRecomb_LDFLAGS = $(EXPFLAGS) $(BUILD_FLAGS) @exp_vtk_lib@ $(CUDA_LDFLAGS)
testIonRecomb_LDADD   = $(BOOST_FILESYSTEM_LIB) $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_SYSTEM_LIB) libtestdsmc.a -lexputil @exp_vtk_ldlib@ $(LAPACK_LIBS) $(BLAS_LIBS) $(CUDA_LIBFLAGS)

genIonRecomb_SOURCES = genIonRecomb.cc
genIonRecomb_DEPENDENCIES = libtestdsmc.a
genIonRecomb_LDFLAGS = $(EXPFLAGS) $(BUILD_FLAGS) @exp_vtk_lib@ $(CUDA_LDFLAGS)
genIonRecomb_LDADD   = $(BOOST_FILESYSTEM_LIB) $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_SYSTEM_LIB) libtestdsmc.a -lexputil @exp_vtk_ldlib@ $(LAPACK_LIBS) $(BLAS_LIBS) $(CUDA_LIBFLAGS)

testQuantile_SOURCES = testQuantile.cc
testQuantile_DEPENDENCIES = libtestdsmc.a
testQuantile_LDFLAGS = $(EXPFLAGS) $(BUILD_FLAGS) $(CUDA_LDFLAGS)
testQuantile_LDADD   = libtestdsmc.a $(BOOST_PROGRAM_OPTIONS_LIB) $(CUDA_LIBFLAGS) $(BOOST_SYSTEM_LIB) 

testCrossCuda_SOURCES = testCrossCuda.cc
testCrossCuda_DEPENDENCIES = libtestdsmc.a
testCrossCuda_LDFLAGS = $(EXTRA_LDFLAGS) $(BOOST_LDFLAGS) $(BUILD_FLAGS) @exp_vtk_lib@ $(CUDA_LDFLAGS)
testCrossCuda_LDADD   = $(BOOST_FILESYSTEM_LIB) $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_SYSTEM_LIB) libtestdsmc.a -lexputil @exp_vtk_ldlib@ $(LAPACK_LIBS) $(BLAS_LIBS) $(CUDA_LIBFLAGS)

testExciteRate_SOURCES = testExciteRate.cc
testExciteRate_DEPENDENCIES = libtestdsmc.a
testExciteRate_LDFLAGS = $(EXTRA_LDFLAGS) $(BOOST_LDFLAGS) $(BUILD_FLAGS) @exp_vtk_lib@ $(CUDA_LDFLAGS)
testExciteRate_LDADD   = $(BOOST_FILESYSTEM_LIB) $(BOOST_SERIALIZATION_LIB) $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_SYSTEM_LIB) libtestdsmc.a -ltestdsmc -lexputil @exp_vtk_ldlib@ $(LAPACK_LIBS) $(BLAS_LIBS) $(CUDA_LIBFLAGS)

testRatio_SOURCES = testRatio.cc
testRatio_DEPENDENCIES = libtestdsmc.a
testRatio_LDFLAGS = $(EXPFLAGS) $(BUILD_FLAGS) @exp_vtk_lib@ $(CUDA_LDFLAGS)
testRatio_LDADD   = $(BOOST_FILESYSTEM_LIB) $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_SYSTEM_LIB) libtestdsmc.a -lexputil @exp_vtk_ldlib@ $(LAPACK_LIBS) $(BLAS_LIBS) $(CUDA_LIBFLAGS)

makeRecombCache_SOURCES = makeRecombCache.cc
makeRecombCache_DEPENDENCIES = libtestdsmc.a
makeRecombCache_LDFLAGS = $(EXPFLAGS) $(BUILD_FLAGS) @exp_vtk_lib@ $(CUDA_LDFLAGS)
makeRecombCache_LDADD   = $(BOOST_FILESYSTEM_LIB) $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_SYSTEM_LIB) libtestdsmc.a -lexputil @exp_vtk_ldlib@ $(LAPACK_LIBS) $(BLAS_LIBS) $(CUDA_LIBFLAGS) $(YAML)

# the rule that can recognize *.cu files and compile the *.cu files to
# *.lo and *.o files

.cu.lo:
	$(top_builddir)/cudalt.py $@ $(NVCC) -c $(NVCC_CXXFLAGS) $<

.cu.o :
	$(NVCC) $(CUDA_CFLAGS) $(NVCC_CXXFLAGS) -c -o $@ $<

CLEANFILES =

## libexpdsmc/Makefile.am ends here
