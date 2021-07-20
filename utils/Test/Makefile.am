## Process this file with automake to produce Makefile.in

# VPATH = $(top_srcdir)/src:.

CC = mpicc
CXX = mpiCC

if ENABLE_DMALLOC
ADD_CFLAGS3 = -DUSE_DMALLOC -g
ADD_CXXFLAGS3 = -DUSE_DMALLOC -g
ADD_LIBS3 = -ldmallocth -ldmallocthcxx
endif

# Need for lib exputil which contains cuda code
#
if ENABLE_CUDA
  CUDA_LIBFLAGS  = -lcudart -lnvToolsExt
endif

AM_CFLAGS = -D_REENTRANT $(ADD_CFLAGS1) $(ADD_CFLAGS2) $(ADD_CFLAGS3)
AM_CXXFLAGS = -D_REENTRANT $(ADD_CXXFLAGS1) $(ADD_CXXFLAGS2) $(ADD_CXXFLAGS3)  @exp_eigen3_include@
AM_CPPFLAGS = -I$(top_srcdir)/include

bin_PROGRAMS = testConfig simpleConfig readSL testBarrier

COMMON_LIBS = -lexputil $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)	\
$(FLIBS) @exp_vtk_ldlib@ @exp_vtk_lib@
EXPFLAGS = -L$(top_srcdir)/exputil/.libs
CXXLINK = mpiCC $(CXXFLAGS) $(EXPFLAGS) -o $@

BOOST_LIBS = $(BOOST_SYSTEM_LIB) $(BOOST_PROGRAM_OPTIONS_LIB)

noinst_HEADERS =

testConfig_SOURCES    = test_config.cc
testConfig_LDFLAGS    = $(BOOST_LDFLAGS)
testConfig_LDADD      = $(BOOST_LIBS)

simpleConfig_SOURCES  = simple_config.cc
simpleConfig_LDFLAGS  = $(BOOST_LDFLAGS)
simpleConfig_LDADD    = $(BOOST_LIBS)

readSL_SOURCES        = read_slsph.cc
readSL_LDFLAGS        = -L$(top_srcdir)/exputil/.libs
readSL_LDADD          = $(COMMON_LIBS) $(CUDA_LDFLAGS) $(CUDA_LIBFLAGS)

testBarrier_SOURCES   = test_barrier.cc
testBarrier_LDFLAGS = -L$(top_srcdir)/exputil/.libs $(BOOST_LDFLAGS)
testBarrier_LDADD = $(COMMON_LIBS) $(BOOST_LIBS) $(CUDA_LDFLAGS) $(CUDA_LIBFLAGS)

CLEANFILES	=

## utils/Test/Makefile.am ends here
