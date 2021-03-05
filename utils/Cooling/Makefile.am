## Process this file with automake to produce Makefile.in

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
AM_CXXFLAGS = -D_REENTRANT -DSTANDALONE $(ADD_CXXFLAGS1) $(ADD_CXXFLAGS2) $(ADD_CXXFLAGS3)

AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/expgnu -I$(top_srcdir)/src -I$(top_srcdir)/src/DSMC

COMMON_LIBS = -lexputil -lexpgnu $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)	\
$(FLIBS) $(CUDA_LDFLAGS) $(CUDA_LIBFLAGS)
MPI_LIBS = -lmpi

EXPLIBS = $(MPI_LIBS) $(COMMON_LIBS) $(ADD_LIBS3)
EXPFLAGS = -L$(top_srcdir)/exputil/.libs -L$(top_srcdir)/expgnu/.libs
CXXLINK = mpiCC $(CXXFLAGS) $(EXPFLAGS) -o $@

noinst_HEADERS =

bin_PROGRAMS = hctest

hctest_SOURCES =  hctest.cc HeatCool.cc
hctest_LDADD   = $(EXPLIBS) @exp_vtk_ldlib@ @exp_vtk_lib@

HeatCool.cc:	$(top_srcdir)/src/DSMC/HeatCool.cc
		cp $(top_srcdir)/src/DSMC/$@ .

CLEANFILES = 	HeatCool.cc

## utils/Cooling/Makefile.am ends here


