## Process this file with automake to produce Makefile.in

CC = mpicc
CXX = mpiCC

if ENABLE_DMALLOC
ADD_CFLAGS3 = -DUSE_DMALLOC -g
ADD_CXXFLAGS3 = -DUSE_DMALLOC -g
ADD_LIBS3 = -ldmallocth -ldmallocthcxx
endif

AM_CFLAGS = -D_REENTRANT $(ADD_CFLAGS1) $(ADD_CFLAGS2) $(ADD_CFLAGS3)
AM_CXXFLAGS = -D_REENTRANT -DSTANDALONE $(ADD_CXXFLAGS1) $(ADD_CXXFLAGS2) $(ADD_CXXFLAGS3)

INCLUDES = -I$(top_srcdir)/include -I$(top_srcdir)/expgnu -I$(top_srcdir)/src -I$(top_srcdir)/src/user

COMMON_LIBS = -lexputil -lexpgnu $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)
MPI_LIBS = -lmpi

EXPLIBS = $(MPI_LIBS) $(COMMON_LIBS) $(ADD_LIBS3)
EXPFLAGS = -L$(top_srcdir)/exputil/.libs -L$(top_srcdir)/expgnu/.libs
CXXLINK = mpiCC $(CXXFLAGS) $(EXPFLAGS) -o $@

noinst_HEADERS =

bin_PROGRAMS = hctest

hctest_SOURCES =  hctest.cc HeatCool.cc
hctest_LDADD   = $(EXPLIBS)

HeatCool.cc:	$(top_srcdir)/src/user/HeatCool.cc
		cp $(top_srcdir)/src/user/$@ .

CLEANFILES = 	HeatCool.cc

## utils/Cooling/Makefile.am ends here


