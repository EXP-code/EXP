## Process this file with automake to produce Makefile.in

# VPATH = $(top_srcdir)/src:.

CC = mpicc
CXX = mpiCC

if ENABLE_DMALLOC
ADD_CFLAGS3 = -DUSE_DMALLOC -g
ADD_CXXFLAGS3 = -DUSE_DMALLOC -g
ADD_LIBS3 = -ldmallocth -ldmallocthcxx
endif

AM_CFLAGS = -D_REENTRANT $(ADD_CFLAGS1) $(ADD_CFLAGS2) $(ADD_CFLAGS3)
AM_CXXFLAGS = -std=c++11 -D_REENTRANT $(ADD_CXXFLAGS1) $(ADD_CXXFLAGS2) $(ADD_CXXFLAGS3)
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/expgnu

bin_PROGRAMS = testparse testparams testparamsMPI

COMMON_LIBS = -lexputil -lexpgnu $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)
EXPFLAGS = -L$(top_srcdir)/exputil/.libs -L$(top_srcdir)/expgnu/.libs
CXXLINK = mpiCC $(CXXFLAGS) $(EXPFLAGS) -o $@

noinst_HEADERS =

testparse_SOURCES = testparse.cc
testparse_LDADD   = $(COMMON_LIBS)

testparams_SOURCES = testparams.cc
testparams_LDADD   = $(COMMON_LIBS)

testparamsMPI_SOURCES = testparamsMPI.cc
testparamsMPI_LDADD   = $(COMMON_LIBS)

CLEANFILES	=

## utils/Test/Makefile.am ends here
