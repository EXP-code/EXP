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

bin_PROGRAMS = testparse testparams testparamsMPI testConfig	\
simpleConfig readSL testBarrier

COMMON_LIBS = -lexputil -lexpgnu $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)
EXPFLAGS = -L$(top_srcdir)/exputil/.libs -L$(top_srcdir)/expgnu/.libs
CXXLINK = mpiCC $(CXXFLAGS) $(EXPFLAGS) -o $@

BOOST_LIBS = $(BOOST_SYSTEM_LIB) $(BOOST_PROGRAM_OPTIONS_LIB)

noinst_HEADERS =

testparse_SOURCES = testparse.cc
testparse_LDADD   = $(COMMON_LIBS)

testparams_SOURCES = testparams.cc
testparams_LDADD   = $(COMMON_LIBS)

testparamsMPI_SOURCES = testparamsMPI.cc
testparamsMPI_LDADD   = $(COMMON_LIBS)

testConfig_SOURCES    = test_config.cc
testConfig_LDFLAGS    = $(BOOST_LDFLAGS)
testConfig_LDADD      = $(BOOST_LIBS)

simpleConfig_SOURCES  = simple_config.cc
simpleConfig_LDFLAGS  = $(BOOST_LDFLAGS)
simpleConfig_LDADD    = $(BOOST_LIBS)

readSL_SOURCES        = read_slsph.cc
readSL_LDFLAGS        = -L$(top_srcdir)/exputil/.libs	\
			-L$(top_srcdir)/expgnu/.libs
readSL_LDADD          = $(COMMON_LIBS)

testBarrier_SOURCES   = test_barrier.cc
testBarrier_LDFLAGS   = -L$(top_srcdir)/exputil/.libs	\
			-L$(top_srcdir)/expgnu/.libs $(BOOST_LDFLAGS)
testBarrier_LDADD     = $(COMMON_LIBS) $(BOOST_LIBS)

CLEANFILES	=

## utils/Test/Makefile.am ends here
