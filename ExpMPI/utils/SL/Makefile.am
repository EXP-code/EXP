## Process this file with automake to produce Makefile.in

CC = mpicc
CXX = mpiCC

if ENABLE_DMALLOC
ADD_CFLAGS3 = -DUSE_DMALLOC -g
ADD_CXXFLAGS3 = -DUSE_DMALLOC -g
ADD_LIBS3 = -ldmallocth -ldmallocthcxx
endif

AM_CFLAGS = -D_REENTRANT $(ADD_CFLAGS1) $(ADD_CFLAGS2) $(ADD_CFLAGS3)
AM_CXXFLAGS = -D_REENTRANT $(ADD_CXXFLAGS1) $(ADD_CXXFLAGS2) $(ADD_CXXFLAGS3)

INCLUDES = -I$(top_srcdir)/include -I$(top_srcdir)/expgnu -I$(top_srcdir)/src

bin_PROGRAMS = slcheck slshift orthochk

COMMON_LIBS = -lexputil -lexpgnu
SYS_LIBS = -lstdc++ -lpthread -ldl -lg2c -lm
MPI_LIBS = -lmpi

LIBS = $(MPI_LIBS) $(SYS_LIBS) $(COMMON_LIBS) $(ADD_LIBS3)
EXPFLAGS = -L$(top_srcdir)/exputil/.libs -L$(top_srcdir)/expgnu/.libs
CXXLINK = mpiCC $(CXXFLAGS) $(EXPFLAGS) -o $@

noinst_HEADERS =

slcheck_SOURCES = slcheck.cc

slshift_SOURCES = slshift.cc SLSphere.cc

orthochk_SOURCES = orthochk.cc biorth1d.cc

biorth1d.h:	$(top_srcdir)/src/biorth1d.h
	cp $(top_srcdir)/src/$@ .

biorth1d.cc:	$(top_srcdir)/src/biorth1d.cc
	cp $(top_srcdir)/src/$@ .

CLEANFILES = 	biorth1d.h biorth1d.cc 

## utils/SL/Makefile.am ends here
