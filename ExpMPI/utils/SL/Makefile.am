## Process this file with automake to produce Makefile.in

CC = mpicc
CXX = mpiCC

CFLAGS = -g
CXXFLAGS = -g

if ENABLE_DMALLOC
ADD_CFLAGS3 = -DUSE_DMALLOC -g
ADD_CXXFLAGS3 = -DUSE_DMALLOC -g
ADD_LIBS3 = -ldmallocth -ldmallocthcxx
endif

AM_CFLAGS = -D_REENTRANT $(ADD_CFLAGS1) $(ADD_CFLAGS2) $(ADD_CFLAGS3)
AM_CXXFLAGS = -D_REENTRANT $(ADD_CXXFLAGS1) $(ADD_CXXFLAGS2) $(ADD_CXXFLAGS3)

CXXLINK = mpiCC $(CXXFLAGS) $(LDFLAGS) -o $@

INCLUDES = -I$(top_srcdir)/include -I$(top_srcdir)/expgnu -I$(top_srcdir)/src

bin_PROGRAMS = slcheck slcheckN

COMMON_LIBS = -lexputil -lexpgnu
SYS_LIBS = -lstdc++ -lpthread -ldl -lm
MPI_LIBS = -lmpi

LIBS = $(MPI_LIBS) $(SYS_LIBS) $(COMMON_LIBS) $(ADD_LIBS3)
EXPFLAGS = -L$(top_srcdir)/exputil/.libs -L$(top_srcdir)/expgnu/.libs
CXXLINK = mpif77 $(CXXFLAGS) $(EXPFLAGS) -o $@

noinst_HEADERS =

slcheck_SOURCES = slcheck.cc SLGridMP2.cc sledge.f
slcheckN_SOURCES = slcheck_nompi.cc SLGridMP2.cc sledge.f

## Copy files

SLGridMP2.h:	$(top_srcdir)/src/SLGridMP2.h
	cp $(top_srcdir)/src/$@ .

SLGridMP2.cc:	$(top_srcdir)/src/SLGridMP2.cc
	cp $(top_srcdir)/src/$@ .

sledge.f:	$(top_srcdir)/src/sledge.f
	cp $(top_srcdir)/src/$@ .

CLEANFILES = 	SLGridMP2.h SLGridMP2.cc sledge.f

## utils/SL/Makefile.am ends here
