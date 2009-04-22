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

INCLUDES = -I$(top_srcdir)/include -I$(top_srcdir)/expgnu -I$(top_srcdir)/src

COMMON_LIBS = -lexputil -lexpgnu
SYS_LIBS = -L -lstdc++ -lpthread -ldl
MPI_LIBS = -lmpi

EXPLIBS = $(MPI_LIBS) $(SYS_LIBS) $(COMMON_LIBS) $(ADD_LIBS3) $(LIBS)
EXPFLAGS = -L$(top_srcdir)/exputil/.libs -L$(top_srcdir)/expgnu/.libs
CXXLINK = mpiCC $(CXXFLAGS) $(EXPFLAGS) -o $@

noinst_HEADERS =

bin_PROGRAMS = diskprof haloprof slabprof

diskprof_SOURCES =  PSP.cc header.cc diskprof.cc EmpOrth9thd.cc localmpi.cc Particle.cc
diskprof_LDADD = $(EXPLIBS)

haloprof_SOURCES =  PSP.cc header.cc haloprof.cc SphereSL.cc localmpi.cc Particle.cc
haloprof_CXXFLAGS = -I$(top_srcdir)/include
haloprof_LDADD = $(EXPLIBS)

slabprof_SOURCES =  slabprof.cc PSP.cc header.cc
slabprof_CXXFLAGS = -I$(top_srcdir)/include
slabprof_LDADD = $(COMMON_LIBS)

#	SphericalSL.cc legendre.cc ddplgndr.cc sinecosine.cc 

EmpOrth9thd.cc:	$(top_srcdir)/src/EmpOrth9thd.cc
	cp $(top_srcdir)/src/$@ .

PSP.cc:	$(top_srcdir)/utils/PhaseSpace/PSP.cc
	cp $(top_srcdir)/utils/PhaseSpace/PSP.cc .
	cp $(top_srcdir)/utils/PhaseSpace/PSP.H .

header.cc:	$(top_srcdir)/utils/PhaseSpace/header.cc
	cp $(top_srcdir)/utils/PhaseSpace/header.cc .
	cp $(top_srcdir)/utils/PhaseSpace/header.H .

localmpi.cc:	$(top_srcdir)/exputil/localmpi.cc
	cp $(top_srcdir)/exputil/$@ .

CLEANFILES = 	EmpOrth9thd.cc localmpi.cc PSP.cc PSP.H header.cc header.H

## utils/ICs/Makefile.am ends here


