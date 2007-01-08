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
SYS_LIBS = -L -L/usr/lib/gcc/i486-linux-gnu/4.0.1 -lstdc++ -lpthread -ldl -lg2c -lm
MPI_LIBS = -lmpi

LIBS = $(MPI_LIBS) $(SYS_LIBS) $(COMMON_LIBS) $(ADD_LIBS3)
EXPFLAGS = -L$(top_srcdir)/exputil/.libs -L$(top_srcdir)/expgnu/.libs
CXXLINK = mpiCC $(CXXFLAGS) $(EXPFLAGS) -o $@

noinst_HEADERS =

bin_PROGRAMS = gendisk

gendisk_SOURCES =  initial.cc DiskHalo.cc EmpOrth9thd.cc \
	SphericalSL.cc legendre.cc Particle.cc \
	ddplgndr.cc sinecosine.cc AddDisk.cc localmpi.cc

Particle.cc:	$(top_srcdir)/src/Particle.cc
	cp $(top_srcdir)/src/$@ .

EmpOrth9thd.cc:	$(top_srcdir)/src/EmpOrth9thd.cc
	cp $(top_srcdir)/src/$@ .

CLEANFILES = 	Particle.cc EmpOrth9thd.cc 

## utils/ICs/Makefile.am ends here


