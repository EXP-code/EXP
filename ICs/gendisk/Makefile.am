

## Process this file with automake to produce Makefile.in

CC = mpicc
CXX = mpiCC
CXXLINK = mpif77 $(CXXFLAGS) $(LDFLAGS) -o $@

# Allows check of libdirs for make dependence
VPATH = .

CXXFLAGS = -g -Wall -DDEBUG -D_REENTRANT

INCLUDES = -I$(top_srcdir)/include -I$(LAMHOME)/include -I$(LAMHOME)/include/mpi2c++ -I$(LAMHOME)/include/lam

bin_PROGRAMS = gendisk

COMMON_LIBS = -lnumerical -lg++ -lmdw
LAPACK_LIBS = -llapack -lblas -lf2c
FFTW_LIBS = -lrfftw -lfftw
SYS_LIBS = -lm
MPI_LIBS = -lmpi

LIBS = $(COMMON_LIBS) $(SYS_LIBS)
LDADD = $(COMMON_LIBS)
LDFLAGS = -L$(top_srcdir)/src  -L$(LAMHOME)/lib

noinst_HEADERS =

gendisk_SOURCES =  initial.cc DiskHalo2.cc SLGridMP2.cc EmpOrth9thd.cc \
	SphericalSL.cc localmpi.cc sledge.f legendre.cc \
	ddplgndr.cc sinecosine.cc AddDisk.cc
gendisk_DEPENDENCIES = 
gendisk_LDADD = 

## Makefile.am ends here
