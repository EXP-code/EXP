

## Process this file with automake to produce Makefile.in

CC = mpicc
CXX = mpiCC

# Allows check of libdirs for make dependence
VPATH = .

CXXFLAGS = -g -Wall -DDEBUG -D_REENTRANT

INCLUDES = -I$(top_srcdir)/include -I$(LAMHOME)/include -I$(LAMHOME)/include/mpi2c++ -I$(LAMHOME)/include/lam

bin_PROGRAMS = dfchk

COMMON_LIBS = -lnumerical -lg++ -lmdw
LAPACK_LIBS = -llapack -lblas -lf2c
FFTW_LIBS = -lrfftw -lfftw
SYS_LIBS = -lm
MPI_LIBS = -lmpi

LIBS = $(COMMON_LIBS) $(SYS_LIBS)
LDADD = $(COMMON_LIBS)
LDFLAGS = -L/home/weinberg/lib -L$(top_srcdir)/src  -L$(LAMHOME)/lib

noinst_HEADERS =

dfchk_SOURCES =  dfcheck2.cc genpoly.cc
dfchk_DEPENDENCIES = 
dfchk_LDADD = 


## Makefile.am ends here
