

## Process this file with automake to produce Makefile.in

CC = $(LAMHOME)/bin/hcc
CXX = $(LAMHOME)/bin/hcp

# Allows check of libdirs for make dependence
VPATH = .

CXXFLAGS = -g -Wall -DDEBUG -D_REENTRANT

INCLUDES = -I$(top_srcdir)/include -I$(LAMHOME)/include -I$(LAMHOME)/include/mpi2c++ -I$(LAMHOME)/include/lam

bin_PROGRAMS = gensph

COMMON_LIBS = -lnumerical -lg++ -lmdw
LAPACK_LIBS = -llapack -lblas -lf2c
FFTW_LIBS = -lrfftw -lfftw
SYS_LIBS = -lm
MPI_LIBS = -lmpi -ltrillium -ltstdio -lt -largs

LIBS = $(COMMON_LIBS) $(SYS_LIBS)
LDADD = $(COMMON_LIBS)
LDFLAGS = -L$(top_srcdir)/src  -L$(LAMHOME)/lib

noinst_HEADERS =

gensph_SOURCES =  gensph.cc genpoly.cc
gensph_DEPENDENCIES = 
gensph_LDADD = 


## Makefile.am ends here
