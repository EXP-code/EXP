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
AM_CXXFLAGS = -D_REENTRANT $(ADD_CXXFLAGS1) $(ADD_CXXFLAGS2) $(ADD_CXXFLAGS3)

INCLUDES = -I$(top_srcdir)/include -I$(top_srcdir)/expgnu

bin_PROGRAMS = pspinfo psp2tipsy pspsplit psp2ascii ascii2psp

COMMON_LIBS = -lexputil -lexpgnu
SYS_LIBS = -lstdc++ -lpthread -ldl -lm
MPI_LIBS = -lmpi

LIBS = $(MPI_LIBS) $(SYS_LIBS) $(COMMON_LIBS) $(ADD_LIBS3)
EXPFLAGS = -L$(top_srcdir)/exputil/.libs -L$(top_srcdir)/expgnu/.libs
CXXLINK = mpiCC $(CXXFLAGS) $(EXPFLAGS) -o $@

noinst_HEADERS =

pspinfo_SOURCES = pspinfo.cc PSP.cc header.cc

psp2tipsy_SOURCES = psp2tipsy.cc PSP.cc header.cc

pspsplit_SOURCES = pspsplit.cc PSP.cc header.cc

psp2ascii_SOURCES = psp2ascii.cc PSP.cc header.cc

ascii2psp_SOURCES = ascii2psp.cc PSP.cc header.cc

*.cc:		header.H

header.H:	$(top_srcdir)/src/header.H
		cp $(top_srcdir)/src/$@ .

header.cc:	$(top_srcdir)/src/header.cc
		cp $(top_srcdir)/src/$@ .

CLEANFILES	= header.cc header.H

## utils/PhaseSpace/Makefile.am ends here
