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

bin_PROGRAMS = pspinfo pspstat psp2tipsy pspsplit psp2ascii ascii2psp \
		tipstd2psp psp2ta pspgas psp2histo psp2histo1d psp2histoE pspreal

COMMON_LIBS =  -lexputil -lexpgnu $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)
BOOST_LIBS = $(BOOST_SYSTEM_LIB) $(BOOST_PROGRAM_OPTIONS_LIB)

EXPLIBS = $(MPI_LIBS) $(COMMON_LIBS) $(ADD_LIBS3)
EXPFLAGS = -L$(top_srcdir)/exputil/.libs -L$(top_srcdir)/expgnu/.libs
CXXLINK = mpiCC $(CXXFLAGS) $(EXPFLAGS) -o $@

noinst_HEADERS =

pspinfo_SOURCES = pspinfo.cc PSP.cc header.cc header.H
pspinfo_LDADD   = $(COMMON_LIBS)

pspstat_SOURCES = pspstat.cc PSP.cc header.cc header.H
pspstat_LDADD   = $(COMMON_LIBS)

psp2tipsy_SOURCES = psp2tipsy.cc PSP.cc header.cc header.H
psp2tipsy_LDADD   = $(COMMON_LIBS)

psp2ta_SOURCES = psp2tarray.cc PSP.cc header.cc header.H
psp2ta_LDADD   = $(COMMON_LIBS)

pspsplit_SOURCES = pspsplit.cc PSP.cc header.cc header.H
pspsplit_LDADD   = $(COMMON_LIBS)

psp2ascii_SOURCES = psp2ascii.cc PSP.cc header.cc header.H
psp2ascii_LDADD   = $(COMMON_LIBS)

ascii2psp_SOURCES = ascii2psp.cc PSP.cc header.cc header.H
ascii2psp_LDADD   = $(COMMON_LIBS)

tipstd2psp_SOURCES = tipstd2psp.cc PSP.cc header.cc header.H tipsyxdr.c
tipstd2psp_LDADD   = $(COMMON_LIBS)

pspgas_SOURCES = pspprocess.cc PSP.cc header.cc header.H
pspgas_LDADD   = $(COMMON_LIBS)

psp2histo_SOURCES = psp2histo.cc PSP.cc header.cc header.H
psp2histo_LDADD   = $(COMMON_LIBS)

psp2histo1d_SOURCES  = psp2histo1d.cc PSP.cc header.cc header.H
psp2histo1d_CPPFLAGS = -I$(top_srcdir)/src
psp2histo1d_LDFLAGS  = $(BOOST_LDFLAGS)
psp2histo1d_LDADD    = $(COMMON_LIBS) $(BOOST_LIBS) $(EXPLIBS)

psp2histoE_SOURCES  = psp2histoE.cc PSP.cc header.cc header.H
psp2histoE_CPPFLAGS = -I$(top_srcdir)/src
psp2histoE_LDFLAGS  = $(BOOST_LDFLAGS)
psp2histoE_LDADD    = $(COMMON_LIBS) $(BOOST_LIBS) $(EXPLIBS)

pspreal_SOURCES = pspreal.cc PSP.cc header.cc header.H
pspreal_LDADD   = $(COMMON_LIBS)

*.cc:		header.H

header.H:	$(top_srcdir)/src/header.H
		cp $(top_srcdir)/src/$@ .

header.cc:	$(top_srcdir)/src/header.cc
		cp $(top_srcdir)/src/$@ .

CLEANFILES	= header.cc header.H

## utils/PhaseSpace/Makefile.am ends here
