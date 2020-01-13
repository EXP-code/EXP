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
AM_CXXFLAGS = -std=c++11 -D_REENTRANT -DSTANDALONE $(ADD_CXXFLAGS1) $(ADD_CXXFLAGS2) $(ADD_CXXFLAGS3) @exp_vtk_include@
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/expgnu

bin_PROGRAMS = pspinfo pspstat psp2tipsy pspsplit psp2ascii ascii2psp	\
		tipstd2psp psp2ta pspgas psp2histo psp2histo1d		\
		psp2histo1D psp2histoE psp2histoH pspreal pspstatS	\
		psp2range psp2eHisto psp2histoT psp2histoTC pspstatT	\
		spl2psp pspinfo2 ascii2psp2

COMMON_LIBS = -lexputil -lexpgnu $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS) @exp_vtk_ldlib@ @exp_vtk_lib@ $(YAML)

BOOST_LIBS = $(BOOST_SYSTEM_LIB) $(BOOST_PROGRAM_OPTIONS_LIB)

EXPLIBS = $(MPI_LIBS) $(COMMON_LIBS) $(ADD_LIBS3)
EXPFLAGS = -L$(top_srcdir)/exputil/.libs -L$(top_srcdir)/expgnu/.libs
CXXLINK = mpiCC -rdynamic $(CXXFLAGS) $(EXPFLAGS) -o $@

noinst_HEADERS =

pspinfo_SOURCES = pspinfo.cc PSP.cc header.cc header.H
pspinfo_LDADD   = $(COMMON_LIBS)

pspinfo2_SOURCES = pspinfo2.cc PSP2.cc header.cc header.H
pspinfo2_LDADD   = $(COMMON_LIBS)

pspstat_SOURCES  = pspstat.cc PSP2.cc header.cc header.H
pspstat_LDADD    = $(COMMON_LIBS)

pspstatS_SOURCES  = pspstatS.cc PSP.cc header.cc header.H
pspstatS_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/include
pspstatS_LDFLAGS  = $(BOOST_LDFLAGS) $(EXPFLAGS)
pspstatS_LDADD    = $(COMMON_LIBS) $(BOOST_LIBS)

pspstatT_SOURCES  = pspstatT.cc PSP2.cc header.cc header.H
pspstatT_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/include
pspstatT_LDFLAGS  = $(BOOST_LDFLAGS) $(EXPFLAGS)
pspstatT_LDADD    = $(COMMON_LIBS) $(BOOST_LIBS)

psp2tipsy_SOURCES = psp2tipsy.cc PSP.cc header.cc header.H
psp2tipsy_LDADD   = $(COMMON_LIBS)

psp2ta_SOURCES = psp2tarray.cc PSP.cc header.cc header.H
psp2ta_LDADD   = $(COMMON_LIBS)

pspsplit_SOURCES = pspsplit.cc PSP.cc header.cc header.H
pspsplit_LDADD   = $(COMMON_LIBS) $(EXPFLAGS)

psp2ascii_SOURCES = psp2ascii.cc PSP2.cc header.cc header.H
psp2ascii_LDADD   = $(COMMON_LIBS)

ascii2psp_SOURCES = ascii2psp.cc PSP.cc header.cc header.H
ascii2psp_LDADD   = $(COMMON_LIBS)

ascii2psp2_SOURCES = ascii2psp2.cc PSP2.cc header.cc header.H
ascii2psp2_LDADD   = $(COMMON_LIBS)

tipstd2psp_SOURCES = tipstd2psp.cc PSP.cc header.cc header.H tipsyxdr.c
tipstd2psp_LDADD   = $(COMMON_LIBS)

pspgas_SOURCES = pspprocess.cc PSP.cc header.cc header.H
pspgas_LDADD   = $(COMMON_LIBS)

psp2histo_SOURCES = psp2histo.cc PSP2.cc header.cc header.H
psp2histo_LDADD   = $(COMMON_LIBS)

psp2histo1D_SOURCES  = psp2histo1D.cc PSP2.cc header.cc header.H
psp2histo1D_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/include
psp2histo1D_LDFLAGS  = $(BOOST_LDFLAGS) $(EXPFLAGS)
psp2histo1D_LDADD    = $(COMMON_LIBS) $(BOOST_LIBS) $(EXPLIBS)

psp2histo1d_SOURCES  = psp2histo1d.cc PSP.cc header.cc header.H
psp2histo1d_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/include
psp2histo1d_LDFLAGS  = $(BOOST_LDFLAGS) $(EXPFLAGS)
psp2histo1d_LDADD    = $(COMMON_LIBS) $(BOOST_LIBS) $(EXPLIBS)

psp2histoE_SOURCES  = psp2histoE.cc PSP.cc header.cc header.H
psp2histoE_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/include
psp2histoE_LDFLAGS  = $(BOOST_LDFLAGS) $(EXPFLAGS)
psp2histoE_LDADD    = $(COMMON_LIBS) $(BOOST_LIBS) $(EXPLIBS)

psp2histoH_SOURCES  = psp2histoH.cc PSP.cc header.cc header.H
psp2histoH_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/include
psp2histoH_LDFLAGS  = $(BOOST_LDFLAGS) $(EXPFLAGS)
psp2histoH_LDADD    = $(COMMON_LIBS) $(BOOST_LIBS) $(EXPLIBS)

psp2histoT_SOURCES  = psp2histoT.cc PSP.cc header.cc header.H
psp2histoT_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/include
psp2histoT_LDFLAGS  = $(BOOST_LDFLAGS) $(EXPFLAGS)
psp2histoT_LDADD    = $(COMMON_LIBS) $(BOOST_LIBS) $(EXPLIBS)

psp2histoTC_SOURCES  = psp2histoTC.cc PSP.cc header.cc header.H
psp2histoTC_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/include
psp2histoTC_LDFLAGS  = $(BOOST_LDFLAGS) $(EXPFLAGS)
psp2histoTC_LDADD    = $(COMMON_LIBS) $(BOOST_LIBS) $(EXPLIBS)

psp2eHisto_SOURCES  = psp2eHisto.cc PSP.cc header.cc header.H
psp2eHisto_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/include
psp2eHisto_LDFLAGS  = $(BOOST_LDFLAGS) $(EXPFLAGS)
psp2eHisto_LDADD    = $(COMMON_LIBS) $(BOOST_LIBS) $(EXPLIBS)
pspreal_SOURCES = pspreal.cc PSP.cc header.cc header.H
pspreal_LDADD   = $(COMMON_LIBS)

psp2range_SOURCES = psp2range.cc PSP2.cc header.cc header.H
psp2range_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/include
psp2range_LDFLAGS  = $(BOOST_LDFLAGS) $(EXPFLAGS)
psp2range_LDADD   = $(COMMON_LIBS) $(BOOST_LIBS)

spl2psp_SOURCES  = spl2psp.cc PSP.cc header.cc header.H Particle.cc
spl2psp_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/include -I$(top_srcdir)/expgnu
spl2psp_LDFLAGS  = $(BOOST_LDFLAGS) $(EXPFLAGS)
spl2psp_LDADD    = $(COMMON_LIBS) $(BOOST_LIBS)

*.cc:		header.H

header.H:	$(top_srcdir)/src/header.H
		cp $(top_srcdir)/src/$@ .

header.cc:	$(top_srcdir)/src/header.cc
		cp $(top_srcdir)/src/$@ .

Particle.cc:	$(top_srcdir)/src/Particle.cc
		cp $(top_srcdir)/src/$@ .

CLEANFILES	= header.cc header.H Particle.cc

## utils/PhaseSpace/Makefile.am ends here
