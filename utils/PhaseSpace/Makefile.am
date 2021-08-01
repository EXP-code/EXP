## Process this file with automake to produce Makefile.in

# VPATH = $(top_srcdir)/src:.

CC = mpicc
CXX = mpiCC

if ENABLE_DMALLOC
ADD_CFLAGS3 = -DUSE_DMALLOC -g
ADD_CXXFLAGS3 = -DUSE_DMALLOC -g
ADD_LIBS3 = -ldmallocth -ldmallocthcxx
endif

# Need for lib exputil which contains cuda code
#
if ENABLE_CUDA
  CUDA_LIBFLAGS  = -lcudart -lnvToolsExt
endif

AM_CFLAGS = -D_REENTRANT $(ADD_CFLAGS1) $(ADD_CFLAGS2) $(ADD_CFLAGS3)
AM_CXXFLAGS = -D_REENTRANT -DSTANDALONE $(ADD_CXXFLAGS1) $(ADD_CXXFLAGS2) $(ADD_CXXFLAGS3) @exp_vtk_include@
AM_CPPFLAGS = -I$(top_srcdir)/include 

bin_PROGRAMS = pspinfo pspstat psp2ascii ascii2psp tipstd2psp		\
		psp2histo psp2histo1d psp2histo1D2 psp2histoE		\
		psp2histoH pspreal pspstatS psp2range psp2eHisto	\
		psp2histoT psp2histoTC pspstatT spl2psp ascii2psp2	\
		psp2vtu psp2rings psp2bess psp2lagu psp2interp

COMMON_LIBS = -lexputil -lexp $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)	\
	$(FLIBS) @exp_vtk_ldlib@ @exp_vtk_lib@ $(YAML)			\
	$(CUDA_LDFLAGS) $(CUDA_LIBFLAGS)

BOOST_LIBS = $(BOOST_SYSTEM_LIB) $(BOOST_RANDOM_LIB) $(BOOST_PROGRAM_OPTIONS_LIB)

EXPLIBS = $(MPI_LIBS) $(COMMON_LIBS) $(ADD_LIBS3)
EXPFLAGS = -L$(top_srcdir)/exputil/.libs -L$(top_srcdir)/src/.libs
CXXLINK = mpiCC -rdynamic $(CXXFLAGS) $(EXPFLAGS) -o $@

noinst_HEADERS  = Grid2D.H

pspinfo_SOURCES = pspinfo.cc PSP.cc
pspinfo_LDADD   = $(COMMON_LIBS)

pspstat_SOURCES  = pspstat.cc PSP.cc
pspstat_LDADD    = $(COMMON_LIBS)

pspstatS_SOURCES  = pspstatS.cc PSP.cc
pspstatS_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/include
pspstatS_LDFLAGS  = $(BOOST_LDFLAGS) $(EXPFLAGS)
pspstatS_LDADD    = $(COMMON_LIBS) $(BOOST_LIBS)

pspstatT_SOURCES  = pspstatT.cc PSP.cc
pspstatT_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/include
pspstatT_LDFLAGS  = $(BOOST_LDFLAGS) $(EXPFLAGS)
pspstatT_LDADD    = $(COMMON_LIBS) $(BOOST_LIBS)

psp2ascii_SOURCES = psp2ascii.cc PSP.cc
psp2ascii_LDADD   = $(COMMON_LIBS)

ascii2psp_SOURCES = ascii2psp.cc PSP.cc
ascii2psp_LDADD   = $(COMMON_LIBS)

ascii2psp2_SOURCES = ascii2psp2.cc PSP.cc
ascii2psp2_LDADD   = $(COMMON_LIBS)

tipstd2psp_SOURCES = tipstd2psp.cc PSP.cc tipsyxdr.c
tipstd2psp_LDADD   = $(COMMON_LIBS)

psp2histo_SOURCES = psp2histo.cc PSP.cc
psp2histo_LDADD   = $(COMMON_LIBS)

psp2histo1D2_SOURCES  = psp2histo1D_2.cc PSP.cc
psp2histo1D2_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/include
psp2histo1D2_LDFLAGS  = $(BOOST_LDFLAGS) $(EXPFLAGS)
psp2histo1D2_LDADD    = $(COMMON_LIBS) $(BOOST_LIBS) $(EXPLIBS)

psp2histo1d_SOURCES  = psp2histo1d.cc PSP.cc
psp2histo1d_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/include
psp2histo1d_LDFLAGS  = $(BOOST_LDFLAGS) $(EXPFLAGS)
psp2histo1d_LDADD    = $(COMMON_LIBS) $(BOOST_LIBS) $(EXPLIBS)

psp2histoE_SOURCES  = psp2histoE.cc PSP.cc
psp2histoE_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/include
psp2histoE_LDFLAGS  = $(BOOST_LDFLAGS) $(EXPFLAGS)
psp2histoE_LDADD    = $(COMMON_LIBS) $(BOOST_LIBS) $(EXPLIBS)

psp2histoH_SOURCES  = psp2histoH.cc PSP.cc
psp2histoH_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/include
psp2histoH_LDFLAGS  = $(BOOST_LDFLAGS) $(EXPFLAGS)
psp2histoH_LDADD    = $(COMMON_LIBS) $(BOOST_LIBS) $(EXPLIBS)

psp2histoT_SOURCES  = psp2histoT.cc PSP.cc
psp2histoT_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/include
psp2histoT_LDFLAGS  = $(BOOST_LDFLAGS) $(EXPFLAGS)
psp2histoT_LDADD    = $(COMMON_LIBS) $(BOOST_LIBS) $(EXPLIBS)

psp2histoTC_SOURCES  = psp2histoTC.cc PSP.cc
psp2histoTC_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/include
psp2histoTC_LDFLAGS  = $(BOOST_LDFLAGS) $(EXPFLAGS)
psp2histoTC_LDADD    = $(COMMON_LIBS) $(BOOST_LIBS) $(EXPLIBS)

psp2eHisto_SOURCES  = psp2eHisto.cc PSP.cc
psp2eHisto_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/include
psp2eHisto_LDFLAGS  = $(BOOST_LDFLAGS) $(EXPFLAGS)
psp2eHisto_LDADD    = $(COMMON_LIBS) $(BOOST_LIBS) $(EXPLIBS)
pspreal_SOURCES = pspreal.cc PSP.cc
pspreal_LDADD   = $(COMMON_LIBS)

psp2range_SOURCES = psp2range.cc PSP.cc
psp2range_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/include
psp2range_LDFLAGS  = $(BOOST_LDFLAGS) $(EXPFLAGS)
psp2range_LDADD   = $(COMMON_LIBS) $(BOOST_LIBS)

spl2psp_SOURCES  = spl2psp.cc PSP.cc
spl2psp_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/include
spl2psp_LDFLAGS  = $(BOOST_LDFLAGS) $(EXPFLAGS)
spl2psp_LDADD    = $(COMMON_LIBS) $(BOOST_LIBS)

psp2vtu_SOURCES  = psp2vtu.cc PSP.cc
psp2vtu_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/include
psp2vtu_LDFLAGS  = $(BOOST_LDFLAGS) $(EXPFLAGS)
psp2vtu_LDADD    = $(COMMON_LIBS) $(BOOST_LIBS)

psp2rings_SOURCES  = psp2rings.cc PSP.cc
psp2rings_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/include
psp2rings_LDFLAGS  = $(BOOST_LDFLAGS) $(EXPFLAGS)
psp2rings_LDADD    = $(COMMON_LIBS) $(BOOST_LIBS) $(EXPLIBS)

psp2bess_SOURCES  = psp2bess.cc PSP.cc
psp2bess_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/include
psp2bess_LDFLAGS  = $(BOOST_LDFLAGS) $(EXPFLAGS)
psp2bess_LDADD    = $(COMMON_LIBS) $(BOOST_LIBS) $(EXPLIBS)

psp2lagu_SOURCES  = psp2lagu.cc PSP.cc
psp2lagu_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/include
psp2lagu_LDFLAGS  = $(BOOST_LDFLAGS) $(EXPFLAGS)
psp2lagu_LDADD    = $(COMMON_LIBS) $(BOOST_LIBS) $(EXPLIBS)

psp2interp_SOURCES  = psp2interp.cc PSP.cc
psp2interp_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/include
psp2interp_LDFLAGS  = $(BOOST_LDFLAGS) $(EXPFLAGS)
psp2interp_LDADD    = $(COMMON_LIBS) $(BOOST_LIBS) $(EXPLIBS)

CLEANFILES =

## utils/PhaseSpace/Makefile.am ends here
