## Process this file with automake to produce Makefile.in

# VPATH = $(top_srcdir)/src:.



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

INCLUDES = -I$(top_srcdir)/include -I$(top_srcdir)/expgnu

bin_PROGRAMS = pspinfo psp2tipsy pspsplit psp2ascii

COMMON_LIBS = -lexputil -lexpgnu
SYS_LIBS = -lstdc++ -lpthread -ldl -lm
MPI_LIBS = -lmpi

LIBS = $(MPI_LIBS) $(SYS_LIBS) $(COMMON_LIBS) $(ADD_LIBS3)
EXPFLAGS = -L$(top_srcdir)/exputil/.libs -L$(top_srcdir)/expgnu/.libs
CXXLINK = mpiCC $(CXXFLAGS) $(EXPFLAGS) -o $@

noinst_HEADERS =

pspinfo_SOURCES = pspinfo.cc header.cc

psp2tipsy_SOURCES = psp2tipsy.cc header.cc

pspsplit_SOURCES = pspsplit.cc header.cc

psp2ascii_SOURCES = psp2ascii.cc header.cc

## src/ParseTest/Makefile.am ends here







