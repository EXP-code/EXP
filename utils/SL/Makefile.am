## Process this file with automake to produce Makefile.in

CC  = mpicc
CXX = mpiCC

if ENABLE_DMALLOC
ADD_CFLAGS3 = -DUSE_DMALLOC -g
ADD_CXXFLAGS3 = -DUSE_DMALLOC -g
ADD_LIBS3 = -ldmallocth -ldmallocthcxx
endif

AM_CFLAGS = -D_REENTRANT $(ADD_CFLAGS1) $(ADD_CFLAGS2) $(ADD_CFLAGS3)
AM_CXXFLAGS = -D_REENTRANT $(ADD_CXXFLAGS1) $(ADD_CXXFLAGS2) $(ADD_CXXFLAGS3) $(BOOST_CPPFLAGS)
AM_CPPFLAGS= -I$(top_srcdir)/include -I$(top_srcdir)/expgnu -I$(top_srcdir)/src $(CUDA_CFLAGS)

if ENABLE_CUDA
CUDA_LIBFLAGS  = -lcudart
endif

bin_PROGRAMS = slcheck slshift orthochk diskpot

COMMON_LIBS = -lexputil -lexpgnu $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS)	\
$(FLIBS) @exp_vtk_ldlib@ @exp_vtk_lib@  $(CUDA_LDFLAGS)
SYS_LIBS = -lstdc++ -lpthread -ldl $(LIBS)
MPI_LIBS = -lmpi

BOOST_LIBS = $(BOOST_SYSTEM_LIB) $(BOOST_FILESYSTEM_LIB) $(BOOST_PROGRAM_OPTIONS_LIB)

EXPLIBS = $(MPI_LIBS) $(SYS_LIBS) $(COMMON_LIBS) $(ADD_LIBS3)
EXPFLAGS = -L$(top_srcdir)/exputil/.libs -L$(top_srcdir)/expgnu/.libs
CXXLINK = mpiCC $(CXXFLAGS) $(EXPFLAGS) -o $@

noinst_HEADERS =

slcheck_SOURCES = slcheck.cc localmpi.cc
slcheck_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS)
slcheck_LDADD = $(EXPLIBS) $(CUDA_LIBFLAGS) $(BOOST_LIBS)

slshift_SOURCES = slshift.cc SLSphere.cc localmpi.cc
slshift_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS)
slshift_LDADD = $(EXPLIBS) $(CUDA_LIBFLAGS) $(BOOST_LIBS)

orthochk_SOURCES = orthochk.cc biorth1d.cc localmpi.cc
orthochk_LDADD = $(EXPLIBS) $(CUDA_LIBFLAGS)

diskpot_SOURCES = diskpot.cc CylindricalDisk.cc SLSphere.cc localmpi.cc
diskpot_LDADD = $(EXPLIBS) $(CUDA_LIBFLAGS)

biorth1d.h:	$(top_srcdir)/src/biorth1d.h
	cp $(top_srcdir)/src/$@ .

biorth1d.cc:	$(top_srcdir)/src/biorth1d.cc
	cp $(top_srcdir)/src/$@ .

localmpi.cc:	$(top_srcdir)/exputil/localmpi.cc
	cp $(top_srcdir)/exputil/$@ .

CLEANFILES = 	biorth1d.h biorth1d.cc localmpi.cc

## utils/SL/Makefile.am ends here
