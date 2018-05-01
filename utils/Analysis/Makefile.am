## Process this file with automake to produce Makefile.in

CC  = mpicc
CXX = mpiCC

if ENABLE_DMALLOC
ADD_CFLAGS3   = -DUSE_DMALLOC -g
ADD_CXXFLAGS3 = -DUSE_DMALLOC -g
ADD_LIBS3     = -ldmallocth -ldmallocthcxx
endif

AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/expgnu -I$(top_srcdir)/src
AM_CFLAGS = -D_REENTRANT $(ADD_CFLAGS3)
AM_CXXFLAGS = -std=c++11 -D_REENTRANT -DSTANDALONE $(INCLUDES) $(ADD_CXXFLAGS3) $(BOOST_CPPFLAGS) @exp_vtk_include@

COMMON_LIBS = $(LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) @exp_vtk_ldlib@ @exp_vtk_lib@ 
MPI_LIBS = -lmpi
EXPLIBS = -lexputil -lexpgnu $(COMMON_LIBS) $(ADD_LIBS3) $(MPI_LIBS)
EXPFLAGS = -L$(top_srcdir)/exputil/.libs -L$(top_srcdir)/expgnu/.libs
CXXLINK = mpiCC $(CXXFLAGS) $(EXPFLAGS) -o $@

noinst_HEADERS =

bin_PROGRAMS = diskprof haloprof slabprof psphisto gas2d gas2dcyl makegrid pspbox eofpeek diskprof2 pcatest

diskprof_SOURCES = PSP.cc header.cc diskprof.cc EmpOrth9thd.cc localmpi.cc Particle.cc
diskprof_LDADD = $(EXPLIBS) @exp_vtk_ldlib@ @exp_vtk_lib@

haloprof_SOURCES = PSP.cc header.cc haloprof.cc SphereSL.cc localmpi.cc Particle.cc
haloprof_LDADD = $(EXPLIBS)

slabprof_SOURCES =  slabprof.cc PSP.cc header.cc
slabprof_LDADD = $(COMMON_LIBS) $(EXPLIBS)

psphisto_SOURCES = psphisto.cc PSP.cc header.cc localmpi.cc Particle.cc
psphisto_LDADD = $(EXPLIBS)

gas2d_SOURCES =  gas2d.cc PSP.cc header.cc localmpi.cc Particle.cc
gas2d_LDADD =    $(EXPLIBS)

gas2dcyl_SOURCES = gas2dcyl.cc PSP.cc header.cc localmpi.cc Particle.cc
gas2dcyl_LDADD = $(EXPLIBS)

makegrid_SOURCES = makegrid.cc header.cc PSP.cc
makegrid_LDADD = $(EXPLIBS) $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB) @exp_vtk_ldlib@ @exp_vtk_lib@ 

pspbox_SOURCES = pspbox.cc header.cc PSP.cc localmpi.cc
pspbox_LDADD = $(EXPLIBS) $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB)

eofpeek_SOURCES = eofpeek.cc

diskprof2_SOURCES = PSP.cc header.cc diskprof2.cc EmpOrth9thd.cc localmpi.cc Particle.cc
diskprof2_LDADD = $(EXPLIBS) $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB) @exp_vtk_ldlib@ @exp_vtk_lib@ 

pcatest_SOURCES = PSP.cc header.cc pcatest.cc EmpOrth9thd.cc localmpi.cc Particle.cc
pcatest_LDADD = $(EXPLIBS) $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB) @exp_vtk_ldlib@ @exp_vtk_lib@ 


#	SphericalSL.cc legendre.cc ddplgndr.cc sinecosine.cc 

EmpOrth9thd.cc:	$(top_srcdir)/src/EmpOrth9thd.cc
	cp $(top_srcdir)/src/$@ .

PSP.H:	$(top_srcdir)/utils/PhaseSpace/PSP.H
	cp $(top_srcdir)/utils/PhaseSpace/PSP.H .

PSP.cc:	$(top_srcdir)/utils/PhaseSpace/PSP.cc header.cc header.H PSP.H
	cp $(top_srcdir)/utils/PhaseSpace/PSP.cc .

header.cc:	$(top_srcdir)/utils/PhaseSpace/header.cc header.H
	cp $(top_srcdir)/utils/PhaseSpace/header.cc .

header.H:	$(top_srcdir)/utils/PhaseSpace/header.H
	cp $(top_srcdir)/utils/PhaseSpace/header.H .

localmpi.cc:	$(top_srcdir)/exputil/localmpi.cc
	cp $(top_srcdir)/exputil/$@ .

# Install and distribute these scripts

dist_bin_SCRIPTS = parse_outfrac ion_collide.py ion_energy.py	\
	ion_rate.py ion_frac.py ion_coll_number.py

CLEANFILES = 	EmpOrth9thd.cc localmpi.cc PSP.cc PSP.H header.cc header.H

## utils/ICs/Makefile.am ends here


