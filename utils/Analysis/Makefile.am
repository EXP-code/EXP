## Process this file with automake to produce Makefile.in

CC  = mpicc
CXX = mpiCC

if ENABLE_DMALLOC
ADD_CFLAGS3   = -DUSE_DMALLOC -g
ADD_CXXFLAGS3 = -DUSE_DMALLOC -g
ADD_LIBS3     = -ldmallocth -ldmallocthcxx
endif

AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/expgnu -I$(top_srcdir)/src $(CUDA_CFLAGS)
AM_CFLAGS = -D_REENTRANT $(ADD_CFLAGS3)
AM_CXXFLAGS = -std=c++11 -D_REENTRANT -DSTANDALONE $(INCLUDES) $(ADD_CXXFLAGS3) $(BOOST_CPPFLAGS) @exp_vtk_include@

if ENABLE_CUDA
CUDA_LIBFLAGS  = -lcudart
endif

COMMON_LIBS = $(LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) @exp_vtk_ldlib@ @exp_vtk_lib@ 
MPI_LIBS = -lmpi
EXPLIBS = -lexputil -lexpgnu $(COMMON_LIBS) $(ADD_LIBS3) $(MPI_LIBS) $(CUDA_LIBFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_SYSTEM_LIB) $(YAML)
EXPFLAGS = -L$(top_srcdir)/exputil/.libs -L$(top_srcdir)/expgnu/.libs $(CUDA_LDFLAGS) @exp_omp_lib@
CXXLINK = mpiCC $(CXXFLAGS) $(EXPFLAGS) -o $@

noinst_HEADERS =

bin_PROGRAMS = diskprof haloprof haloprof2 slabprof psphisto gas2d    \
gas2dcyl makegrid pspbox eofpeek diskprof2 pcatest diskprof3 mssaprofD \
mssaprofH readcoefs readcoefs2 pcatest diskprof4 sphprof yaml_diff

diskprof_SOURCES = PSP.cc header.cc diskprof.cc EmpCylSL.cc localmpi.cc Particle.cc
diskprof_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS) 
diskprof_LDADD   = $(EXPLIBS) @exp_vtk_ldlib@ @exp_vtk_lib@ 

haloprof_SOURCES = PSP.cc header.cc haloprof.cc SphereSL.cc localmpi.cc Particle.cc
haloprof_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS)
haloprof_LDADD = $(EXPLIBS)

haloprof2_SOURCES = PSP2.cc header.cc haloprof2.cc SphereSL.cc localmpi.cc Particle.cc
haloprof2_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS)
haloprof2_LDADD = $(EXPLIBS)

slabprof_SOURCES =  slabprof.cc PSP.cc header.cc
slabprof_LDFLAGS =  $(EXPFLAGS) $(BOOST_LDFLAGS)
slabprof_LDADD   = $(COMMON_LIBS) $(EXPLIBS) $(CUDA_LIBFLAGS)

psphisto_SOURCES = psphisto.cc PSP.cc header.cc localmpi.cc Particle.cc
psphisto_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS)
psphisto_LDADD   = $(EXPLIBS)

gas2d_SOURCES =  gas2d.cc PSP.cc header.cc localmpi.cc Particle.cc
gas2d_LDFLAGS =  $(EXPFLAGS) $(BOOST_LDFLAGS)
gas2d_LDADD   =  $(EXPLIBS)

gas2dcyl_SOURCES = gas2dcyl.cc PSP.cc header.cc localmpi.cc Particle.cc
gas2dcyl_LDFLAGS = $(EXPFLAGS)  $(BOOST_LDFLAGS)
gas2dcyl_LDADD   = $(EXPLIBS)

makegrid_SOURCES = makegrid.cc header.cc PSP.cc
makegrid_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS)
makegrid_LDADD   = $(EXPLIBS) $(BOOST_PROGRAM_OPTIONS_LIB) @exp_vtk_ldlib@ @exp_vtk_lib@

pspbox_SOURCES = pspbox.cc header.cc PSP.cc localmpi.cc
pspbox_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS) 
pspbox_LDADD   = $(EXPLIBS) $(BOOST_PROGRAM_OPTIONS_LIB)

eofpeek_SOURCES = eofpeek.cc

diskprof2_SOURCES = PSP.cc header.cc diskprof2.cc EmpCylSL.cc localmpi.cc Particle.cc
diskprof2_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS) 
diskprof2_LDADD   = $(EXPLIBS) $(BOOST_PROGRAM_OPTIONS_LIB) @exp_vtk_ldlib@ @exp_vtk_lib@

mssaprofD_SOURCES  = PSP.cc header.cc mssaprof_disk.cc EmpCylSL.cc localmpi.cc Particle.cc
mssaprofD_LDFLAGS  = $(EXPFLAGS) $(BOOST_LDFLAGS) 
mssaprofD_LDADD    = $(EXPLIBS) $(BOOST_PROGRAM_OPTIONS_LIB) @exp_vtk_ldlib@ @exp_vtk_lib@

mssaprofH_SOURCES  = PSP.cc header.cc mssaprof_halo.cc SphereSL.cc localmpi.cc Particle.cc
mssaprofH_LDFLAGS  = $(EXPFLAGS) $(BOOST_LDFLAGS) 
mssaprofH_LDADD    = $(EXPLIBS) $(BOOST_PROGRAM_OPTIONS_LIB) @exp_vtk_ldlib@ @exp_vtk_lib@

pcatest_SOURCES = PSP.cc header.cc pcatest.cc EmpCylSL.cc localmpi.cc Particle.cc
pcatest_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS) 
pcatest_LDADD   = $(EXPLIBS) $(BOOST_PROGRAM_OPTIONS_LIB) @exp_vtk_ldlib@ @exp_vtk_lib@

diskprof3_SOURCES = diskprof3.cc EmpCylSL.cc localmpi.cc
diskprof3_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS) 
diskprof3_LDADD   = $(EXPLIBS) $(BOOST_PROGRAM_OPTIONS_LIB) @exp_vtk_ldlib@ @exp_vtk_lib@

readcoefs_SOURCES = readcoefs.cc
readcoefs_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS) 
readcoefs_LDADD   = $(BOOST_PROGRAM_OPTIONS_LIB)

readcoefs2_SOURCES = readcoefs2.cc
readcoefs2_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS) 
readcoefs2_LDADD   = $(BOOST_PROGRAM_OPTIONS_LIB)

diskprof4_SOURCES = PSP2.cc header.cc diskprof4.cc EmpCylSL.cc localmpi.cc Particle.cc
diskprof4_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS) 
diskprof4_LDADD   = $(EXPLIBS) $(BOOST_PROGRAM_OPTIONS_LIB) @exp_vtk_ldlib@ @exp_vtk_lib@

sphprof_SOURCES = PSP2.cc header.cc sphprof.cc SphereSL.cc	\
			localmpi.cc Particle.cc
sphprof_LDFLAGS   = $(EXPFLAGS) $(BOOST_LDFLAGS)
sphprof_LDADD     = $(EXPLIBS)

yaml_diff_SOURCES = yaml_diff.cc
yaml_diff_LDFLAGS   = $(BOOST_LDFLAGS)
yaml_diff_LDADD   = $(BOOST_PROGRAM_OPTIONS_LIB) $(YAML)

EmpCylSL.cc:	$(top_srcdir)/src/EmpCylSL.cc
	cp $(top_srcdir)/src/$@ .

PSP.H:	$(top_srcdir)/utils/PhaseSpace/PSP.H
	cp $(top_srcdir)/utils/PhaseSpace/PSP.H .

PSP.cc:	$(top_srcdir)/utils/PhaseSpace/PSP.cc header.cc header.H PSP.H
	cp $(top_srcdir)/utils/PhaseSpace/PSP.cc .

PSP2.H:	$(top_srcdir)/utils/PhaseSpace/PSP2.H
	cp $(top_srcdir)/utils/PhaseSpace/PSP2.H .

PSP2.cc:	$(top_srcdir)/utils/PhaseSpace/PSP2.cc header.cc header.H PSP2.H #
	cp $(top_srcdir)/utils/PhaseSpace/PSP2.cc .

header.cc:	$(top_srcdir)/utils/PhaseSpace/header.cc header.H
	cp $(top_srcdir)/utils/PhaseSpace/header.cc .

header.H:	$(top_srcdir)/utils/PhaseSpace/header.H
	cp $(top_srcdir)/utils/PhaseSpace/header.H .

localmpi.cc:	$(top_srcdir)/exputil/localmpi.cc
	cp $(top_srcdir)/exputil/$@ .

# Install and distribute these scripts

dist_bin_SCRIPTS = parse_outfrac ion_collide.py ion_energy.py	\
	ion_rate.py ion_frac.py ion_coll_number.py

CLEANFILES = 	EmpCylSL.cc localmpi.cc PSP.cc PSP.H header.cc header.H

## utils/ICs/Makefile.am ends here


