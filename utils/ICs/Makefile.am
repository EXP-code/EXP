## Process this file with automake to produce Makefile.in

CC  = mpicc
CXX = mpiCC

if ENABLE_DMALLOC
ADD_CFLAGS3 = -DUSE_DMALLOC -g
ADD_CXXFLAGS3 = -DUSE_DMALLOC -g
ADD_LIBS3 = -ldmallocth -ldmallocthcxx
endif

AM_CFLAGS = -D_REENTRANT $(ADD_CFLAGS1) $(ADD_CFLAGS2) $(ADD_CFLAGS3)
AM_CXXFLAGS = -fPIC -std=c++11 -D_REENTRANT -DSTANDALONE $(ADD_CXXFLAGS1) $(ADD_CXXFLAGS2) $(ADD_CXXFLAGS3) @exp_vtk_include@ $(BOOST_CPPFLAGS)

AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/expgnu -I$(top_srcdir)/src -I$(top_srcdir)/utils/SL

if ENABLE_CUDA
CUDA_LIBFLAGS  = -lcudart
endif

COMMON_LIBS = -lexputil -lexpgnu $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) @exp_vtk_ldlib@ @exp_vtk_lib@
MPI_LIBS = -lmpi

BOOST_LIBS = $(BOOST_SYSTEM_LIB) $(BOOST_FILESYSTEM_LIB) $(BOOST_PROGRAM_OPTIONS_LIB)

EXPLIBS = $(MPI_LIBS) $(COMMON_LIBS) $(ADD_LIBS3) $(CUDA_LIBFLAGS)
EXPFLAGS = -L$(top_srcdir)/exputil/.libs -L$(top_srcdir)/expgnu/.libs $(CUDA_LDFLAGS)
CXXLINK = mpiCC $(CXXFLAGS) $(EXPFLAGS) -o $@

noinst_HEADERS =

bin_PROGRAMS = gendisk shrinkics gendisk2 gsphere pstmod gendisk3	\
	empinfo empdump makeIon phIon phConv

shrinkics_SOURCES = shrinkics.cc
shrinkics_LDADD = $(EXPLIBS)

gendisk_SOURCES = initial.cc DiskHalo.cc EmpOrth9thd.cc \
	SphericalSL.cc legendre.cc Particle.cc SParticle.cc \
	ddplgndr.cc sinecosine.cc AddDisk.cc localmpi.cc
gendisk_LDADD = $(COMMON_LIBS) $(EXPLIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

gendisk2_SOURCES =  initial2.cc DiskHalo2.cc EmpOrth9thd.cc \
	SphericalSL.cc legendre.cc Particle.cc SParticle.cc \
	ddplgndr.cc sinecosine.cc AddDisk.cc localmpi.cc
gendisk2_LDFLAGS =  $(EXPFLAGS) $(BOOST_LDFLAGS)
gendisk2_LDADD = $(COMMON_LIBS) $(EXPLIBS) $(BOOST_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

gendisk3_SOURCES =  initial3.cc DiskHalo3.cc EmpOrth9thd.cc \
	CylDisk.cc SphericalSL.cc legendre.cc Particle.cc SParticle.cc \
	ddplgndr.cc sinecosine.cc AddDisk.cc localmpi.cc
gendisk3_LDFLAGS = $(EXPFLAGS)
gendisk3_LDADD = $(COMMON_LIBS) $(EXPLIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

gsphere_SOURCES = bonnerebert.cc
gsphere_LDADD = $(COMMON_LIBS) $(EXPLIBS)

pstmod_SOURCES = pst_model.cc CylindricalDisk.cc EllipsoidForce.cc SLSphere.cc
pstmod_LDFLAGS = $(EXPFLAGS)
pstmod_LDADD = $(COMMON_LIBS) $(EXPLIBS)

empinfo_SOURCES = empinfo.cc EmpOrth9thd.cc localmpi.cc
empinfo_LDFLAGS = $(EXPFLAGS)
empinfo_LDADD = $(COMMON_LIBS) $(EXPLIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

empdump_SOURCES = empdump.cc EmpOrth9thd.cc localmpi.cc
empdump_LDFLAGS = $(EXPFLAGS)
empdump_LDADD = $(COMMON_LIBS) $(EXPLIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

makeIon_SOURCES = makeIonIC.cc Particle.cc Particle.H Species.H \
	atomic_constants.H atomic_constants.cc
makeIon_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS)
makeIon_LDADD = $(COMMON_LIBS) $(EXPLIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(BOOST_LIBS)	\
	$(LIBS) $(FLIBS)

phIon_SOURCES = ph_ion.cc
phIon_LDFLAGS = $(COMMON_LIBS) $(EXPFLAGS) $(BOOST_LDFLAGS)
phIon_LDADD = $(BLAS_LIBS) $(BOOST_LIBS) $(LIBS) $(FLIBS)

phConv_SOURCES = ph_conv.cc
phConv_LDFLAGS = $(COMMON_LIBS) $(EXPFLAGS) $(BOOST_LDFLAGS)
phConv_LDADD = $(BLAS_LIBS) $(BOOST_LIBS) $(LIBS)

initial.cc:		Particle.H

initial2.cc:		Particle.H

initial3.cc:		Particle.H

empinfo.cc:		Particle.H

empdump.cc:		Particle.H

makeIonIC.cc:		atomic_constants.H

CylindricalDisk.cc:	$(top_srcdir)/utils/SL/CylindricalDisk.cc
			cp $(top_srcdir)/utils/SL/$@ .

SLSphere.cc:		$(top_srcdir)/utils/SL/SLSphere.cc
			cp $(top_srcdir)/utils/SL/$@ .

EmpOrth9thd.cc:		$(top_srcdir)/src/EmpOrth9thd.cc Particle.H
			cp $(top_srcdir)/src/$@ .

localmpi.cc:		$(top_srcdir)/exputil/localmpi.cc
			cp $(top_srcdir)/exputil/$@ .

Particle.cc:		$(top_srcdir)/src/Particle.cc Particle.H
			cp $(top_srcdir)/src/$@ .

Particle.H:		$(top_srcdir)/src/Particle.H
			cp $(top_srcdir)/src/$@ .
			ln -s Particle.H Particle.h

Species.H:		$(top_srcdir)/src/Species.H
			cp $(top_srcdir)/src/$@ .

makeIonIC.cc:		atomic_constants.H

atomic_constants.cc:	$(top_srcdir)/src/user/atomic_constants.cc atomic_constants.H
			cp $(top_srcdir)/src/user/$@ .

atomic_constants.H:	$(top_srcdir)/src/user/atomic_constants.H
			cp $(top_srcdir)/src/user/$@ .


CLEANFILES = EmpOrth9thd.cc localmpi.cc CylindricalDisk.cc	\
	SLSphere.cc Particle.cc Particle.H Particle.h Species.H \
	atomic_constants.H atomic_constants.cc

## utils/ICs/Makefile.am ends here


