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

COMMON_LIBS = -lexputil -lexpgnu $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) @exp_vtk_ldlib@ @exp_vtk_lib@
MPI_LIBS = -lmpi

BOOST_LIBS = $(BOOST_SYSTEM_LIB) $(BOOST_FILESYSTEM_LIB) $(BOOST_PROGRAM_OPTIONS_LIB)

EXPLIBS = $(MPI_LIBS) $(COMMON_LIBS) $(ADD_LIBS3)
EXPFLAGS = -L$(top_srcdir)/exputil/.libs -L$(top_srcdir)/expgnu/.libs
CXXLINK = mpiCC $(CXXFLAGS) $(EXPFLAGS) -o $@

noinst_HEADERS =

if ENABLE_CUDA
CUDA_LIBFLAGS  = -lcudart
endif

bin_PROGRAMS = gendisk shrinkics gendisk2 gsphere pstmod gendisk3	\
	empinfo empdump makeIon phIon phConv

shrinkics_SOURCES = shrinkics.cc
shrinkics_LDADD = $(EXPLIBS)

gendisk_SOURCES = initial.cc DiskHalo.cc EmpOrth9thd.cc \
	SphericalSL.cc legendre.cc Particle.cc SParticle.cc \
	ddplgndr.cc sinecosine.cc AddDisk.cc localmpi.cc
gendisk_LDADD = $(COMMON_LIBS) $(EXPLIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS) $(CUDA_LIBFLAGS)

gendisk2_SOURCES =  initial2.cc DiskHalo2.cc EmpOrth9thd.cc \
	SphericalSL.cc legendre.cc Particle.cc SParticle.cc \
	ddplgndr.cc sinecosine.cc AddDisk.cc localmpi.cc
gendisk2_LDFLAGS = $(BOOST_LDFLAGS)
gendisk2_LDADD = $(COMMON_LIBS) $(EXPLIBS) $(BOOST_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS) $(CUDA_LIBFLAGS)

gendisk3_SOURCES =  initial3.cc DiskHalo3.cc EmpOrth9thd.cc \
	CylDisk.cc SphericalSL.cc legendre.cc Particle.cc SParticle.cc \
	ddplgndr.cc sinecosine.cc AddDisk.cc localmpi.cc
gendisk3_LDADD = $(COMMON_LIBS) $(EXPLIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS) $(CUDA_LIBFLAGS)

gsphere_SOURCES = bonnerebert.cc
gsphere_LDADD = $(COMMON_LIBS) $(EXPLIBS)

pstmod_SOURCES = pst_model.cc CylindricalDisk.cc EllipsoidForce.cc SLSphere.cc
pstmod_LDADD = $(COMMON_LIBS) $(EXPLIBS) $(CUDA_LIBFLAGS)

empinfo_SOURCES = empinfo.cc EmpOrth9thd.cc localmpi.cc
empinfo_LDADD = $(COMMON_LIBS) $(EXPLIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS) $(CUDA_LIBFLAGS)

empdump_SOURCES = empdump.cc EmpOrth9thd.cc localmpi.cc
empdump_LDADD = $(COMMON_LIBS) $(EXPLIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS) $(CUDA_LIBFLAGS)

makeIon_SOURCES = makeIonIC.cc Particle.cc Particle.H Species.H	\
	atomic_constants.H atomic_constants.cc
makeIon_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS)
makeIon_LDADD = $(COMMON_LIBS) $(EXPLIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(BOOST_LIBS)	\
	$(LIBS) $(FLIBS) $(CUDA_LIBFLAGS)

phIon_SOURCES = ph_ion.cc
phIon_LDFLAGS = $(COMMON_LIBS) $(EXPFLAGS) $(BOOST_LDFLAGS)
phIon_LDADD = $(BLAS_LIBS) $(BOOST_LIBS) $(LIBS) $(FLIBS) $(CUDA_LIBFLAGS)

phConv_SOURCES = ph_conv.cc
phConv_LDFLAGS = $(COMMON_LIBS) $(EXPFLAGS) $(BOOST_LDFLAGS)
phConv_LDADD = $(BLAS_LIBS) $(BOOST_LIBS) $(LIBS) $(CUDA_LIBFLAGS)

initial.cc:		Particle.H

initial2.cc:		Particle.H

initial3.cc:		Particle.H

empinfo.cc:		Particle.H

empdump.cc:		Particle.H

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


