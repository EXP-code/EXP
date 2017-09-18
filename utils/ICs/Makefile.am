## Process this file with automake to produce Makefile.in

CC  = mpicc
CXX = mpiCC

if ENABLE_DMALLOC
ADD_CFLAGS3 = -DUSE_DMALLOC -g
ADD_CXXFLAGS3 = -DUSE_DMALLOC -g
ADD_LIBS3 = -ldmallocth -ldmallocthcxx
endif

AM_CFLAGS = -D_REENTRANT $(ADD_CFLAGS1) $(ADD_CFLAGS2) $(ADD_CFLAGS3)
AM_CXXFLAGS = -std=c++11 -D_REENTRANT -DSTANDALONE $(ADD_CXXFLAGS1) $(ADD_CXXFLAGS2) $(ADD_CXXFLAGS3)

AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/expgnu -I$(top_srcdir)/src -I$(top_srcdir)/utils/SL

COMMON_LIBS = -lexputil -lexpgnu $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)
MPI_LIBS = -lmpi

BOOST_LIBS = $(BOOST_SYSTEM_LIB) $(BOOST_PROGRAM_OPTIONS_LIB)

EXPLIBS = $(MPI_LIBS) $(COMMON_LIBS) $(ADD_LIBS3)
EXPFLAGS = -L$(top_srcdir)/exputil/.libs -L$(top_srcdir)/expgnu/.libs
CXXLINK = mpiCC $(CXXFLAGS) $(EXPFLAGS) -o $@

noinst_HEADERS =

bin_PROGRAMS = gendisk shrinkics gendisk2 gsphere pstmod gendisk3	\
	empinfo empdump makeIon phIon phConv

shrinkics_SOURCES = shrinkics.cc
shrinkics_LDADD = $(EXPLIBS)

gendisk_SOURCES = initial.cc DiskHalo.cc EmpOrth9thd.cc \
	SphericalSL.cc legendre.cc Particle.cc SParticle.cc \
	ddplgndr.cc sinecosine.cc AddDisk.cc localmpi.cc
gendisk_LDADD = $(EXPLIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

gendisk2_SOURCES =  initial2.cc DiskHalo2.cc EmpOrth9thd.cc \
	SphericalSL.cc legendre.cc Particle.cc SParticle.cc \
	ddplgndr.cc sinecosine.cc AddDisk.cc localmpi.cc
gendisk2_LDADD = $(EXPLIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

gendisk3_SOURCES =  initial3.cc DiskHalo3.cc EmpOrth9thd.cc \
	CylDisk.cc SphericalSL.cc legendre.cc Particle.cc SParticle.cc \
	ddplgndr.cc sinecosine.cc AddDisk.cc localmpi.cc
gendisk3_LDADD = $(EXPLIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

gsphere_SOURCES = bonnerebert.cc
gsphere_LDADD = $(EXPLIBS)

pstmod_SOURCES = pst_model.cc CylindricalDisk.cc EllipsoidForce.cc SLSphere.cc
pstmod_LDADD = $(EXPLIBS)

empinfo_SOURCES = empinfo.cc EmpOrth9thd.cc localmpi.cc
empinfo_LDADD = $(EXPLIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

empdump_SOURCES = empdump.cc EmpOrth9thd.cc localmpi.cc
empdump_LDADD = $(EXPLIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

makeIon_SOURCES = makeIonIC.cc Particle.cc Particle.H Species.H
makeIon_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS)
makeIon_LDADD = $(EXPLIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(BOOST_LIBS)	\
	$(LIBS) $(FLIBS)

phIon_SOURCES = ph_ion.cc
phIon_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS)
phIon_LDADD = $(BLAS_LIBS) $(BOOST_LIBS) $(LIBS) $(FLIBS)

phConv_SOURCES = ph_conv.cc
phConv_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS)
phConv_LDADD = $(BLAS_LIBS) $(BOOST_LIBS) $(LIBS) $(FLIBS)

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


CLEANFILES = EmpOrth9thd.cc localmpi.cc CylindricalDisk.cc	\
	SLSphere.cc Particle.cc Particle.H Particle.h Species.H

## utils/ICs/Makefile.am ends here


