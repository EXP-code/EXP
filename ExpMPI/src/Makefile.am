## Process this file with automake to produce Makefile.in

VPATH = $(top_srcdir)/src:.

CC = mpicc
CXX = mpiCC

CFLAGS = -O3
CXXFLAGS = -O3

if ENABLE_SEQCHECK
ADD_CFLAGS1 = -DSEQCHECK
ADD_CXXFLAGS1 = -DSEQCHECK
endif

if ENABLE_DEBUG
ADD_CFLAGS2 = -DDEBUG
ADD_CXXFLAGS2 = -DDEBUG
endif

if ENABLE_DMALLOC
ADD_CFLAGS3 = -DUSE_DMALLOC -g
ADD_CXXFLAGS3 = -DUSE_DMALLOC -g
ADD_LIBS3 = -ldmallocth -ldmallocthcxx
endif

AM_CFLAGS = -D_REENTRANT $(ADD_CFLAGS1) $(ADD_CFLAGS2) $(ADD_CFLAGS3)
AM_CXXFLAGS = -D_REENTRANT $(ADD_CXXFLAGS1) $(ADD_CXXFLAGS2) $(ADD_CXXFLAGS3)

CXXLINK = hf77 -rdynamic $(CXXFLAGS) $(EXPFLAGS) -o $@

exp_SOURCES = \
     Basis.cc Bessel.cc CBrock.cc Component.cc CBrockDisk.cc		     \
     Cube.cc Cylinder.cc EmpOrth9thd.cc ExternalCollection.cc                \
     ExternalForce.cc Hernquist.cc Orient.cc Particle.cc PotAccel.cc         \
     SLGridMP2.cc ScatterMFP.cc SphericalBasis.cc AxisymmetricBasis.cc       \
     Sphere.cc SphericalBasisMixtureSL.cc SphereTwoCenter.cc SphereEJCOM.cc  \
     global.cc begin.cc biorth1d.cc ddplgndr.cc Direct.cc end.cc             \
     OutputContainer.cc OutPS.cc OutPSN.cc OutCHKPT.cc Output.cc header.cc   \
     expand.cc externalShock.cc generateRelaxation.cc HaloBulge.cc           \
     incpos.cc incvel.cc initvel.cc ComponentContainer.cc                    \
     OutRelaxation.cc OrbTrace.cc OutDiag.cc OutLog.cc OutCoef.cc            \
     parse.cc Slab.cc SlabSL.cc step.cc syncvel.cc sbessj.cc sbessz.cc       \
     tidalField.cc ultra.cc ultrasphere.cc MPL.cc sledge.f zbrent.cc

INCLUDES = -I$(top_srcdir)/include -I$(top_srcdir)/expgnu -I$(top_srcdir)/src

bin_PROGRAMS = exp

COMMON_LIBS = -lexputil -lexpgnu -lstdc++
SYS_LIBS = -lpthread -ldl -lm
MPI_LIBS = -lmpi

LIBS = $(COMMON_LIBS) $(MPI_LIBS) $(SYS_LIBS)
LDADD = $(COMMON_LIBS) $(ADD_LIBS3) 
EXPFLAGS = -L$(top_srcdir)/exputil/.libs -L$(top_srcdir)/expgnu/.libs

noinst_HEADERS = Basis.H Bessel.H CBrock.H Component.H Cube.H Cylinder.H \
	ExternalCollection.H ExternalForce.H Hernquist.H Particle.H \
	Particle.H PotAccel.H ScatterMFP.H SphericalBasis.H externalShock.H \
	generateRelaxation.H tidalField.H EmpOrth9thd.h Orient.H OutPS.H \
	SLGridMP2.h biorth1d.h \
	coef.H exp_thread.h expand.h sltableMP2.h tipsydefs.h version.h

# special case for f77 for now

sledge.o:	sledge.f
		f77 -O3 -o sledge.o -c sledge.f

SUBDIRS = user

## src/Makefile.am ends here







