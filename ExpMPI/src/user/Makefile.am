CCLD  = --tag=F77  mpif77
CXXLD = --tag=F77 mpif77

INCLUDES = -I$(top_srcdir)/include -I$(top_srcdir)/expgnu -I$(top_srcdir)/src

AM_CFLAGS = -fPIC -Wall -D_REENTRANT
AM_CXXFLAGS = -fPIC -Wall -D_REENTRANT -Wno-non-virtual-dtor
AM_LDFLAGS = -L$(top_srcdir)/expgnu -L$(top_srcdir)/exputil

libdir = $(prefix)/lib/user

# bin_PROGRAMS = testcoord

lib_LTLIBRARIES = libuser.la libbar.la libebar.la libsat.la libdiffuse.la \
	libsatfix.la libdiffrot.la liblogpot.la libdisk.la libtorque.la \
	libebarp.la librespot.la librespotn.la librespotorb.la \
	librptest.la libhalo.la libebars.la libwake.la libsatfixorb.la \
	libdsmc.la libebarn.la libtreegas.la

#
# Convenience library for duplicated source components (static)
#
noinst_LIBRARIES = libcommon.a libdsmc.a

noinst_HEADERS = BarForcing.H Perturbation.H UserDiffuse.H UserResPot.H \
	CircularOrbit.H ResPot.H UserDisk.H UserResPotN.H Collide.H \
	ResPotOrb.H UserEBar.H UserResPotOrb.H SatFix.H UserEBarN.H \
	UserTest.H SatFixOrb.H UserEBarP.H UserTreeDSMC.H \
	CollideLTE.H TwoBodyDiffuse.H UserEBarS.H UserWake.H EllipForce.H \
	UserBar.H UserHalo.H HeatCool.H UserDSMC.H UserLogPot.H \
	UserDiffRot.H UserRPtest.H

libcommon_a_SOURCES = ResPot.cc ResPotOrb.cc BarForcing.cc Perturbation.cc \
	EllipForce.cc SatelliteOrbit.cc CircularOrbit.cc TwoBodyDiffuse.cc

libdsmc_a_SOURCES = Collide.cc CollideLTE.cc HeatCool.cc

#
# Shared (loadable) libraries
#
libuser_la_SOURCES = UserTest.cc

libbar_la_SOURCES = UserBar.cc

libebar_la_SOURCES = UserEBar.cc
libebar_la_LIBADD = libcommon.a

libebarn_la_SOURCES = UserEBarN.cc
libebarn_la_LIBADD = libcommon.a

libebarp_la_SOURCES = UserEBarP.cc
libebarp_la_LIBADD = libcommon.a

libebars_la_SOURCES = UserEBarS.cc
libebars_la_LIBADD = libcommon.a

libsat_la_SOURCES = UserSat.cc
libsat_la_LIBADD = libcommon.a

libdiffuse_la_SOURCES = UserDiffuse.cc

libsatfix_la_SOURCES = SatFix.cc

libsatfixorb_la_SOURCES = SatFixOrb.cc
libsatfixorb_la_LIBADD = libcommon.a

libdiffrot_la_SOURCES = UserDiffRot.cc

liblogpot_la_SOURCES = UserLogPot.cc

libdisk_la_SOURCES = UserDisk.cc

libhalo_la_SOURCES = UserHalo.cc

libtorque_la_SOURCES = UserTorque.cc

librespot_la_SOURCES = UserResPot.cc 
librespot_la_LIBADD = libcommon.a

librespotn_la_SOURCES = UserResPotN.cc
librespotn_la_LIBADD = libcommon.a

librespotorb_la_SOURCES = UserResPotOrb.cc
librespotorb_la_LIBADD = libcommon.a

librptest_la_SOURCES = UserRPtest.cc
librptest_la_LIBADD = libcommon.a

libwake_la_SOURCES = UserWake.cc
libwake_la_LIBADD = libcommon.a

libdsmc_la_SOURCES = UserDSMC.cc
libdsmc_la_LIBADD = libcommon.a

libtreegas_la_SOURCES = UserTreeDSMC.cc 
libtreegas_la_LIBADD = libcommon.a libdsmc.a

#
# Test routine
#
# testcoord_SOURCES = testcoord.cc getkey.c
# testcoord_LDADD = libcommon.a
# testcoord_LDFLAGS = -L$(top_srcdir)/exputil -L$(top_srcdir)/expgnu -lexputil -lexpgnu -lstdc++ -lmpi_cxx -lmpi
# testcoord_LDFLAGS = -L$(top_srcdir)/exputil -L$(top_srcdir)/expgnu -lexputil -lexpgnu -L/usr/lib/gcc/i486-linux-gnu/4.0.1 -lstdc++ -lmpi_cxx -lmpi
