CCLD  = --tag=F77  mpif77
CXXLD = --tag=F77 mpif77

INCLUDES = -I$(top_srcdir)/include -I$(top_srcdir)/expgnu -I$(top_srcdir)/src

AM_CFLAGS = -fPIC -Wall -D_REENTRANT
AM_CXXFLAGS = -fPIC -Wall -D_REENTRANT
AM_LDFLAGS = -L$(top_srcdir)/expgnu -L$(top_srcdir)/exputil

libdir = $(prefix)/lib/user

bin_PROGRAMS = testcoord

lib_LTLIBRARIES = libuser.la libbar.la libebar.la libsat.la libdiffuse.la \
	libsatfix.la libdiffrot.la liblogpot.la libdisk.la libtorque.la \
	libebarp.la librespot.la librespotn.la librptest.la libhalo.la

libuser_la_SOURCES = UserTest.cc

libbar_la_SOURCES = UserBar.cc

libebar_la_SOURCES = UserEBar.cc EllipForce.cc

libebarp_la_SOURCES = UserEBarP.cc EllipForce.cc

libsat_la_SOURCES = UserSat.cc SatelliteOrbit.cc
#libsat_la_LDFLAGS = -L$(top_srcdir)/exputil/.libs \
#	-L$(top_srcdir)/expgnu/.libs -lstdc++

libdiffuse_la_SOURCES = UserDiffuse.cc

libsatfix_la_SOURCES = SatFix.cc

libdiffrot_la_SOURCES = UserDiffRot.cc

liblogpot_la_SOURCES = UserLogPot.cc

libdisk_la_SOURCES = UserDisk.cc

libhalo_la_SOURCES = UserHalo.cc

libtorque_la_SOURCES = UserTorque.cc

librespot_la_SOURCES = UserResPot.cc ResPot.cc BarForcing.cc Perturbation.cc
librespot_la_CFLAGS = $(AM_CFLAGS) 	#Force removal of .o

librespotn_la_SOURCES = UserResPotN.cc ResPot.cc BarForcing.cc Perturbation.cc
librespotn_la_CFLAGS = $(AM_CFLAGS)	#Force removal of .o

librptest_la_SOURCES = UserRPtest.cc ResPot.cc
librptest_la_CFLAGS = $(AM_CFLAGS)	#Force removal of .o

testcoord_SOURCES = testcoord.cc ResPot.cc BarForcing.cc Perturbation.cc getkey.c
testcoord_LDFLAGS = -L$(top_srcdir)/exputil -L$(top_srcdir)/expgnu -lexputil -lexpgnu -lstdc++ -llammpi++ -lmpi
testcoord_CFLAGS = $(AM_CFLAGS)	#Force removal of .o
