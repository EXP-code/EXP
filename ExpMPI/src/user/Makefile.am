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

#
# User utility library for duplicated source components
#
COMMONSRC = ResPot.cc BarForcing.cc Perturbation.cc \
	EllipForce.cc SatelliteOrbit.cc

COMMONOBJ = $(COMMONSRC:.cc=.o)


libutils.a:	$(COMMONOBJ)
	ar rcs libutils.a $(COMMONOBJ)


#
# Shared (loadable) libraries
#

libuser_la_SOURCES = UserTest.cc

libbar_la_SOURCES = UserBar.cc

libebar_la_SOURCES = UserEBar.cc
libebar_la_DEPENDS = libutils.a
libebar_la_LIBADD = -lutils

libebarp_la_SOURCES = UserEBarP
libebarp_la_DEPENDS = libutils.a
libebarp_la_LIBADD = -lutils

libsat_la_SOURCES = UserSat.cc
libsat_la_DEPENDS = libutils.a
libsat_la_LIBADD = -lutils

libdiffuse_la_SOURCES = UserDiffuse.cc

libsatfix_la_SOURCES = SatFix.cc

libdiffrot_la_SOURCES = UserDiffRot.cc

liblogpot_la_SOURCES = UserLogPot.cc

libdisk_la_SOURCES = UserDisk.cc

libhalo_la_SOURCES = UserHalo.cc

libtorque_la_SOURCES = UserTorque.cc

librespot_la_SOURCES = UserResPot.cc 
librespot_la_CFLAGS = $(AM_CFLAGS) 	#Force removal of .o
librespot_la_DEPENDS = libutils.a
librespot_la_LIBADD = -lutils

librespotn_la_SOURCES = UserResPotN.cc
librespotn_la_CFLAGS = $(AM_CFLAGS)	#Force removal of .o
librespotn_la_DEPENDS = libutils.a
librespotn_la_LIBADD = -lutils

librptest_la_SOURCES = UserRPtest.cc
librptest_la_CFLAGS = $(AM_CFLAGS)	#Force removal of .o
librptest_la_DEPENDS = libutils.a
librptest_la_LIBADD = -lutils

testcoord_SOURCES = testcoord.cc getkey.c
testcoord_LDFLAGS = -L$(top_srcdir)/exputil -L$(top_srcdir)/expgnu -lutil -lexputil -lexpgnu -lstdc++ -llammpi++ -lmpi
testcoord_CFLAGS = $(AM_CFLAGS)	#Force removal of .o
