CCLD  = --tag=F77  mpif77
CXXLD = --tag=F77 mpif77

INCLUDES = -I$(top_srcdir)/include -I$(top_srcdir)/expgnu -I$(top_srcdir)/src

AM_CFLAGS = -fPIC -Wall -D_REENTRANT
AM_CXXFLAGS = -fPIC -Wall -D_REENTRANT
AM_LDFLAGS = -L$(top_srcdir)/expgnu -L$(top_srcdir)/exputil

libdir = $(prefix)/lib/user

lib_LTLIBRARIES = libuser.la libbar.la libebar.la libsat.la libdiffuse.la \
	libsatfix.la libdiffrot.la liblogpot.la libdisk.la

libuser_la_SOURCES = UserTest.cc

libbar_la_SOURCES = UserBar.cc

libebar_la_SOURCES = UserEBar.cc EllipForce.cc

libsat_la_SOURCES = UserSat.cc SatelliteOrbit.cc ParamParse.cc
libsat_la_LDFLAGS = -L$(top_srcdir)/exputil/.libs \
	-L$(top_srcdir)/expgnu/.libs -lstdc++

libdiffuse_la_SOURCES = UserDiffuse.cc

libsatfix_la_SOURCES = SatFix.cc

libdiffrot_la_SOURCES = UserDiffRot.cc

liblogpot_la_SOURCES = UserLogPot.cc

libdisk_la_SOURCES = UserDisk.cc

