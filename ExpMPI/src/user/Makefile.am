CC= mpicc
CXX = mpiCC

INCLUDES = -I$(top_srcdir)/include -I$(top_srcdir)/expgnu -I$(top_srcdir)/src

CFLAGS = -O3 -fPIC -Wall
CXXFLAGS = -O3 -fPIC -Wall
LDFLAGS = -L$(top_srcdir)/expgnu -L$(top_srcdir)/exputil

libdir = $(prefix)/lib/user

lib_LTLIBRARIES = libuser.la libbar.la libsat.la

libuser_la_SOURCES = UserTest.cc

libbar_la_SOURCES = UserBar.cc

libsat_la_SOURCES = UserSat.cc SatelliteOrbit.cc ParamParse.cc
libsat_la_LDFLAGS = -L$(top_srcdir)/exputil/.libs \
	-L$(top_srcdir)/expgnu/.libs -lstdc++
