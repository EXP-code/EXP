## Process this file with automake to produce Makefile.in

LIBS = -lpthread -lmpi -lm

lib_LTLIBRARIES = libexpgnu.la

libexpgnu_la_SOURCES = ACG.cc Binomial.cc DiscUnif.cc Erlang.cc Geom.cc \
	HypGeom.cc LogNorm.cc MLCG.cc NegExp.cc Normal.cc Poisson.cc \
	RNG.cc Random.cc Regex.cc RndInt.cc Uniform.cc Weibull.cc \
	Gamma.cc String.cc SmplHist.cc SmplStat.cc \
	error.cc timer.c

noinst_HEADERS = ACG.h Binomial.h DiscUnif.h Erlang.h Gamma.h Geom.h \
	HypGeom.h LogNorm.h MLCG.h NegExp.h Normal.h Poisson.h RNG.h \
	Random.h Regex.h RndInt.h SmplHist.h SmplStat.h String.h \
	Uniform.h Weibull.h builtin.h rx.h std.h


# AM_CFLAGS = -Wall -fPIC
# AM_CXXFLAGS = -Wall -fPIC
AM_CFLAGS = -fPIC -D_REENTRANT
AM_CXXFLAGS = -fPIC -D_REENTRANT
# CFLAGS = -g -DDEBUG
# CXXFLAGS = -g -DDEBUG
CFLAGS = -O3
CXXFLAGS = -O3

## libexpgnu/Makefile.am ends here
