## Makefile for the libsrnd directory of galaxy.
## Copyright (C) 2000 Aaron Culich
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

## Process this file with automake to produce Makefile.in

noinst_LIBRARIES = libsrnd.a
libsrnd_a_SOURCES = ACG.cc Binomial.cc DiscUnif.cc Erlang.cc Geom.cc \
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

## libsrnd/Makefile.am ends here
