noinst_LIBRARIES = libmdw.a

libmdw_a_SOURCES = brent.c cel.c cutil.c dgammln.c dgammp.c \
	dgammq.c dgcf.c dgser.c dplgndr.c el2.c exponi.c \
	f1dim.c fpeieee.c frprmn.c gasdev.c gauss.c \
	getkey.c gknots.c indexx.c lagknots.c linmin.c \
	locate.c lubksb.c ludcmp.c mnbrak.c odd2.c plgndr.c \
	poidev.c rombe2d.c sbessj.c sbessz.c sort.c sort2.c \
	spline.c splint.c splint1.c splint1d.c splsum.c tensor.c \
	xrand.c zbrent.c

noinst_HEADERS = cutil.h

# AM_CFLAGS = -Wall -fPIC
# AM_CXXFLAGS = -Wall -fPIC
AM_CFLAGS = -fPIC -D_REENTRANT
AM_CXXFLAGS = -fPIC -D_REENTRANT
# CFLAGS = -g -DDEBUG
# CXXFLAGS = -g -DDEBUG
CFLAGS = -O3
CXXFLAGS = -O3
