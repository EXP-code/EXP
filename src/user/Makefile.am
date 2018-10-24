AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/expgnu	\
-I$(top_srcdir)/src -I$(top_srcdir)/src/DSMC -I/usr/include/mpi	\
$(CUDA_CFLAGS)

BUILD_FLAGS =  -L$(top_srcdir)/expgnu/.libs -L$(top_srcdir)/exputil/.libs

if ENABLE_GPTL
  ADD_CFLAGS4  = -DUSE_GPTL $(gptl_include)
endif

AM_CFLAGS = -fPIC -Wall -D_REENTRANT $(ADD_CFLAGS4) -fopenmp
AM_FFLAGS = -fPIC -Wall -D_REENTRANT $(ADD_CFLAGS4) -fopenmp
AM_CXXFLAGS = -std=c++11 -fPIC -fno-inline-small-functions -Wall	\
-D_REENTRANT -Wno-non-virtual-dtor $(ADD_CFLAGS4) -fopenmp
AM_LDFLAGS = -L$(top_srcdir)/expgnu -L$(top_srcdir)/exputil -lgomp

libdir = $(prefix)/lib/user

lib_LTLIBRARIES = libuser.la libbar.la libebar.la libsat.la	\
	libdiffuse.la libsatfix.la libdiffrot.la liblogpot.la	\
	libdisk.la libtorque.la libebarp.la librespot.la	\
	librespotn.la librespotorb.la librptest.la libhalo.la	\
	libebars.la libwake.la libsatfixorb.la libebarn.la	\
	libperiodic.la libatmos.la libreflect.la libsatwake.la	\
	libpst.la libuserslab.la libtrigger.la libshear.la	\
	libsnheat.la libprofile.la

#
# Convenience library for duplicated source components (static)
#
noinst_LIBRARIES = libcommon.a libswake.a
noinst_HEADERS = BarForcing.H Perturbation.H UserDiffuse.H		\
	UserResPot.H CircularOrbit.H ResPot.H UserDisk.H		\
	UserResPotN.H ResPotOrb.H UserEBar.H UserResPotOrb.H SatFix.H	\
	UserEBarN.H UserTest.H SatFixOrb.H UserEBarP.H			\
	TwoBodyDiffuse.H UserEBarS.H UserWake.H EllipForce.H		\
	UserBar.H UserHalo.H UserLogPot.H UserDiffRot.H UserRPtest.H	\
	UserAtmos.H UserReflect.H UserPST.H UserShear.H UserSNheat.H	\
	UserProfile.H

libcommon_a_SOURCES = ResPot.cc ResPotOrb.cc BarForcing.cc	\
	Perturbation.cc EllipForce.cc SatelliteOrbit.cc		\
	CircularOrbit.cc TwoBodyDiffuse.cc TimeSeriesCoefs.cc	\
	LinearOrbit.cc UnboundOrbit.cc

libswake_a_SOURCES = respmat6.cc cratint.cc rombeveccmp.cc

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

libperiodic_la_SOURCES = UserPeriodic.cc 
libperiodic_la_LIBADD  = libcommon.a

libsnheat_la_SOURCES = UserSNheat.cc 
libsnheat_la_LIBADD  = libcommon.a

libatmos_la_SOURCES = UserAtmos.cc 
libatmos_la_LIBADD  = libcommon.a

libreflect_la_SOURCES = UserReflect.cc 

libsatwake_la_SOURCES = UserSatWake.cc invsvd.cc
libsatwake_la_LIBADD  = libswake.a libcommon.a $(LAPACK_LIBS) $(BLAS_LIBS)

libpst_la_SOURCES  = UserPST.cc CylindricalDisk.cc EllipsoidForce.cc SLSphere.cc
libpst_la_CXXFLAGS = $(AM_CXXFLAGS) -I$(top_srcdir)/utils/SL -I$(top_srcdir)/utils/ICs
libpst_la_LIBADD = libcommon.a

libuserslab_la_SOURCES = UserSlabHalo.cc 

libtrigger_la_SOURCES = UserBarTrigger.cc

libshear_la_SOURCES = UserShear.cc

libprofile_la_SOURCES = UserProfile.cc


CylindricalDisk.cc:	$(top_srcdir)/utils/SL/CylindricalDisk.cc
			cp $(top_srcdir)/utils/SL/$@ .

SLSphere.cc:		$(top_srcdir)/utils/SL/SLSphere.cc
			cp $(top_srcdir)/utils/SL/$@ .

EllipsoidForce.cc:	$(top_srcdir)/utils/ICs/EllipsoidForce.cc
			cp $(top_srcdir)/utils/ICs/$@ .

CLEANFILES = 	CylindricalDisk.cc EllipsoidForce.cc SLSphere.cc
