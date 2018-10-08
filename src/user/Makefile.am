AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/expgnu	\
-I$(top_srcdir)/src -I/usr/include/mpi
BUILD_FLAGS =  -L$(top_srcdir)/expgnu/.libs -L$(top_srcdir)/exputil/.libs

if ENABLE_GPTL
  ADD_CFLAGS4  = -DUSE_GPTL $(gptl_include)
endif

AM_CFLAGS   = -fPIC -Wall -D_REENTRANT $(ADD_CFLAGS4) -fopenmp
AM_FFLAGS   = -fPIC -Wall -D_REENTRANT $(ADD_CFLAGS4) -fopenmp
AM_CXXFLAGS = -std=c++11 -fPIC -fno-inline-small-functions -Wall \
	-D_REENTRANT -Wno-non-virtual-dtor $(ADD_CFLAGS4) -fopenmp
AM_LDFLAGS = -L$(top_srcdir)/expgnu -L$(top_srcdir)/exputil -lgomp

if ENABLE_CUDA
  ADD_CUDA = cudaIon.cu cudaCollideIon.cu cudaElastic.cu
  OBJ_CUDA = cudaIon.o cudaCollideIon.o cudaElastic.o
  ifeq ($(GENCODE_FLAGS),)
    $(foreach sm,@sm_vers@,$(eval GENCODE_FLAGS += -gencode arch=compute_$(sm),code=sm_$(sm)))
  endif
  ifeq ($(LNKCODE_FLAGS),)
    $(foreach sm,@sm_vers@,$(eval LNKCODE_FLAGS += -arch=sm_$(sm)))
  endif
  CUDA_LIBS  = -lcurand -lcudart -lcuda
  NVCC_CXXFLAGS  = -I../../ -I. $(GENCODE_FLAGS) $(NVCC_EXTRA) -dc -Xcompiler '-fPIC' -Xcompiler $(AM_CPPFLAGS)
endif

AM_CXXFLAGS += $(CUDA_CFLAGS)

libdir      = $(prefix)/lib/user

bin_PROGRAMS = testCrossSection testIonRecomb genIonRecomb testQuantile \
	testCrossCuda

lib_LTLIBRARIES = libuser.la libbar.la libebar.la libsat.la		\
	libdiffuse.la libsatfix.la libdiffrot.la liblogpot.la		\
	libdisk.la libtorque.la libebarp.la librespot.la		\
	librespotn.la librespotorb.la librptest.la libhalo.la		\
	libebars.la libwake.la libsatfixorb.la libdsmc.la libebarn.la	\
	libtreegas.la libperiodic.la libatmos.la libreflect.la		\
	libsatwake.la libpst.la libuserslab.la libtrigger.la		\
	libshear.la libsnheat.la libprofile.la

#
# Convenience library for duplicated source components (static)
#
noinst_LIBRARIES = libcommon.a libdsmc.a libswake.a
noinst_HEADERS = BarForcing.H Perturbation.H UserDiffuse.H		\
	UserResPot.H CircularOrbit.H ResPot.H UserDisk.H UserResPotN.H	\
	Collide.H ResPotOrb.H UserEBar.H UserResPotOrb.H SatFix.H	\
	UserEBarN.H UserTest.H SatFixOrb.H UserEBarP.H UserTreeDSMC.H	\
	CollideLTE.H TwoBodyDiffuse.H UserEBarS.H UserWake.H		\
	EllipForce.H UserBar.H UserHalo.H HeatCool.H UserDSMC.H		\
	UserLogPot.H UserDiffRot.H UserRPtest.H UserAtmos.H		\
	UserReflect.H UserPST.H UserShear.H UserSNheat.H UserProfile.H  \
	CollideIon.H Ion.H interactSelect.H TopBase.H

libcommon_a_SOURCES = ResPot.cc ResPotOrb.cc BarForcing.cc	\
	Perturbation.cc EllipForce.cc SatelliteOrbit.cc		\
	CircularOrbit.cc TwoBodyDiffuse.cc TimeSeriesCoefs.cc	\
	LinearOrbit.cc UnboundOrbit.cc

libdsmc_a_SOURCES = Collide.cc CollideLTE.cc HeatCool.cc		\
	CollideIon.cc Ion.cc interactSelect.cc Elastic.cc TopBase.cc	\
	NTC.cc Quantile.cc QuantileBag.cc phfit2.f			\
	atomic_constants.cc spline.cc $(ADD_CUDA)
libdsmc_a_LIBADD = libdsmc-link.o

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

libdsmc_la_SOURCES = UserDSMC.cc
libdsmc_la_LIBADD  = libcommon.a

libtreegas_la_SOURCES = UserTreeDSMC.cc
libtreegas_la_LDFLAGS = $(BOOST_LDFLAGS) $(CUDA_LDFLAGS)
libtreegas_la_LIBADD  = $(BOOST_FILESYSTEM_LIB) $(BOOST_SERIALIZATION_LIB) libcommon.a libdsmc.a

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

testCrossSection_SOURCES = testCross.cc phfit2.f cudaIon.cu
testCrossSection_LDFLAGS = $(BOOST_LDFLAGS) $(BUILD_FLAGS) @exp_vtk_lib@ $(CUDA_LDFLAGS)
testCrossSection_LDADD   = cudaIon-link.o $(BOOST_FILESYSTEM_LIB) $(BOOST_PROGRAM_OPTIONS_LIB) libdsmc.a -lexputil -lexpgnu @exp_vtk_ldlib@ $(LAPACK_LIBS) $(BLAS_LIBS) $(CUDA_LIBS)

testIonRecomb_SOURCES = testIonRecomb.cc phfit2.f cudaIon.cu
testIonRecomb_LDFLAGS = $(BOOST_LDFLAGS) $(BUILD_FLAGS) @exp_vtk_lib@ $(CUDA_LDFLAGS)
testIonRecomb_LDADD   = cudaIon-link.o $(BOOST_FILESYSTEM_LIB) $(BOOST_PROGRAM_OPTIONS_LIB) libdsmc.a -lexputil -lexpgnu @exp_vtk_ldlib@ $(LAPACK_LIBS) $(BLAS_LIBS) $(CUDA_LIBS)

genIonRecomb_SOURCES = genIonRecomb.cc phfit2.f cudaIon.cu
genIonRecomb_LDFLAGS = $(BOOST_LDFLAGS) $(BUILD_FLAGS) @exp_vtk_lib@ $(CUDA_LDFLAGS)
genIonRecomb_LDADD   = cudaIon-link.o $(BOOST_FILESYSTEM_LIB) $(BOOST_PROGRAM_OPTIONS_LIB) libdsmc.a -lexputil -lexpgnu @exp_vtk_ldlib@ $(LAPACK_LIBS) $(BLAS_LIBS) $(CUDA_LIBS)

testQuantile_SOURCES = testQuantile.cc Quantile.cc 
testQuantile_LDFLAGS = $(BOOST_LDFLAGS) $(BUILD_FLAGS)
testQuantile_LDADD   = $(BOOST_PROGRAM_OPTIONS_LIB)

testCrossCuda_SOURCES = testCrossCuda.cc cudaIon.cu
testCrossCuda_LDFLAGS = $(BOOST_LDFLAGS) $(BUILD_FLAGS) @exp_vtk_lib@  $(CUDA_LDFLAGS)
testCrossCuda_LDADD   = cudaIon-link.o $(BOOST_FILESYSTEM_LIB) $(BOOST_PROGRAM_OPTIONS_LIB) libdsmc.a -lexputil -lexpgnu @exp_vtk_ldlib@ $(LAPACK_LIBS) $(BLAS_LIBS) $(CUDA_LIBS) 

CylindricalDisk.cc:	$(top_srcdir)/utils/SL/CylindricalDisk.cc
			cp $(top_srcdir)/utils/SL/$@ .

SLSphere.cc:		$(top_srcdir)/utils/SL/SLSphere.cc
			cp $(top_srcdir)/utils/SL/$@ .

EllipsoidForce.cc:	$(top_srcdir)/utils/ICs/EllipsoidForce.cc
			cp $(top_srcdir)/utils/ICs/$@ .

# the rule that can recognize *.cu files and compile the *.cu files to
# *.lo and *.o files
.cu.lo :
	$(top_builddir)/cudalt.py $@ $(NVCC) -c $(NVCC_CXXFLAGS) $<

.cu.o:
	$(NVCC) $(CUDA_CFLAGS) $(NVCC_CXXFLAGS) -c -o $@ $<

# trigger cuda compilation to make device linkage
#

# standlone for testing
cudaIon-link.o:	cudaIon.o
	$(NVCC) $(LNKCODE_FLAGS) -Xcompiler '-fPIC' -dlink -o $@ $*.o

# linkage for all cuda objects for library
libdsmc-link.o:	$(OBJ_CUDA)
	$(NVCC) $(LNKCODE_FLAGS) -Xcompiler '-fPIC' -dlink -o $@ $(OBJ_CUDA)

CLEANFILES = 	CylindricalDisk.cc EllipsoidForce.cc SLSphere.cc
