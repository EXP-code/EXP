## Process this file with automake to produce Makefile.in

CC  = mpicc
CXX = mpiCC

if ENABLE_DMALLOC
ADD_CFLAGS3 = -DUSE_DMALLOC -g
ADD_CXXFLAGS3 = -DUSE_DMALLOC -g
ADD_LIBS3 = -ldmallocth -ldmallocthcxx
endif

AM_CFLAGS = -D_REENTRANT $(ADD_CFLAGS1) $(ADD_CFLAGS2) $(ADD_CFLAGS3)
AM_CXXFLAGS = -fPIC -std=c++11 -D_REENTRANT -DSTANDALONE $(OPENMP_CXXFLAGS) $(ADD_CXXFLAGS1) $(ADD_CXXFLAGS2) $(ADD_CXXFLAGS3) @exp_vtk_include@ $(BOOST_CPPFLAGS) $(OPENMP_CXXFLAGS)
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/expgnu -I$(top_srcdir)/src -I$(top_srcdir)/utils/SL $(CUDA_CFLAGS)

if ENABLE_CUDA
CUDA_LIBFLAGS  = -lcudart
endif

COMMON_LIBS = -lexputil -lexpgnu $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) @exp_vtk_ldlib@ @exp_vtk_lib@
MPI_LIBS = -lmpi

BOOST_LIBS = $(BOOST_SYSTEM_LIB) $(BOOST_FILESYSTEM_LIB) $(BOOST_PROGRAM_OPTIONS_LIB)

EXPLIBS = $(MPI_LIBS) $(COMMON_LIBS) $(ADD_LIBS3) $(CUDA_LIBFLAGS) @exp_omp_lib@
EXPFLAGS = -L$(top_srcdir)/exputil/.libs -L$(top_srcdir)/expgnu/.libs $(CUDA_LDFLAGS)
CXXLINK = mpiCC $(CXXFLAGS) $(EXPFLAGS) -o $@

noinst_HEADERS =

bin_PROGRAMS = gendisk shrinkics gendisk2 gsphere pstmod gendisk3	\
	empinfo empdump eofbasis makeIon phIon phConv eofcomp testcoefs \
	testcoefs2 testdeval forcetest hdf52accel forcetest2

shrinkics_SOURCES = shrinkics.cc
shrinkics_LDADD = $(EXPLIBS)

gendisk_SOURCES = initial.cc DiskHalo.cc EmpCylSL.cc \
	SphericalSL.cc legendre.cc Particle.cc SParticle.cc \
	ddplgndr.cc sinecosine.cc AddDisk.cc localmpi.cc
gendisk_LDADD = $(COMMON_LIBS) $(EXPLIBS) $(YAML) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

gendisk2_SOURCES =  initial2.cc DiskHalo2.cc EmpCylSL.cc \
	SphericalSL.cc legendre.cc Particle.cc SParticle.cc \
	ddplgndr.cc sinecosine.cc AddDisk.cc localmpi.cc
gendisk2_LDFLAGS =  $(EXPFLAGS) $(BOOST_LDFLAGS)
gendisk2_LDADD = $(COMMON_LIBS) $(EXPLIBS) $(YAML) $(BOOST_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

gendisk3_SOURCES =  initial3.cc DiskHalo3.cc EmpCylSL.cc \
	CylDisk.cc SphericalSL.cc legendre.cc Particle.cc SParticle.cc \
	ddplgndr.cc sinecosine.cc AddDisk.cc localmpi.cc
gendisk3_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS)
gendisk3_LDADD = $(COMMON_LIBS) $(EXPLIBS) $(YAML) $(BOOST_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

gsphere_SOURCES = bonnerebert.cc
gsphere_LDADD = $(COMMON_LIBS) $(EXPLIBS)

pstmod_SOURCES = pst_model.cc CylindricalDisk.cc EllipsoidForce.cc SLSphere.cc
pstmod_LDFLAGS = $(EXPFLAGS)
pstmod_LDADD = $(COMMON_LIBS) $(EXPLIBS)

empinfo_SOURCES = empinfo.cc EmpCylSL.cc localmpi.cc
empinfo_LDFLAGS = $(EXPFLAGS)
empinfo_LDADD = $(COMMON_LIBS) $(EXPLIBS) $(YAML) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

empdump_SOURCES = empdump.cc EmpCylSL.cc localmpi.cc
empdump_LDFLAGS = $(EXPFLAGS)
empdump_LDADD = $(COMMON_LIBS) $(EXPLIBS) $(YAML) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

eofbasis_SOURCES = eof_basis.cc EmpCylSL.cc localmpi.cc
eofbasis_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS)
eofbasis_LDADD = $(COMMON_LIBS) $(EXPLIBS) $(YAML) $(BOOST_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

eofcomp_SOURCES = eof_compare.cc EmpCylSL.cc localmpi.cc
eofcomp_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS)
eofcomp_LDADD = $(COMMON_LIBS) $(EXPLIBS) $(YAML) $(BOOST_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

makeIon_SOURCES = makeIonIC.cc Particle.cc Particle.H Species.H \
	atomic_constants.H atomic_constants.cc
makeIon_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS)
makeIon_LDADD = $(COMMON_LIBS) $(EXPLIBS) $(LAPACK_LIBS) $(BLAS_LIBS)	\
	$(BOOST_LIBS) $(LIBS) $(FLIBS) $(YAML)

phIon_SOURCES = ph_ion.cc
phIon_LDFLAGS = $(COMMON_LIBS) $(EXPFLAGS) $(BOOST_LDFLAGS)
phIon_LDADD = $(BLAS_LIBS) $(BOOST_LIBS) $(LIBS) $(FLIBS)

phConv_SOURCES = ph_conv.cc
phConv_LDFLAGS = $(COMMON_LIBS) $(EXPFLAGS) $(BOOST_LDFLAGS)
phConv_LDADD = $(BLAS_LIBS) $(BOOST_LIBS) $(LIBS)

testcoefs_SOURCES = check_coefs.cc EmpCylSL.cc SphericalSL.cc	\
	legendre.cc sinecosine.cc localmpi.cc
testcoefs_LDFLAGS = $(COMMON_LIBS) $(EXPFLAGS) $(BOOST_LDFLAGS)
testcoefs_LDADD = $(BLAS_LIBS) $(BOOST_LIBS) $(YAML) $(LIBS) @exp_omp_lib@

testcoefs2_SOURCES = check_coefs2.cc DiskEval.cc EmpCylSL.cc	\
	SphericalSL.cc legendre.cc sinecosine.cc localmpi.cc
testcoefs2_LDFLAGS = $(COMMON_LIBS) $(EXPFLAGS) $(YAML) $(BOOST_LDFLAGS)
testcoefs2_LDADD = $(BLAS_LIBS) $(BOOST_LIBS) $(LIBS) @exp_omp_lib@

testdeval_SOURCES = test_deval.cc DiskEval.cc
testdeval_LDFLAGS = $(COMMON_LIBS) $(EXPFLAGS) $(BOOST_LDFLAGS)
testdeval_LDADD = $(BLAS_LIBS) $(BOOST_LIBS) $(LIBS) @exp_omp_lib@

forcetest_SOURCES = force_test.cc DiskEval.cc
forcetest_LDFLAGS = $(COMMON_LIBS) $(EXPFLAGS) $(BOOST_LDFLAGS)
forcetest_LDADD = $(BLAS_LIBS) $(BOOST_LIBS) $(LIBS) @exp_omp_lib@

forcetest2_SOURCES = force_test2.cc DiskEval.cc
forcetest2_LDFLAGS = $(COMMON_LIBS) $(EXPFLAGS) $(BOOST_LDFLAGS)
forcetest2_LDADD = $(BLAS_LIBS) $(BOOST_LIBS) $(LIBS) @exp_omp_lib@

hdf52accel_SOURCES = hdf52accel.cc
hdf52accel_CXXFLAGS = $(AM_CXXFLAGS) $(HDF5_CPPFLAGS)
hdf52accel_LDFLAGS = $(HDF5_LDFLAGS) $(BOOST_LDFLAGS)
hdf52accel_LDADD = $(BOOST_LIBS) $(LIBS) @exp_omp_lib@ $(HDF5_LDFLAGS) -lhdf5_cpp -lhdf5_serial

initial.cc:		Particle.H Particle.h

initial2.cc:		Particle.H Particle.h

initial3.cc:		Particle.H Particle.h

empinfo.cc:		Particle.H Particle.h

empdump.cc:		Particle.H Particle.h

makeIonIC.cc:		Particle.H Particle.h atomic_constants.H

CylindricalDisk.cc:	$(top_srcdir)/utils/SL/CylindricalDisk.cc
			cp $(top_srcdir)/utils/SL/$@ .

SLSphere.cc:		$(top_srcdir)/utils/SL/SLSphere.cc
			cp $(top_srcdir)/utils/SL/$@ .

EmpCylSL.cc:		$(top_srcdir)/src/EmpCylSL.cc Particle.h
			cp $(top_srcdir)/src/$@ .

localmpi.cc:		$(top_srcdir)/exputil/localmpi.cc
			cp $(top_srcdir)/exputil/$@ .

Particle.cc:		$(top_srcdir)/src/Particle.cc Particle.H
			cp $(top_srcdir)/src/$@ .

Particle.H:		$(top_srcdir)/src/Particle.H
			cp $(top_srcdir)/src/$@ .

Particle.h:		Particle.H
			rm -rf Particle.h
			ln -s Particle.H Particle.h

Species.H:		$(top_srcdir)/src/Species.H
			cp $(top_srcdir)/src/$@ .

makeIonIC.cc:		atomic_constants.H

atomic_constants.cc:	$(top_srcdir)/src/DSMC/atomic_constants.cc atomic_constants.H
			cp $(top_srcdir)/src/DSMC/$@ .

atomic_constants.H:	$(top_srcdir)/src/DSMC/atomic_constants.H
			cp $(top_srcdir)/src/DSMC/$@ .


CLEANFILES = EmpCylSL.cc localmpi.cc CylindricalDisk.cc	\
	SLSphere.cc Particle.cc Particle.H Particle.h Species.H \
	atomic_constants.H atomic_constants.cc

## utils/ICs/Makefile.am ends here


