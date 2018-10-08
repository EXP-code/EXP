## Process this file with automake to produce Makefile.in

VPATH = $(top_srcdir)/src:.

CC = mpicc
CXX = mpiCC

if ENABLE_SEQCHECK
ADD_CFLAGS1  = -DSEQCHECK
endif

if ENABLE_DEBUG
ADD_CFLAGS2  = -DDEBUG
endif

if ENABLE_DMALLOC
ADD_CFLAGS3  = -DUSE_DMALLOC -g
ADD_LIBS3    = -ldmallocth -ldmallocthcxx
endif

if ENABLE_GPTL
ADD_CFLAGS4  = -DUSE_GPTL $(gptl_include)
ADD_LIBS4    = $(gptl_lib) -lgptl
endif

AM_CFLAGS    = -D_REENTRANT -fPIC -fopenmp $(ADD_CFLAGS1) $(ADD_CFLAGS2) $(ADD_CFLAGS3) $(ADD_CFLAGS4) $(CUDA_CFLAGS)

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/include -I$(top_srcdir)/exputil -I$(top_srcdir)/expgnu -I$(top_srcdir)/src -I/usr/include/mpi

AM_CXXFLAGS  = -std=c++11 -D_REENTRANT -fopenmp $(ADD_CFLAGS1) $(ADD_CFLAGS2) $(ADD_CFLAGS3) $(ADD_CFLAGS4) @exp_vtk_include@  $(CUDA_CFLAGS)

if ENABLE_CUDA
  ADD_CUDA = cudaParticle.cu cudaSphericalBasis.cu cudaSphere.cu	\
	cudaCylinder.cu cudaEmpOrth.cu cudaComponent.cu NVTX.cc
  ifeq ($(GENCODE_FLAGS),)
    $(foreach sm,@sm_vers@,$(eval GENCODE_FLAGS += -gencode arch=compute_$(sm),code=sm_$(sm)))
  endif
  CUDA_LIBFLAGS  = -lcudart -lnvToolsExt
  NVCC_CXXFLAGS  = -I. $(GENCODE_FLAGS) $(NVCC_EXTRA) -Xcompiler $(AM_CPPFLAGS)
endif

CXXLINK = mpiCC -rdynamic $(CXXFLAGS) $(EXPFLAGS) -o $@

exp_SOURCES = Basis.cc Bessel.cc CBrock.cc Component.cc CBrockDisk.cc	\
     Cube.cc Cylinder.cc EmpOrth9thd.cc ExternalCollection.cc		\
     ExternalForce.cc Hernquist.cc Orient.cc Particle.cc PotAccel.cc	\
     ScatterMFP.cc SphericalBasis.cc AxisymmetricBasis.cc Sphere.cc	\
     TwoCenter.cc EJcom.cc global.cc begin.cc ddplgndr.cc Direct.cc	\
     Shells.cc NoForce.cc end.cc OutputContainer.cc OutPS.cc		\
     OutPSN.cc OutCHKPT.cc Output.cc header.cc expand.cc		\
     externalShock.cc generateRelaxation.cc HaloBulge.cc incpos.cc	\
     incvel.cc ComponentContainer.cc OutAscii.cc OutMulti.cc		\
     OutRelaxation.cc OrbTrace.cc OutDiag.cc OutLog.cc OutCoef.cc	\
     multistep.cc parse.cc Slab.cc SlabSL.cc step.cc sbessj.cc		\
     sbessz.cc tidalField.cc ultra.cc ultrasphere.cc MPL.cc zbrent.cc	\
     OutFrac.cc OutCalbr.cc ParticleFerry.cc pCell.cc chkSlurm.c	\
     chkTimer.cc pHOT.cc VtkPCA.cc $(ADD_CUDA)


bin_PROGRAMS = exp

BOOST_LIBS = $(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB)
COMMON_LIBS = -lexputil -lexpgnu $(BOOST_LIBS)

SYS_LIBS = -lpthread -ldl -lgomp
MPI_LIBS = -lmpi

EXPLIBS  = $(COMMON_LIBS) @LAPACK_LIBS@ $(MPI_LIBS) $(SYS_LIBS) $(LIBS)
LDADD    = $(COMMON_LIBS) @LAPACK_LIBS@ $(ADD_LIBS) $(CUDA_LDFLAGS)
EXPFLAGS = -L$(top_srcdir)/exputil/.libs -L$(top_srcdir)/expgnu/.libs $(BOOST_LDFLAGS)

exp_LDADD = $(EXPLIBS) $(ADD_LIBS3) $(ADD_LIBS4) $(LAPACK_LIBS)		\
	$(BLAS_LIBS) $(LIBS) $(FLIBS) @exp_vtk_ldlib@ @exp_vtk_lib@	\
	$(CUDA_LDFLAGS) $(CUDA_LIBFLAGS)


noinst_HEADERS = Basis.H Bessel.H CBrock.H Component.H Cube.H Cylinder.H     \
	ExternalCollection.H ExternalForce.H Hernquist.H Particle.H          \
	Particle.H PotAccel.H ScatterMFP.H SphericalBasis.H externalShock.H  \
	generateRelaxation.H tidalField.H EmpOrth9thd.h Orient.H OutPS.H     \
	OutMulti.H coef.H exp_thread.h expand.h sltableMP2.h tipsydefs.h     \
        version.h pHOT.H pCell.H pHOT_types.H ParticleFerry.H                \
	pHOT.H pHOT_types.H

# the rule that can recognize *.cu files and compile the *.cu files to
# *.lo and *.o files
.cu.lo:
	$(top_builddir)/cudalt.py $@ $(NVCC) -c $(NVCC_CXXFLAGS) $<

.cu.o :
	$(NVCC) $(CUDA_CFLAGS) $(NVCC_CXXFLAGS) -c -o $@ $<

# special case for f77 for now

sledge.o:	sledge.f
		f77 -g -o sledge.o -c sledge.f

SUBDIRS = user

## src/Makefile.am ends here
