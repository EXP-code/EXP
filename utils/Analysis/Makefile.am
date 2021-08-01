## Process this file with automake to produce Makefile.in

CC  = mpicc
CXX = mpiCC

if ENABLE_DMALLOC
ADD_CFLAGS3   = -DUSE_DMALLOC -g
ADD_CXXFLAGS3 = -DUSE_DMALLOC -g
ADD_LIBS3     = -ldmallocth -ldmallocthcxx
endif

AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/src $(CUDA_CFLAGS) $(HDF5_CPPFLAGS) `libpng-config --cflags`
AM_CFLAGS = -D_REENTRANT $(ADD_CFLAGS3)
AM_CXXFLAGS = -D_REENTRANT -DSTANDALONE $(INCLUDES) $(ADD_CXXFLAGS3) $(BOOST_CPPFLAGS) @exp_vtk_include@ @exp_eigen3_include@ `libpng-config --ldflags`

if ENABLE_CUDA
CUDA_LIBFLAGS  = -lcudart
endif

COMMON_LIBS = $(LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) @exp_vtk_ldlib@ @exp_vtk_lib@ 
MPI_LIBS = -lmpi
EXPLIBS = -lexputil -lexp $(COMMON_LIBS) $(ADD_LIBS3) $(MPI_LIBS)	\
	$(CUDA_LIBFLAGS) $(BOOST_RANDOM_LIB)				\
	$(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_SYSTEM_LIB) $(YAML)
EXPFLAGS = -L$(top_srcdir)/exputil/.libs -L$(top_srcdir)/src/.libs	\
$(CUDA_LDFLAGS) @exp_omp_lib@
CXXLINK = mpiCC $(CXXFLAGS) $(EXPFLAGS) -o $@

noinst_HEADERS =

bin_PROGRAMS = diskprof haloprof slabprof psphisto gas2d gas2dcyl	\
makegrid pspbox eofpeek pcatest mssaprofD mssaprofH readcoefs		\
readcoefs2 pcatest sphprof yaml_diff diskfreqs pspxval_sph		\
pspxval_cyl pspxval_cyl2 diskeof pspKL_sph pspKL_cyl kdtest		\
readcoefs2F pspxval_sph2 diskcoef halocoef

diskcoef_SOURCES = diskprof_coef.cc Coefs.cc
diskcoef_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS)
diskcoef_LDADD   = $(EXPLIBS) $(YAML) @exp_vtk_ldlib@ @exp_vtk_lib@ $(HDF5_LIBS)

halocoef_SOURCES = haloprof_coef.cc Coefs.cc SphereSL.cc
halocoef_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS) 
halocoef_LDADD   = $(EXPLIBS) $(YAML) @exp_vtk_ldlib@ @exp_vtk_lib@ $(HDF5_LIBS)


haloprof_SOURCES = haloprof2.cc SphereSL.cc localmpi.cc	\
	ParticleReader.cc
haloprof_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS)
haloprof_LDADD = $(EXPLIBS) $(HDF5_LIBS)

slabprof_SOURCES =  slabprof.cc ParticleReader.cc
slabprof_LDFLAGS =  $(EXPFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS)
slabprof_LDADD   = $(COMMON_LIBS) $(EXPLIBS) $(CUDA_LIBFLAGS) $(HDF5_LIBS)

psphisto_SOURCES = psphisto.cc ParticleReader.cc
psphisto_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS)
psphisto_LDADD   = $(EXPLIBS) $(HDF5_LIBS)

gas2d_SOURCES =  gas2d.cc ParticleReader.cc
gas2d_LDFLAGS =  $(EXPFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS)
gas2d_LDADD   =  $(EXPLIBS) $(HDF5_LIBS)

gas2dcyl_SOURCES = gas2dcyl.cc ParticleReader.cc
gas2dcyl_LDFLAGS = $(EXPFLAGS)  $(BOOST_LDFLAGS) $(HDF5_LDFLAGS)
gas2dcyl_LDADD   = $(EXPLIBS) $(HDF5_LIBS)

makegrid_SOURCES = makegrid.cc ParticleReader.cc
makegrid_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS)
makegrid_LDADD = $(EXPLIBS) $(BOOST_PROGRAM_OPTIONS_LIB)	\
	@exp_vtk_ldlib@ @exp_vtk_lib@ $(HDF5_LIBS)

pspbox_SOURCES = pspbox.cc ParticleReader.cc
pspbox_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS) 
pspbox_LDADD   = $(EXPLIBS) $(BOOST_PROGRAM_OPTIONS_LIB) $(HDF5_LIBS)

eofpeek_SOURCES = eofpeek.cc

diskprof_SOURCES = diskprof2.cc ParticleReader.cc
diskprof_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS) $(HDF5_LDFLAGS)
diskprof_LDADD = $(EXPLIBS) $(YAML) $(BOOST_PROGRAM_OPTIONS_LIB)	\
	$(HDF5_LIBS) @exp_vtk_ldlib@ @exp_vtk_lib@ $(HDF5_LIBS)

diskfreqs_SOURCES = ParticleReader.cc diskfreqs.cc Coefs.cc SphereSL.cc
diskfreqs_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS) 
diskfreqs_LDADD = $(EXPLIBS) $(YAML) $(BOOST_PROGRAM_OPTIONS_LIB)	\
	@exp_vtk_ldlib@ @exp_vtk_lib@ $(HDF5_LIBS)

mssaprofD_SOURCES  = ParticleReader.cc mssaprof_disk.cc
mssaprofD_LDFLAGS  = $(EXPFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS) 
mssaprofD_LDADD = $(EXPLIBS) $(YAML) $(BOOST_PROGRAM_OPTIONS_LIB)	\
	@exp_vtk_ldlib@ @exp_vtk_lib@ $(HDF5_LIBS)

mssaprofH_SOURCES  = ParticleReader.cc mssaprof_halo.cc SphereSL.cc
mssaprofH_LDFLAGS  = $(EXPFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS) 
mssaprofH_LDADD = $(EXPLIBS) $(BOOST_PROGRAM_OPTIONS_LIB)	\
	@exp_vtk_ldlib@ @exp_vtk_lib@ $(HDF5_LIBS)

pcatest_SOURCES = ParticleReader.cc pcatest.cc
pcatest_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS) 
pcatest_LDADD = $(EXPLIBS) $(YAML) $(BOOST_PROGRAM_OPTIONS_LIB)	\
	@exp_vtk_ldlib@ @exp_vtk_lib@ $(HDF5_LIBS)

readcoefs_SOURCES = readcoefs.cc Coefs.cc
readcoefs_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS) 
readcoefs_LDADD   = $(BOOST_PROGRAM_OPTIONS_LIB) $(YAML) $(HDF5_LIBS)

readcoefs2_SOURCES = readcoefs2.cc Coefs.cc
readcoefs2_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS) 
readcoefs2_LDADD   = $(BOOST_PROGRAM_OPTIONS_LIB) $(YAML) $(HDF5_LIBS)

readcoefs2F_SOURCES = readcoefs2F.cc Coefs.cc
readcoefs2F_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS) 
readcoefs2F_LDADD   = $(BOOST_PROGRAM_OPTIONS_LIB) $(YAML) $(HDF5_LIBS)

sphprof_SOURCES = ParticleReader.cc sphprof.cc SphereSL.cc
sphprof_LDFLAGS   = $(EXPFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS)
sphprof_LDADD     = $(EXPLIBS) $(HDF5_LIBS)

pspxval_sph_SOURCES = cross_validation_sph.cc SphereSL.cc	\
			ParticleReader.cc
pspxval_sph_LDFLAGS  = $(EXPFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS)
pspxval_sph_LDADD    = $(EXPLIBS) $(HDF5_LIBS)

pspxval_sph2_SOURCES = cross_validation_sph2.cc SphereSL.cc	\
			ParticleReader.cc
pspxval_sph2_LDFLAGS  = $(EXPFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS)
pspxval_sph2_LDADD    = $(EXPLIBS) $(HDF5_LIBS)

pspxval_cyl_SOURCES = cross_validation_cyl.cc ParticleReader.cc	\
			dtorh1.f rout.f
pspxval_cyl_LDFLAGS  = $(EXPFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS)
pspxval_cyl_LDADD    = $(EXPLIBS) $(HDF5_LIBS)

pspxval_cyl2_SOURCES = cross_validation_cyl2.cc ParticleReader.cc
pspxval_cyl2_LDFLAGS  = $(EXPFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS)
pspxval_cyl2_LDADD    = $(EXPLIBS) $(HDF5_LIBS)

pspKL_sph_SOURCES = KL_sph.cc SphereSL.cc ParticleReader.cc
pspKL_sph_LDFLAGS  = $(EXPFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS)
pspKL_sph_LDADD    = $(EXPLIBS) $(HDF5_LIBS)

pspKL_cyl_SOURCES = KL_cyl.cc ParticleReader.cc
pspKL_cyl_LDFLAGS  = $(EXPFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS)
pspKL_cyl_LDADD    = $(EXPLIBS) $(HDF5_LIBS)

diskeof_SOURCES = diskeof.cc ParticleReader.cc
diskeof_LDFLAGS  = $(EXPFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS)
diskeof_LDADD    = $(EXPLIBS) $(HDF5_LIBS)

yaml_diff_SOURCES = yaml_diff.cc
yaml_diff_LDFLAGS = $(BOOST_LDFLAGS) $(HDF5_LDFLAGS)
yaml_diff_LDADD   = $(BOOST_PROGRAM_OPTIONS_LIB) $(YAML)

kdtest_SOURCES = KDcyltest.cc ParticleReader.cc
kdtest_LDFLAGS = $(EXPFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS)
kdtest_LDADD   = $(EXPLIBS) $(BOOST_PROGRAM_OPTIONS_LIB) $(YAML) $(HDF5_LIBS)

# Install and distribute these scripts

dist_bin_SCRIPTS = parse_outfrac ion_collide.py ion_energy.py	\
	ion_rate.py ion_frac.py ion_coll_number.py

CLEANFILES = header.cc header.H

## utils/ICs/Makefile.am ends here


