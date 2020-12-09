dnl ==========================================================================
dnl Process this file with autoconf to produce a configure script.
dnl ==========================================================================
AC_INIT(EXP, 7.612, [mdw@umass.edu], [], [https://mdweinberg@bitbucket.org/mdweinberg/exp.git])

AM_CONFIG_HEADER(config.h)
AC_CONFIG_SRCDIR(src/expand.cc)
AM_INIT_AUTOMAKE
LT_INIT

AC_DEFINE([NEED_MAIN], 1, [Define if a function `main' is needed.])

dnl ==========================================================================
dnl Add info from git repo
dnl ==========================================================================

AC_REVISION([m4_esyscmd_s([git describe --always])])

dnl ==========================================================================
dnl Checks for programs
dnl ==========================================================================

AC_PROG_AWK
AC_PROG_CC(mpicc gcc cc)
AC_PROG_F77(mpif77 mpiF77 g77 f77)
AC_PROG_CPP
AC_PROG_CXX(mpic++ mpiCC g++ c++)
AC_PROG_CXXCPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_LIBTOOL

AC_CONFIG_MACRO_DIR([m4])

AC_LANG_PUSH(C++)

AX_CXX_COMPILE_STDCXX_14([noext],[mandatory])

dnl ==========================================================================
dnl Check for hash
dnl ==========================================================================
AC_CHECK_HEADER(tr1/unordered_map,[AC_DEFINE([HAVE_TR1],[],[Have tr1])],[])
AC_CHECK_HEADER(unordered_map,[AC_DEFINE([HAVE_CXX0X],[],["Have C++0x])],[])
AC_CHECK_HEADER(hash_map,[AC_DEFINE([HAVE_HASH],[],[Have GNU C++])],[])

dnl ==========================================================================
dnl Check for libpng++
dnl ==========================================================================
AC_CHECK_HEADER(png++/png.hpp,[AC_DEFINE([HAVE_LIBPNGPP],[],[Defined if you have the have the png++ wrapper])],
[AC_MSG_RESULT([Couldn't find or include png++/png.hpp.  Please install the png++ package if you want png support. This is optional.])])

dnl ==========================================================================
dnl Check for OpenMP
dnl ==========================================================================

AC_OPENMP

haveomp = no
AC_CHECK_HEADERS(omp.h,[haveomp=[yes]])

if test x$haveomp = xyes; then
   exp_omp_lib="-lgomp"
fi

AC_SUBST(exp_omp_lib)

dnl ==========================================================================
dnl Check for HDF support
dnl ==========================================================================

AX_LIB_HDF5()

dnl ==========================================================================
dnl Look for BOOST libraries
dnl ==========================================================================

AX_BOOST_BASE([1.35.0])
AX_BOOST_SYSTEM
AX_BOOST_FILESYSTEM
AX_BOOST_SERIALIZATION
AX_BOOST_REGEX
AX_BOOST_PROGRAM_OPTIONS


dnl ==========================================================================
dnl Check for YAML-CPP
dnl ==========================================================================
AX_YAML([], [AC_MSG_ERROR([Could not find the yaml-cpp library!  Install the libyaml-cpp-dev package in Ubuntu.])])


dnl ==========================================================================
dnl Look for lapack and blas
dnl ==========================================================================

AX_BLAS
AX_LAPACK

dnl ==========================================================================
dnl Look for CUDA SDK
dnl ==========================================================================

AX_CHECK_CUDA
AM_CONDITIONAL([ENABLE_CUDA], [test "x$CUDA_CFLAGS" != "x"])

dnl ==========================================================================
dnl Added macro to allow list of CUDA compute capabilities
dnl ==========================================================================

sm_vers="35 60 61 70"

AC_ARG_WITH([cuda-capabilities],
	[AC_HELP_STRING([--with-cuda-capabilities=LIST],
	[desired SM compute capabilities (default: 35 60 61 70)])],
	[sm_vers="$withval"])

AC_SUBST(sm_vers)

dnl ==========================================================================
dnl Added macro to allow non-standard location of VTK installation
dnl ==========================================================================
havevtk=no;

AC_ARG_WITH([vtk-lib],
	[AC_HELP_STRING([--with-vtk-lib=DIR],
	[location of VTK libraries is DIR (default: not set)])],
	[VTK_LIB_LOCATION="$withval"; havevtk=yes])

AC_ARG_WITH([vtk-path],
	[AC_HELP_STRING([--with-vtk-path=DIR],
	[location of VTK is DIR (default: not set)])],
	[VTK_PREFIX="$withval"; havevtk=yes])

AC_ARG_WITH([vtk-version],
	[AC_HELP_STRING([--with-vtk-version=VAL],
	[VTK's include directory name is vtk-suffix, e.g. vtk-5.8/. What's the suffix? (Default: not set)])],
	[vtk_suffix=[-$withval]; vtk_version=[$withval]],
	[vtk_suffix=""; vtk_version=5.8])

AC_ARG_WITH([vtk-lib-suffix],
	[AC_HELP_STRING([--with-vtk-lib-suffix],
	[VTK's libraries includes the version as library suffix (e.g. libvtkalglib-6.1). Debian/Ubuntu needs this for vtk>=6.0. (Default: not used)])],
	[vtk_libval=$vtk_suffix],
	[vtk_libval=""])

# Look in some standard places . . .

AC_CHECK_FILE([/usr/include/vtk/vtkMapper.h],
	[havevtk=[yes]; VTK_HEADER_LOCATION=[/usr/include/vtk]])
AC_CHECK_FILE([/usr/include/vtk$vtk_suffix/vtkMapper.h],
	[havevtk=[yes]; VTK_HEADER_LOCATION=[/usr/include/vtk$vtk_suffix]])
AC_CHECK_FILE([/usr/local/include/vtk/vtkMapper.h],
	[havevtk=[yes]; VTK_HEADER_LOCATION=[/usr/local/include/vtk]])
AC_CHECK_FILE([/usr/local/include/vtk$vtk_suffix/vtkMapper.h],
	[havevtk=[yes]; VTK_HEADER_LOCATION=[/usr/local/include/vtk$vtk_suffix]])
AC_CHECK_FILE([$VTK_PREFIX/include/vtk$vtk_suffix/vtkMapper.h],
	[havevtk=[yes]; VTK_HEADER_LOCATION=[$VTK_PREFIX/include/vtk$vtk_suffix]])
AC_CHECK_FILE([/vtkMapper.h],
	[havevtk=[yes]; VTK_HEADER_LOCATION=[$VTK_PREFIX/include/vtk$vtk_suffix]])

AC_SUBST(exp_vtk_include)
AC_SUBST(exp_vtk_lib)
AC_SUBST(exp_vtk_ldlib)

AC_MSG_CHECKING([for VTK])
if test x$havevtk = xyes; then

   AC_DEFINE(HAVE_VTK,,[define if VTK is available])

   if test "$VTK_LIB_LOCATION" != ""; then
      exp_vtk_lib="-L$VTK_LIB_LOCATION -L$VTK_LIB_LOCATION/lib/vtk"
   else
      exp_vtk_lib="-L$VTK_PREFIX/lib/vtk$vtk_suffix -L$VTK_PREFIX/lib"
   fi

   if test "$VTK_HEADER_LOCATION" != ""; then
      exp_vtk_include="-I$VTK_HEADER_LOCATION -I$VTK_HEADER_LOCATION/include/vtk"
   fi
   
   OLD_CFLAGS=$CFLAGS
   OLD_CXXFLAGS=$CXXFLAGS
   OLD_LDFLAGS=$LDFLAGS
   CFLAGS="$exp_vtk_include $CFLAGS"
   CXXFLAGS="$exp_vtk_include $CXXFLAGS"
   LDFLAGS="$exp_vtk_lib $LDFLAGS"

   dnl Compare supplied version to check whether <7.1 or >=7.1 libraries are needed
   CMPV="7.1"
   AS_VERSION_COMPARE(${vtk_version}, ${CMPV}, [VTK_POST_7_1=-1], [VTK_POST_7_1=0],[VTK_POST_7_1=1])

   AS_IF([test "${VTK_POST_7_1}" != "-1"], AC_MSG_RESULT([VTK post 7.1: yes]), AC_MSG_RESULT([VTK post 7.1: no]))

   exp_vtk_ldlib="-lvtkCommonCore$vtk_libval -lvtkCommonDataModel$vtk_libval -lvtkIOXML$vtk_libval"

   if test "${VTK_POST_7_1}" == "-1"; then

      if test x$ac_cv_lib_png_main = xyes; then
      	 exp_vtk_ldlib="$exp_vtk_ldlib -lpng"
      fi

      if test x$ac_cv_lib_z_main = xyes; then
      	 exp_vtk_ldlib="$exp_vtk_ldlib -lz"
      fi

      if test x$ac_cv_lib_freetype_main = xyes; then
      	 exp_vtk_ldlib="$exp_vtk_ldlib -lfreetype"
      fi

      if test x$ac_cv_lib_ftgl_main = xyes; then
      	 exp_vtk_ldlib="$exp_vtk_ldlib -lftgl"
      else
	 exp_vtk_ldlib="$exp_vtk_ldlib -lvtkftgl$vtk_libval"
      fi
   fi

   CFLAGS=$OLD_CFLAGS
   CXXFLAGS=$OLD_CXXFLAGS
   LDFLAGS=$OLD_LDFLAGS

fi
AC_MSG_RESULT([$havevtk])

dnl ==========================================================================
dnl Check for eigen3
dnl ==========================================================================

# Look in some standard places . . .

AC_CHECK_FILE([/usr/include/eigen3/Eigen/Eigen],
	[haveeigen3=[yes]; EIGEN3_LOCATION=[/usr/include/eigen3]])

AC_CHECK_FILE([/usr/local/include/eigen3/Eigen/Eigen],
	[haveeigen3=[yes]; EIGEN3_LOCATION=[/usr/local/include/eigen3]])

AC_ARG_WITH([eigen3],
	[AC_HELP_STRING([--with-eigen3=DIR],
	[location of Eigen3 headers is DIR (default: not set)])],
	[EIGEN3_LOCATION="$withval"; haveeigen3=yes])

AC_MSG_CHECKING([for Eigen3])
if test x$haveeigen3 = xyes; then

   AC_DEFINE(HAVE_EIGEN3,,[define if EIGEN3 is available])

   if test "$EIGEN3_LOCATION" != ""; then
      exp_eigen3_include="-I$EIGEN3_LOCATION"
   fi
fi					

AC_SUBST(exp_eigen3_include)

dnl ==========================================================================
dnl Arguments
dnl ==========================================================================

AC_ARG_ENABLE(gptl,
    [  --enable-gptl           add General Purpose Timing Library],
    enable_gptl=yes,
    enable_gptl=no)

AC_ARG_WITH(gptl,
	[  --with-gptl=DIR         location of GPTL library],
	[   GPTL_LOCATION="$withval"],
	[   GPTL_LOCATION=""])

if test "$GPTL_LOCATION" != ""; then
   LIBS="-L$GPTL_LOCATION/lib $LIBS"
   gptl_include="-I$GPTL_LOCATION/include"
   gptl_lib="-L$GPTL_LOCATION/lib"
fi

AC_SUBST(gptl_include) 
AC_SUBST(gptl_lib) 


AC_ARG_ENABLE(seqcheck,
    [  --enable-seqcheck       add code to check particle sequence on repartition],
    enable_seqcheck=yes,
    enable_seqcheck=no)

AC_ARG_ENABLE(debug,
    [  --enable-debug          compile with verbose output and debugging flags],
    enable_debug=yes,
    enable_debug=no)

AC_ARG_ENABLE(dmalloc,
    [  --enable-dmalloc        link in the dmalloc package for memory debugging],
    enable_dmalloc=yes,
    enable_dmalloc=no)


# ==========================================================================
# CHECKING FOR DSMC
# ==========================================================================
AC_MSG_CHECKING([for DSMC])

AC_ARG_WITH([dsmc],
	[AC_HELP_STRING([--with-dsmc=[yes|no]],
	[compile the DSMC directory (default=no)])],
	[dsmc_flag="$withval"],[dsmc_flag="no"])

dsmc_result=no
have_dsmc=0
dsmc_make=""

AC_SUBST(LDFLAGS_DSMC)
AC_SUBST(INCLUDE_DSMC)
AC_SUBST(SUBDIR_DSMC)
AC_SUBST(EXPADD_DSMC)

LDFLAGS_DSMC=""
INCLUDE_DSMC=""
SUBDIR_DSMC=""
EXPADD_DSMC=""

if test x$dsmc_flag = xyes; then
  dsmc_result=yes
  have_dsmc=1
  dsmc_make="src/DSMC/Makefile"
  LDFLAGS_DSMC="-L\$(top_srcdir)/src/DSMC -lexpdsmc"
  INCLUDE_DSMC="-I\$(top_srcdir)/src/DSMC"
  EXPADD_DSMC="./DSMC/libexpdsmc.a"
  SUBDIR_DSMC="DSMC"
  AC_CONFIG_FILES([src/DSMC/Makefile])
fi
  
AC_MSG_RESULT([$dsmc_result])

AC_DEFINE_UNQUOTED(DSMC_ENABLED, [$have_dsmc], [True if DSMC is enabled])
AM_CONDITIONAL(DSMC_ENABLED, test x$dsmc_result = xyes)


dnl ==========================================================================
dnl Look for doxygen
dnl ==========================================================================
AC_PATH_PROG(DOXYGEN, doxygen, notfound)
AM_CONDITIONAL(DOXYGEN_FOUND, test x$DOXYGEN != xnotfound)

AC_SUBST(SUBDIR_DOXY)

SUBDIR_DOXY=""
doxygen_make=""
if test x$DOXYGEN != xnotfound; then
   SUBDIR_DOXY="doc"
   doxygen_make="doc/Makefile"
fi

dnl ==========================================================================
dnl Checks for libraries.
dnl ==========================================================================

dnl Replace 'main' with a function in -largs:
AC_CHECK_LIB(args, main)

dnl Replace 'main' with a function in -lpthread:
AC_CHECK_LIB(pthread, main)

dnl Replace 'main' with a function in -lstdc:
AC_CHECK_LIB(stdc, main)

dnl Replace 'main' with a function in -lg2c:
AC_CHECK_LIB(g2c, main)

dnl Replace 'main' with a function in -lgfortran:
AC_CHECK_LIB(gfortran, main)

dnl Replace `main' with a function in -lmpi:
AC_CHECK_LIB(mpi, main, , [have_mpi=no])

dnl Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS(malloc.h unistd.h mpi.h fpu_control.h)

dnl Must have an MPI installation
if test x$have_mpi = xno; then
   AC_MSG_FAILURE([Error! You need to have MPI installed!])
fi

dnl Check for slurm, used for checkpoint timing
AC_CHECK_LIB(slurm, main)

dnl ==========================================================================
dnl Checks for typedefs, structures, and compiler characteristics.
dnl ==========================================================================
AC_C_CONST
AC_C_INLINE

dnl ==========================================================================
dnl Set variables for Makefile
dnl ==========================================================================
AM_CONDITIONAL(ENABLE_GPTL, test "$enable_gptl" = "yes")
AM_CONDITIONAL(ENABLE_SEQCHECK, test "$enable_seqcheck" = "yes")
AM_CONDITIONAL(ENABLE_DEBUG, test "$enable_debug" = "yes")
AM_CONDITIONAL(ENABLE_DMALLOC, test "$enable_dmalloc" = "yes")

dnl ==========================================================================
dnl Make the makefiles
dnl ==========================================================================
AC_OUTPUT(expgnu/Makefile exputil/Makefile src/Makefile Makefile $doxygen_make src/user/Makefile $dsmc_make utils/Makefile utils/PhaseSpace/Makefile utils/SL/Makefile utils/Test/Makefile utils/ICs/Makefile utils/Cooling/Makefile utils/Analysis/Makefile)

echo "----------------------------------------------------------------"
echo "Now type: "
echo ""
echo "    make"
echo ""
echo "or"
echo ""
echo "    make CXXFLAGS=\"-O3\" CFLAGS=\"-O3\""
echo ""
echo "to make the package for production work.  The default flags are"
echo "\"-g -O2\".  You may automatically install the executable and"
echo "libraries with the command:"
echo ""
echo "    make install prefix=/home/user"
echo ""
echo "You may wish/need to add the option \"-Wno-literal-suffix\" to"
echo "suppress compilation complaints from g++ generated by the OpenMPI"
echo "headers, which are not c++11 compliant."
echo "----------------------------------------------------------------"
 
if test "x$CUDA_CFLAGS" != "x"; then
echo ""
echo "CUDA compiler flags are specified by the NVCC_EXTRA command line"
echo "flag to \"make\" in case one needs to add various CUDA specific"
echo "options such as \"-G\" for cuda-gdb debugging.  In addition, you"
echo "may specify"
echo ""
echo "    CXX_FLAGS=\"-DO_SINGLE\" ..."
echo ""
echo "to compile with single precision.  Double precision is default. "
echo "The compute capability list is '35 60 61 70' by default and may "
echo "be changed by using --with-cuda-capabilities=LIST flag."
echo "----------------------------------------------------------------"
fi

if test x$DOXYGEN != xnotfound; then
echo "The documentation will be in: doc/html.  You can make the"
echo "documentation before the compilation, if you like."
echo "----------------------------------------------------------------"
else
echo "You need the doxygen package to build the documentation."
echo "You can get this at http://www.stack.nl/~dimitri/doxygen"
echo "----------------------------------------------------------------"
fi
