dnl ==========================================================================
dnl Process this file with autoconf to produce a configure script.
dnl ==========================================================================
AC_INIT(EXP, 7.11)
AM_CONFIG_HEADER(config.h)
AC_CONFIG_SRCDIR(src/expand.cc)
AM_INIT_AUTOMAKE

AC_DEFINE([NEED_MAIN], 1, [Define if a function `main' is needed.])


dnl ==========================================================================
dnl Checks for programs.
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

dnl ==========================================================================
dnl Check for hash
dnl ==========================================================================
AC_CHECK_HEADER(tr1/unordered_map,[AC_DEFINE([HAVE_TR1],[],["Have tr1"])],[])
AC_CHECK_HEADER(unordered_map,[AC_DEFINE([HAVE_CXX0X],[],["Have C++0x"])],[])
AC_CHECK_HEADER(hash_map,[AC_DEFINE([HAVE_HASH],[],["Have GNU C++"])],[])

dnl ==========================================================================
dnl Look for BOOST libraries
dnl ==========================================================================

AX_BOOST_BASE([1.35.0])
AX_BOOST_FILESYSTEM
AX_BOOST_SERIALIZATION
AX_BOOST_REGEX
AX_BOOST_PROGRAM_OPTIONS
AX_BOOST_SYSTEM


dnl ==========================================================================
dnl Look for lapack and blas
dnl ==========================================================================

ACX_BLAS
ACX_LAPACK

dnl ==========================================================================
dnl Added macro to allow non-standard location of VTK installation
dnl ==========================================================================
havevtk=no;
AC_ARG_WITH(vtk-headers,
	[  --with-vtk-headers=DIR	location of VTK headers is DIR],
	[VTK_HEADER_LOCATION="$withval"; havevtk=yes])

AC_ARG_WITH(vtk-lib,
	[  --with-vtk-lib=DIR		location of VTK libraries is DIR],
	[VTK_LIB_LOCATION="$withval"; havevtk=yes])

AC_ARG_WITH(vtk-path,
	[  --with-vtk-path=DIR		location of VTK is DIR],
	[VTK_LIB_LOCATION="$withval"; VTK_HEADER_LOCATION="$withval"; havevtk=yes])

dnl Look in some standard places . . .

if test -f "/usr/include/vtkMapper.h"		||
   test -f "/usr/include/vtk/vtkMapper.h"	||
   test -f "/usr/local/include/vtkMapper.h"	||
   test -f "/usr/local/include/vtk/vtkMapper.h"; then havevtk=yes;
fi

AC_SUBST(exp_vtk_include)
AC_SUBST(exp_vtk_lib)
AC_SUBST(exp_vtk_ldlib)

if test x$havevtk = xyes; then
   if test "$VTK_LIB_LOCATION" != ""; then
      exp_vtk_lib="-L$VTK_LIB_LOCATION -L$VTK_LIB_LOCATION/lib/vtk"
   fi
   if test "$VTK_HEADER_LOCATION" != ""; then
      exp_vtk_include="-I$VTK_HEADER_LOCATION -I$VTK_HEADER_LOCATION/include/vtk"
   fi
   exp_vtk_ldlib="-lvtkHybrid -lvtkRendering -lvtkGraphics -lvtkFiltering -lvtkImaging -lvtkCommon -lvtkIO"
fi

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


dnl ==========================================================================
dnl Look for doxygen
dnl ==========================================================================
AC_PATH_PROG(DOXYGEN, doxygen, notfound)
AM_CONDITIONAL(DOXYGEN_FOUND, test x$DOXYGEN != xnotfound)

dnl ==========================================================================
dnl Checks for libraries.
dnl ==========================================================================

dnl Replace 'main' with a function in -largs:
AC_CHECK_LIB(args, main)

dnl Replace 'main' with a function in -lexputil:
AC_CHECK_LIB(exputil, main)

dnl Replace 'main' with a function in -lpthread:
AC_CHECK_LIB(pthread, main)

dnl Replace 'main' with a function in -lexpgnu:
AC_CHECK_LIB(expgnu, main)

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
AC_CHECK_HEADERS(malloc.h unistd.h mpi.h)

if test x$have_mpi = xno; then
   AC_MSG_FAILURE([Error! You need to have MPI installed!])
fi

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
AC_OUTPUT(expgnu/Makefile exputil/Makefile src/Makefile Makefile doc/Makefile src/user/Makefile utils/Makefile utils/PhaseSpace/Makefile utils/SL/Makefile utils/Test/Makefile utils/ICs/Makefile utils/Cooling/Makefile utils/Analysis/Makefile)

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
echo "----------------------------------------------------------------"
 
if test x$DOXYGEN != xnotfound; then
echo "The documentation will be in: doc/html.  You can make the"
echo "documentation before the compilation, if you like."
echo "----------------------------------------------------------------"
else
echo "You need the doxygen package to build the documentation."
echo "You can get this at http://www.stack.nl/~dimitri/doxygen"
echo "----------------------------------------------------------------"
fi
