# ===========================================================================
#
# SYNOPSIS
#
#   AX_BOOST_MPI
#
# DESCRIPTION
#
#   Test for MPI library from the Boost C++ libraries. The macro
#   requires a preceding call to AX_BOOST_BASE. Further documentation is
#   available at <http://randspringer.de/boost/index.html>.
#
#   This macro calls:
#
#     AC_SUBST(BOOST_MPI_LIB)
#
#   And sets:
#
#     HAVE_BOOST_MPI
#
# LAST MODIFICATION
#
#   2008-07-17
#
# COPYLEFT
#
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#   Copyright (c) 2008 Michael Tindal
#
#   Cloned from ax_boost_filesystem by Martin Weinberg 2008 to handle
#   boost:mpi configuration
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_BOOST_MPI],
[
	AC_ARG_WITH([boost-mpi],
	AS_HELP_STRING([--with-boost-mpi@<:@=special-lib@:>@],
                   [use the MPI library from boost - it is possible to specify a certain library for the linker
                        e.g. --with-boost-mpi=boost_mpi-gcc-mt ]),
        [
        if test "$withval" = "no"; then
			want_boost="no"
        elif test "$withval" = "yes"; then
            want_boost="yes"
            ax_boost_user_mpi_lib=""
        else
		    want_boost="yes"
        	ax_boost_user_mpi_lib="$withval"
		fi
        ],
        [want_boost="yes"]
	)

	if test "x$want_boost" = "xyes"; then
        AC_REQUIRE([AC_PROG_CC])
		CPPFLAGS_SAVED="$CPPFLAGS"
		CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
		export CPPFLAGS

		LDFLAGS_SAVED="$LDFLAGS"
		LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
		export LDFLAGS

        AC_CACHE_CHECK(whether the Boost::MPI library is available,
					   ax_cv_boost_mpi,
        [AC_LANG_PUSH([C++])
         AC_COMPILE_IFELSE(AC_LANG_PROGRAM([[@%:@include <boost/mpi.hpp>]],
                                   [[using namespace boost::mpi;
                                   boost::mpi::communicator world;
                                   return 0;]]),
            				       ax_cv_boost_mpi=yes, ax_cv_boost_mpi=no)
         AC_LANG_POP([C++])
		])
		if test "x$ax_cv_boost_mpi" = "xyes"; then
			AC_DEFINE(HAVE_BOOST_MPI,,[define if the Boost::MPI library is available])
            BOOSTLIBDIR=`echo $BOOST_LDFLAGS | sed -e 's/@<:@^\/@:>@*//'`
            if test "x$ax_boost_user_mpi_lib" = "x"; then
                for libextension in `ls $BOOSTLIBDIR/libboost_mpi*.{so,a}* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^lib\(boost_mpi.*\)\.so.*$;\1;' -e 's;^lib\(boost_mpi.*\)\.a*$;\1;'` ; do
                     ax_lib=${libextension}
				    AC_CHECK_LIB($ax_lib, exit,
                                 [BOOST_MPI_LIB="-l$ax_lib"; AC_SUBST(BOOST_MPI_LIB) link_mpi="yes"; break],
                                 [link_mpi="no"])
  				done
                if test "x$link_program_options" != "xyes"; then
                for libextension in `ls $BOOSTLIBDIR/boost_mpi*.{dll,a}* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^\(boost_mpi.*\)\.dll.*$;\1;' -e 's;^\(boost_mpi.*\)\.a*$;\1;'` ; do
                     ax_lib=${libextension}
				    AC_CHECK_LIB($ax_lib, exit,
                                 [BOOST_MPI_LIB="-l$ax_lib"; AC_SUBST(BOOST_MPI_LIB) link_mpi="yes"; break],
                                 [link_mpi="no"])
  				done
	            fi
            else
               for ax_lib in $ax_boost_user_mpi_lib boost_mpi-$ax_boost_user_mpi_lib; do
				      AC_CHECK_LIB($ax_lib, exit,
                                   [BOOST_MPI_LIB="-l$ax_lib"; AC_SUBST(BOOST_MPI_LIB) link_mpi="yes"; break],
                                   [link_mpi="no"])
                  done

            fi
			if test "x$link_mpi" != "xyes"; then
				AC_MSG_ERROR(Could not link against $ax_lib !)
			fi
		fi

		CPPFLAGS="$CPPFLAGS_SAVED"
    	LDFLAGS="$LDFLAGS_SAVED"
	fi
])
