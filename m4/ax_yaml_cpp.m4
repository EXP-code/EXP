#
# SYNOPSIS
#
#   ACX_YAML([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro looks for a version of the YAML library.  The YAML_CPPFLAGS
#   and YAML output variables hold the compile and link flags.
#
#   To link an application with YAML, you should link with:
#
#   	$YAML
#
#   The user may use:
# 
#       --with-yaml-cpp=<flags> --with-yaml-libs=<flags> 
#
#   to manually specify the YAML include and linking flags.
#
#   ACTION-IF-FOUND is a list of shell commands to run if an YAML library is
#   found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it is
#   not found. If ACTION-IF-FOUND is not specified, the default action will
#   define HAVE_YAML and set output variables above.
#
#   This macro requires autoconf 2.50 or later.
#
# LAST MODIFICATION
#
#   2010-05-12
#

AC_DEFUN([AX_YAML], [

ax_yaml_ok=no
ax_yaml_default="-lyaml-cpp"

YAML_CPPFLAGS=""
YAML=""

AC_ARG_WITH(yaml-cpp, [AC_HELP_STRING([--with-yaml-cpp=<flags>], [use YAML preprocessing flags <flags>.  Set to "no" to disable.])])

AC_ARG_WITH(yaml-libs, [AC_HELP_STRING([--with-yaml-libs=<flags>], [use YAML linking flags <flags>.  Set to "no" to disable.])])

if test x"$with_yaml_cpp" != x; then
   if test x"$with_yaml_cpp" != xno; then
      YAML_CPPFLAGS="$with_yaml_cpp"
   else
      ax_yaml_ok=disable
   fi
fi

if test x"$with_yaml_libs" != x; then
   if test x"$with_yaml_libs" != xno; then
      YAML="$with_yaml_libs"
   else
      ax_yaml_ok=disable
   fi
fi

if test $ax_yaml_ok = disable; then
   echo "**** YAML explicitly disabled by configure."
else

   # Save environment

   ax_yaml_save_CXX="$CXX"
   ax_yaml_save_CPPFLAGS="$CPPFLAGS"
   ax_yaml_save_LIBS="$LIBS"

   # Test serial compile and linking

   CPPFLAGS="$CPPFLAGS $YAML_CPPFLAGS"
   LIBS="$YAML $ax_yaml_save_LIBS -lm"

   AC_CHECK_HEADERS([yaml-cpp/yaml.h])

   AX_CXX_CHECK_LIB([yaml-cpp], [YAML::Parser], [ax_yaml_ok=yes;AC_DEFINE(HAVE_YAML,1,[Define if you have the YAML library.])])

   if test $ax_yaml_ok = no; then
      YAML="$ax_yaml_default"
      LIBS="$ax_yaml_default $ax_yaml_save_LIBS -lm"
      AX_CXX_CHECK_LIB([yaml-cpp], [YAML::Parser], [ax_yaml_ok=yes;AC_DEFINE(HAVE_YAML,1,[Define if you have the YAML library.])])
   fi

   # Restore environment

   CXX="$ax_yaml_save_CXX"
   LIBS="$ax_yaml_save_LIBS"
   CPPFLAGS="$ax_yaml_save_CPPFLAGS"

fi

# Define exported variables

AC_SUBST(YAML_CPPFLAGS)
AC_SUBST(YAML)

# Execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
   
if test x"$ax_yaml_ok" = xyes; then
   ifelse([$1],,[echo "**** Enabling support for YAML."],[$1])
else
   ifelse([$2],,[echo "**** YAML not found - disabling support."],[$2])
fi

])

