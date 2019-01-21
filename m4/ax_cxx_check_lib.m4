AC_DEFUN([AX_CXX_CHECK_LIB],
[m4_ifval([$3], , [AH_CHECK_LIB([$1])])dnl
AS_LITERAL_IF([$1],
	      [AS_VAR_PUSHDEF([ac_Lib], [ac_cv_lib_$1_$2])],
	      [AS_VAR_PUSHDEF([ac_Lib], [ac_cv_lib_$1''_$2])])dnl
AC_CACHE_CHECK([for $2 in -l$1], ac_Lib,
[ac_check_lib_save_LIBS=$LIBS
LIBS="-l$1 $5 $LIBS"
case "$2" 
in *::*::*\(*)
AC_LINK_IFELSE([AC_LANG_PROGRAM([
 namespace `echo "$2" | sed -e "s/::.*//"` 
 { class `echo "$2" | sed -e "s/.*::\\(.*\\)::.*/\\1/" -e "s/(.*//"` 
   { public: int `echo "$2" | sed -e "s/.*:://" -e "/(/!s/..*/&amp;()/"`;
   };
 }
],[`echo "$2" | sed  -e "s/(.*//" -e "s/\\(.*\\)::\\(.*\\)/((\\1*)(0))-&gt;\\2/g"`()])],
	       [AS_VAR_SET(ac_Lib, yes)],
	       [AS_VAR_SET(ac_Lib, no)])
;; *::*::*)
AC_LINK_IFELSE([AC_LANG_PROGRAM([
 namespace `echo "$2" | sed -e "s/::.*//"` 
 { namespace `echo "$2" | sed -e "s/.*::\\(.*\\)::.*/\\1/"` 
   { class `echo "$2" | sed -e "s/.*:://"` 
      { public: `echo "$2" | sed -e "s/.*:://"` ();
      };
   }
 }
],[new $2()])],
	       [AS_VAR_SET(ac_Lib, yes)],
	       [AS_VAR_SET(ac_Lib, no)])
;; *::*\(*)
AC_LINK_IFELSE([AC_LANG_PROGRAM([
 class `echo "$2" | sed -e "s/\\(.*\\)::.*/\\1/" -e "s/(.*//"` 
   { public: int `echo "$2" | sed -e "s/.*:://" -e "/(/!s/..*/&amp;()/"`;
   };
],[`echo "$2" | sed  -e "s/(.*//" -e "s/\\(.*\\)::\\(.*\\)/((\\1*)(0))-&gt;\\2/g"`()])],
	       [AS_VAR_SET(ac_Lib, yes)],
	       [AS_VAR_SET(ac_Lib, no)])
;; *::*)
AC_LINK_IFELSE([AC_LANG_PROGRAM([
 namespace `echo "$2" | sed -e "s/::.*//"` 
 { class `echo "$2" | sed -e "s/.*:://"`
   { public: `echo "$2" | sed -e "s/.*:://"` ();
   };
 }
],[new $2()])],
	       [AS_VAR_SET(ac_Lib, yes)],
	       [AS_VAR_SET(ac_Lib, no)])
;; *)
AC_LINK_IFELSE([AC_LANG_CALL([], [$2])],
	       [AS_VAR_SET(ac_Lib, yes)],
	       [AS_VAR_SET(ac_Lib, no)])
;; esac
LIBS=$ac_check_lib_save_LIBS])
AS_IF([test AS_VAR_GET(ac_Lib) = yes],
      [m4_default([$3], [AC_DEFINE_UNQUOTED(AS_TR_CPP(HAVE_LIB$1))
  LIBS="-l$1 $LIBS"
])],
      [$4])dnl
AS_VAR_POPDEF([ac_Lib])dnl
])# AC_CHECK_LIB
