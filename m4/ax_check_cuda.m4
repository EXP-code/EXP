##### 
#
# SYNOPSIS
#
# AX_CHECK_CUDA
#
# DESCRIPTION
#
# Use existence of 'nvcc', 'cuda.h' and 'libcuda.so' as a proxy for
# the existence of the CUDA SDK.
#
# Locations of these are included in 
#   CUDA_CFLAGS and 
#   CUDA_LDFLAGS.
# Path to nvcc is included as
#   NVCC_PATH
# in config.h.
# 
##### 

AC_DEFUN([AX_CHECK_CUDA], [

# Result flag
ax_cuda_ok=no

# Provide your CUDA path with this		
AC_ARG_WITH(cuda, [  --with-cuda=PREFIX      Prefix of your CUDA installation], [cuda_prefix=$withval], [cuda_prefix="/usr/local/cuda"])

# Setting the prefix to the default if only --with-cuda was given
if test "$cuda_prefix" == "yes"; then
	if test "$withval" == "yes"; then
		cuda_prefix="/usr/local/cuda"
	fi
fi

# Checking for nvcc
AC_MSG_CHECKING([nvcc in $cuda_prefix/bin])
if test -x "$cuda_prefix/bin/nvcc"; then
	AC_MSG_RESULT([found])
	AC_DEFINE_UNQUOTED([NVCC_PATH], ["$cuda_prefix/bin/nvcc"], [Path to nvcc binary])
	ax_cuda_ok=yes

else
	AC_MSG_RESULT([not found!])
fi

if test $ax_cuda_ok = yes; then

# We need to add the CUDA search directories for header and lib searches

# Saving the current flags
ax_save_CFLAGS="${CFLAGS}"
ax_save_LDFLAGS="${LDFLAGS}"

# CUDA variables for make
AC_SUBST([CUDA_CFLAGS])
AC_SUBST([CUDA_LDFLAGS])
AC_SUBST([NVCC])

CUDA_CFLAGS="-I$cuda_prefix/include"
CFLAGS="$CUDA_CFLAGS $CFLAGS"
CUDA_LDFLAGS="-L$cuda_prefix/lib"
LDFLAGS="$CUDA_LDFLAGS $LDFLAGS"
NVCC="$cuda_prefix/bin/nvcc"

# And the header and the lib
AC_MSG_CHECKING([for cuda.h])
AC_CHECK_HEADER([cuda.h], [], ax_cuda_ok=no, [#include <cuda.h>])

if test $ax_cuda_ok = no; then
   AC_MSG_RESULT([Couldn't find cuda.h])
fi

AC_MSG_CHECKING([for libcuda])
AC_CHECK_LIB([cuda], [cuInit], [], ax_cuda_ok=no)

if test $ax_cuda_ok = no; then
   AC_MSG_RESULT([Couldn't find libcuda])
fi

# Returning to the original flags
CFLAGS=${ax_save_CFLAGS}
LDFLAGS=${ax_save_LDFLAGS}

fi

])
