#!/bin/sh

if test ! -f install-sh ; then touch install-sh ; fi

NPROC=`cat /proc/cpuinfo | grep processor | wc -l`

MAKE=`which gnumake`
if test ! -x "$MAKE" ; then MAKE=`which gmake` ; fi
if test ! -x "$MAKE" ; then MAKE=`which make` ; fi
HAVE_GNU_MAKE=`$MAKE --version|grep -c "Free Software Foundation"`

if test "$HAVE_GNU_MAKE" != "1"; then
echo !!!! Warning: not tested with non Gnu-Make $MAKE
else
echo Found GNU Make at $MAKE ... good.
fi

echo This script runs git, cmake and make...
echo You did remember necessary arguments for cmake, right?

if test ! -x `which cmake`
then echo you need cmake to generate Makefiles
fi

mkdir -p build
cd build
CC=gcc CXX=g++ cmake -DCUDA_USE_STATIC_CUDA_RUNTIME=off -DCMAKE_CUDA_ARCHITECTURES="72" -DENABLE_CUDA=1 -DEigen3_DIR=$EIGEN_BASE/share/eigen3/cmake -Wno-dev ..
make -j $((2*$(NPROC)))

