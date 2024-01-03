# Configuring and building EXP

We are now using git submodules to provide `yaml-cpp`, which is not
common in the HPC environments.  So, from the top-level directory, do
the following:

```
   git submodule update --init --recursive
```

This will install `yaml-cpp` in the `extern` directory.  The png++ C++
wrappers to the png library are also installed in `extern`.  Note: png
support is optional.

## EXP uses CMake for building a configuration

I recommend building "out of source" to allow for multiple
configurations.  This allows one to have build various versions
available from the same source, such as `Release` and `Debug`.  To
begin, make a build directory and change to that:

```
   mkdir -p build
   cd build
```

CMake is designed to detect commonly used utilities and libraries
automatically, but sometimes needs help and hints.  For example, if
CMake does not find a library, you can add the location of the library
to the `CMAKE_PREFIX_PATH` environment variable.  I will give
customization examples of Eigen3 and CUDA below.  For example, one
often needs to use the CMake `Eigen3_DIR` variable set to the install
location so CMake can find the Eigen3 CMake rules.

If you do not want png support or do not have png installed, then you
can build without by leaving the `ENABLE_PNG=NO` flag in CMake; this
is the default.  If you want png support, add `-DENABLE_PNG=YES` to
the CMake call.

Similarly, CMake will do its best to find your FFTW package by
default.  If that fails, and if you need or want FFTW and your FFTW is
installed in an unusual location, you can define that location using
`-DFFTW_ROOT=/path/to/fftw/location`.  Or you can provide the location
in the `FFTWDIR` environment variable.

The CMake default for the install location is `/usr/local`.
Generally, the install location will need to be changed in the example
below.  E.g. I would use `-DCMAKE_INSTALL_PREFIX=/home/mdw_umass_edu`
on the UMass Unity cluster to install in my home directory.

## EXP options

There are a number of EXP-specific options that control the build.
The most important of these are:

1. EXP has two independent interfaces: the EXP N-body code, `exp`, and
   the Python bindings to the core library, `pyEXP`.  Both are `ON` by
   default.  If you only want `pyEXP`, add `-DENABLE_NBODY=NO` to your
   CMake call.

2. CUDA support is off by default.  If you want it, add
   `-DENABLE_CUDA=YES` to your CMake call.  This is only useful for
   the N-body code.  CUDA support for the Python bindings is planned.

3. A subset of user modules are ON by default.  These can be disabled
   by `-DENABLE_USER=NO`.  This compiles a subset of user modules.
   The full set can be compiled using `-DENABLE_USER_ALL=YES`.

4. The `BUILD_DOCS` variable controls the build of the Doxygen
   'online' manual.  This provides details of all EXP classes at the
   header level, useful if you are interested in developing.  To get
   the docs, set `-DBUILD_DOCS=ON` and Doxygen is found, you'll get
   the HTML docs installed in shared/EXP

5. Slurm is used at many HPC centers for job control.  The N-body code
   is Slurm aware and can smoothly end your job when the next time
   step will exceed the wall-clock limit.  If you would like Slurm to
   be aware of this, set `-DENABLE_SLURM=ON`

6. EXP can use both VTK and PNG for diagnostic graphical output.
   Neither are essential.  Graphical output has an ASCII fall back
   that is designed to be read by Numpy and visualized with Pyplot.
   If you would like either of the native graphic support enabled set
   `-DENABLE_PNG=ON` and/or `-DENABLE_VTK=ON`.


## Configuring Eigen3

Some installations provide an `EIGEN_BASE` environment variable that
locates the install directory what contains 'include' and 'share'.
Alternatively, replace `EIGEN_BASE` with that path or set `EIGEN_BASE`
manually.

## Configuring CUDA support

Different versions of CMake seem to treat CUDA architecture
specification differently.  Since CMake version >= 3.18, the default
is the lowest `nvcc`-supported compute capability for CUDA.  You are
therefore **required** to set `CUDAARCHS` to a semi-colon-separated list
of custom compute capabilites or enter your desired string manually
using `ccmake`.

The CUDA real size is double (real\*8) by default.  You can configure
EXP to use real\*4 with the `-DENABLE_CUDA_SINGLE=on` flag to CMake.  This
will save some GPU memory if you are close to your hardware limit, but
I don't recommend this generally.

The CUDA particle structure can carry a fixed number real attributes.
This is configurable at compile time using the -DCUDA_EXP_DATTRIB=X flag
for X attributes.  This is 4 by default.

Putting these together so far, your CMake call would be:

```
export CUDAARCHS="75;80;86"
cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_CUDA=YES -DENABLE_USER=YES -DEigen3_DIR=$EIGEN_BASE/share/eigen3/cmake -DCMAKE_INSTALL_PREFIX=/home/user -Wno-dev ..
````

## Configuring without CUDA

Without CUDA support, the same CMake call becomes:
```
cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_USER=YES -DEigen3_DIR=$EIGEN_BASE/share/eigen3/cmake -DCMAKE_INSTALL_PREFIX=/home/user -Wno-dev ..
```

Many users will like configuring with one the CMake GUI tools, such as
`ccmake` or `cmake-gui` instead of the command-line `cmake`.  The GUI
will allow you to change the parameters interactively and display the
help info for each parameter.  For example:

```
   ccmake ..
```

and then enter your preferred build type and other options
interactively.  This provides a nice view of the configuration as a
bonus.  If you are not familiar either of these, I recommend `ccmake`
rather than the Qt `cmake-gui` for portability.

You can use CMake build type Debug for debugging and etc. or use None
or 'empty' and set your own CFLAGS and CXXFLAGS.  See the CMake manual
for additional details.  EXP also has a few custom build types for
debugging and profiling that use Google's code sanitizer package.

## Building

Now you are ready to build.

Make the package `make -j N`. Here, `N` is the number of jobs to run
simultaneously.  I often use N=2*<number of cores> to take advantage
of hyper threading.  So a typical build command would be:

```
   make -j 8
```
Depending on your configuration, this may take some minutes.

## Installing

Finally, install to the target location.  You can select the target
install location using the CMAKE_INSTALL_PREFIX variable in CMake.
Then:

```
   make install
```

## A note on multiple builds

CMake workflow is designed to permit multiple build types (e.g. Debug,
Release) in separate directories with the same source.  However, EXP
generates a 'config.h' based on the available packages. For example,
if you want to generate a build hierarchy like this:

```
   build/
   build/debug
   build/release
```
with the different build types alone, the multiple build strategy
will work perfectly.

## Testing

EXP using `CTest` for basic unit testing.  The tests are configured by
default, but if you really don't want them, you can change the
`ENABLE_TESTS` variable in CMake (e.g. `-DENABLE_TESTS=OFF`).  There
is no overhead for leaving them enabled.

Once you have successfully build and installed EXP, you can run `make
test` or the `ctest` command to run the unit tests.  They will take
about 1 minute to run.  I recommend this as a first step to ensure
that everything is working as expected.
