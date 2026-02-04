# Install script for directory: /home/runner/work/EXP/EXP/expui

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libexpui.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libexpui.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libexpui.so"
         RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu/openmpi/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/runner/work/EXP/EXP/_codeql_build_dir/expui/libexpui.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libexpui.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libexpui.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libexpui.so"
         OLD_RPATH "/home/runner/work/EXP/EXP/_codeql_build_dir/exputil:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/runner/work/EXP/EXP/_codeql_build_dir/extern/yaml-cpp:"
         NEW_RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu/openmpi/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libexpui.so")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/EXP" TYPE FILE FILES
    "/home/runner/work/EXP/EXP/expui/Coefficients.H"
    "/home/runner/work/EXP/EXP/expui/CoefStruct.H"
    "/home/runner/work/EXP/EXP/expui/BasisFactory.H"
    "/home/runner/work/EXP/EXP/expui/BiorthBasis.H"
    "/home/runner/work/EXP/EXP/expui/BiorthBess.H"
    "/home/runner/work/EXP/EXP/expui/FieldGenerator.H"
    "/home/runner/work/EXP/EXP/expui/UnitValidator.H"
    "/home/runner/work/EXP/EXP/_codeql_build_dir/config_exp.h"
    "/home/runner/work/EXP/EXP/expui/../include/DiskDensityFunc.H"
    "/home/runner/work/EXP/EXP/expui/../include/ParticleReader.H"
    "/home/runner/work/EXP/EXP/expui/../include/header.H"
    "/home/runner/work/EXP/EXP/expui/../include/Particle.H"
    "/home/runner/work/EXP/EXP/expui/../include/gadget.H"
    "/home/runner/work/EXP/EXP/expui/../include/tipsy.H"
    "/home/runner/work/EXP/EXP/expui/../include/BiorthCube.H"
    "/home/runner/work/EXP/EXP/expui/../include/BiorthCyl.H"
    "/home/runner/work/EXP/EXP/expui/../include/sltableMP2.H"
    "/home/runner/work/EXP/EXP/expui/../include/SLGridMP2.H"
    "/home/runner/work/EXP/EXP/expui/../include/YamlCheck.H"
    "/home/runner/work/EXP/EXP/expui/../include/EmpCylSL.H"
    "/home/runner/work/EXP/EXP/expui/../include/localmpi.H"
    "/home/runner/work/EXP/EXP/expui/../include/exputils.H"
    "/home/runner/work/EXP/EXP/expui/../include/massmodel.H"
    "/home/runner/work/EXP/EXP/expui/../include/interp.H"
    "/home/runner/work/EXP/EXP/expui/../include/orbit.H"
    "/home/runner/work/EXP/EXP/expui/../include/StringTok.H"
    "/home/runner/work/EXP/EXP/expui/../include/OrthoFunction.H"
    "/home/runner/work/EXP/EXP/expui/../include/gaussQ.H"
    "/home/runner/work/EXP/EXP/expui/../include/PseudoAccel.H"
    "/home/runner/work/EXP/EXP/expui/../include/QuadLS.H"
    "/home/runner/work/EXP/EXP/expui/../include/MonotCubicInterpolator.H"
    "/home/runner/work/EXP/EXP/expui/../include/QPDistF.H"
    "/home/runner/work/EXP/EXP/expui/../include/DiskWithHalo.H"
    "/home/runner/work/EXP/EXP/expui/../include/EmpCyl2d.H"
    "/home/runner/work/EXP/EXP/expui/../include/EXPmath.H"
    "/home/runner/work/EXP/EXP/expui/../include/libvars.H"
    "/home/runner/work/EXP/EXP/expui/../include/Timer.H"
    "/home/runner/work/EXP/EXP/expui/../include/coef.H"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libexpui.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libexpui.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libexpui.so"
         RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu/openmpi/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/runner/work/EXP/EXP/_codeql_build_dir/expui/libexpui.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libexpui.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libexpui.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libexpui.so"
         OLD_RPATH "/home/runner/work/EXP/EXP/_codeql_build_dir/exputil:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/runner/work/EXP/EXP/_codeql_build_dir/extern/yaml-cpp:"
         NEW_RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu/openmpi/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libexpui.so")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/runner/work/EXP/EXP/_codeql_build_dir/expui/CMakeFiles/expui.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/nativetoh5" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/nativetoh5")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/nativetoh5"
         RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu/openmpi/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/runner/work/EXP/EXP/_codeql_build_dir/expui/nativetoh5")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/nativetoh5" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/nativetoh5")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/nativetoh5"
         OLD_RPATH "/home/runner/work/EXP/EXP/_codeql_build_dir/expui:/home/runner/work/EXP/EXP/_codeql_build_dir/exputil:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/runner/work/EXP/EXP/_codeql_build_dir/extern/yaml-cpp:"
         NEW_RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu/openmpi/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/nativetoh5")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/runner/work/EXP/EXP/_codeql_build_dir/expui/CMakeFiles/nativetoh5.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/h5compare" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/h5compare")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/h5compare"
         RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu/openmpi/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/runner/work/EXP/EXP/_codeql_build_dir/expui/h5compare")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/h5compare" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/h5compare")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/h5compare"
         OLD_RPATH "/home/runner/work/EXP/EXP/_codeql_build_dir/expui:/home/runner/work/EXP/EXP/_codeql_build_dir/exputil:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/runner/work/EXP/EXP/_codeql_build_dir/extern/yaml-cpp:"
         NEW_RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu/openmpi/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/h5compare")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/runner/work/EXP/EXP/_codeql_build_dir/expui/CMakeFiles/h5compare.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/viewcoefs" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/viewcoefs")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/viewcoefs"
         RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu/openmpi/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/runner/work/EXP/EXP/_codeql_build_dir/expui/viewcoefs")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/viewcoefs" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/viewcoefs")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/viewcoefs"
         OLD_RPATH "/home/runner/work/EXP/EXP/_codeql_build_dir/expui:/home/runner/work/EXP/EXP/_codeql_build_dir/exputil:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/runner/work/EXP/EXP/_codeql_build_dir/extern/yaml-cpp:"
         NEW_RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu/openmpi/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/viewcoefs")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/runner/work/EXP/EXP/_codeql_build_dir/expui/CMakeFiles/viewcoefs.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/h5power" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/h5power")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/h5power"
         RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu/openmpi/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/runner/work/EXP/EXP/_codeql_build_dir/expui/h5power")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/h5power" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/h5power")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/h5power"
         OLD_RPATH "/home/runner/work/EXP/EXP/_codeql_build_dir/expui:/home/runner/work/EXP/EXP/_codeql_build_dir/exputil:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/runner/work/EXP/EXP/_codeql_build_dir/extern/yaml-cpp:"
         NEW_RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu/openmpi/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/h5power")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/runner/work/EXP/EXP/_codeql_build_dir/expui/CMakeFiles/h5power.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/makecoefs" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/makecoefs")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/makecoefs"
         RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu/openmpi/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/runner/work/EXP/EXP/_codeql_build_dir/expui/makecoefs")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/makecoefs" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/makecoefs")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/makecoefs"
         OLD_RPATH "/home/runner/work/EXP/EXP/_codeql_build_dir/expui:/home/runner/work/EXP/EXP/_codeql_build_dir/exputil:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/runner/work/EXP/EXP/_codeql_build_dir/extern/yaml-cpp:"
         NEW_RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu/openmpi/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/makecoefs")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/runner/work/EXP/EXP/_codeql_build_dir/expui/CMakeFiles/makecoefs.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/testread" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/testread")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/testread"
         RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu/openmpi/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/runner/work/EXP/EXP/_codeql_build_dir/expui/testread")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/testread" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/testread")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/testread"
         OLD_RPATH "/home/runner/work/EXP/EXP/_codeql_build_dir/expui:/home/runner/work/EXP/EXP/_codeql_build_dir/exputil:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/runner/work/EXP/EXP/_codeql_build_dir/extern/yaml-cpp:"
         NEW_RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu/openmpi/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/testread")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/runner/work/EXP/EXP/_codeql_build_dir/expui/CMakeFiles/testread.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/testunits" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/testunits")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/testunits"
         RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu/openmpi/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/runner/work/EXP/EXP/_codeql_build_dir/expui/testunits")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/testunits" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/testunits")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/testunits"
         OLD_RPATH "/home/runner/work/EXP/EXP/_codeql_build_dir/expui:/home/runner/work/EXP/EXP/_codeql_build_dir/exputil:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/runner/work/EXP/EXP/_codeql_build_dir/extern/yaml-cpp:"
         NEW_RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu/openmpi/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/testunits")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/runner/work/EXP/EXP/_codeql_build_dir/expui/CMakeFiles/testunits.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/home/runner/work/EXP/EXP/_codeql_build_dir/expui/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
