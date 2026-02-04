# Install script for directory: /home/runner/work/EXP/EXP/utils/SL

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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/slcheck" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/slcheck")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/slcheck"
         RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/runner/work/EXP/EXP/_codeql_build_dir/utils/SL/slcheck")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/slcheck" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/slcheck")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/slcheck"
         OLD_RPATH "/home/runner/work/EXP/EXP/_codeql_build_dir/exputil:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/runner/work/EXP/EXP/_codeql_build_dir/extern/yaml-cpp:/usr/lib/x86_64-linux-gnu/hdf5/serial:"
         NEW_RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/slcheck")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/runner/work/EXP/EXP/_codeql_build_dir/utils/SL/CMakeFiles/slcheck.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/slshift" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/slshift")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/slshift"
         RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/runner/work/EXP/EXP/_codeql_build_dir/utils/SL/slshift")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/slshift" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/slshift")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/slshift"
         OLD_RPATH "/home/runner/work/EXP/EXP/_codeql_build_dir/exputil:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/runner/work/EXP/EXP/_codeql_build_dir/extern/yaml-cpp:/usr/lib/x86_64-linux-gnu/hdf5/serial:"
         NEW_RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/slshift")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/runner/work/EXP/EXP/_codeql_build_dir/utils/SL/CMakeFiles/slshift.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/orthochk" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/orthochk")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/orthochk"
         RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/runner/work/EXP/EXP/_codeql_build_dir/utils/SL/orthochk")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/orthochk" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/orthochk")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/orthochk"
         OLD_RPATH "/home/runner/work/EXP/EXP/_codeql_build_dir/exputil:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/runner/work/EXP/EXP/_codeql_build_dir/extern/yaml-cpp:/usr/lib/x86_64-linux-gnu/hdf5/serial:"
         NEW_RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/orthochk")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/runner/work/EXP/EXP/_codeql_build_dir/utils/SL/CMakeFiles/orthochk.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/diskpot" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/diskpot")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/diskpot"
         RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/runner/work/EXP/EXP/_codeql_build_dir/utils/SL/diskpot")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/diskpot" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/diskpot")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/diskpot"
         OLD_RPATH "/home/runner/work/EXP/EXP/_codeql_build_dir/exputil:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/runner/work/EXP/EXP/_codeql_build_dir/extern/yaml-cpp:/usr/lib/x86_64-linux-gnu/hdf5/serial:"
         NEW_RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/diskpot")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/runner/work/EXP/EXP/_codeql_build_dir/utils/SL/CMakeFiles/diskpot.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/qtest" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/qtest")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/qtest"
         RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/runner/work/EXP/EXP/_codeql_build_dir/utils/SL/qtest")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/qtest" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/qtest")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/qtest"
         OLD_RPATH "/home/runner/work/EXP/EXP/_codeql_build_dir/exputil:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/runner/work/EXP/EXP/_codeql_build_dir/extern/yaml-cpp:/usr/lib/x86_64-linux-gnu/hdf5/serial:"
         NEW_RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/qtest")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/runner/work/EXP/EXP/_codeql_build_dir/utils/SL/CMakeFiles/qtest.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/eoftest" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/eoftest")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/eoftest"
         RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/runner/work/EXP/EXP/_codeql_build_dir/utils/SL/eoftest")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/eoftest" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/eoftest")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/eoftest"
         OLD_RPATH "/home/runner/work/EXP/EXP/_codeql_build_dir/exputil:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/runner/work/EXP/EXP/_codeql_build_dir/extern/yaml-cpp:/usr/lib/x86_64-linux-gnu/hdf5/serial:"
         NEW_RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/eoftest")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/runner/work/EXP/EXP/_codeql_build_dir/utils/SL/CMakeFiles/eoftest.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/oftest" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/oftest")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/oftest"
         RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/runner/work/EXP/EXP/_codeql_build_dir/utils/SL/oftest")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/oftest" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/oftest")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/oftest"
         OLD_RPATH "/home/runner/work/EXP/EXP/_codeql_build_dir/exputil:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/runner/work/EXP/EXP/_codeql_build_dir/extern/yaml-cpp:/usr/lib/x86_64-linux-gnu/hdf5/serial:"
         NEW_RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/oftest")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/runner/work/EXP/EXP/_codeql_build_dir/utils/SL/CMakeFiles/oftest.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/slabchk" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/slabchk")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/slabchk"
         RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/runner/work/EXP/EXP/_codeql_build_dir/utils/SL/slabchk")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/slabchk" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/slabchk")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/slabchk"
         OLD_RPATH "/home/runner/work/EXP/EXP/_codeql_build_dir/exputil:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/runner/work/EXP/EXP/_codeql_build_dir/extern/yaml-cpp:/usr/lib/x86_64-linux-gnu/hdf5/serial:"
         NEW_RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/slabchk")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/runner/work/EXP/EXP/_codeql_build_dir/utils/SL/CMakeFiles/slabchk.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/expontst" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/expontst")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/expontst"
         RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/runner/work/EXP/EXP/_codeql_build_dir/utils/SL/expontst")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/expontst" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/expontst")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/expontst"
         OLD_RPATH "/home/runner/work/EXP/EXP/_codeql_build_dir/exputil:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/runner/work/EXP/EXP/_codeql_build_dir/extern/yaml-cpp:/usr/lib/x86_64-linux-gnu/hdf5/serial:"
         NEW_RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/expontst")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/runner/work/EXP/EXP/_codeql_build_dir/utils/SL/CMakeFiles/expontst.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/home/runner/work/EXP/EXP/_codeql_build_dir/utils/SL/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
