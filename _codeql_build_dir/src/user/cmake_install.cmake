# Install script for directory: /home/runner/work/EXP/EXP/src/user

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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libuser.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libuser.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libuser.so"
         RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/user" TYPE SHARED_LIBRARY FILES "/home/runner/work/EXP/EXP/_codeql_build_dir/src/user/libuser.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libuser.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libuser.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libuser.so"
         OLD_RPATH "/home/runner/work/EXP/EXP/_codeql_build_dir/src:/home/runner/work/EXP/EXP/_codeql_build_dir/expui:/home/runner/work/EXP/EXP/_codeql_build_dir/exputil:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/runner/work/EXP/EXP/_codeql_build_dir/extern/yaml-cpp:/usr/lib/x86_64-linux-gnu/hdf5/serial:"
         NEW_RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libuser.so")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/runner/work/EXP/EXP/_codeql_build_dir/src/user/CMakeFiles/user.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libbar.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libbar.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libbar.so"
         RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/user" TYPE SHARED_LIBRARY FILES "/home/runner/work/EXP/EXP/_codeql_build_dir/src/user/libbar.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libbar.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libbar.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libbar.so"
         OLD_RPATH "/home/runner/work/EXP/EXP/_codeql_build_dir/src:/home/runner/work/EXP/EXP/_codeql_build_dir/expui:/home/runner/work/EXP/EXP/_codeql_build_dir/exputil:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/runner/work/EXP/EXP/_codeql_build_dir/extern/yaml-cpp:/usr/lib/x86_64-linux-gnu/hdf5/serial:"
         NEW_RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libbar.so")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/runner/work/EXP/EXP/_codeql_build_dir/src/user/CMakeFiles/bar.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/liblogpot.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/liblogpot.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/liblogpot.so"
         RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/user" TYPE SHARED_LIBRARY FILES "/home/runner/work/EXP/EXP/_codeql_build_dir/src/user/liblogpot.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/liblogpot.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/liblogpot.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/liblogpot.so"
         OLD_RPATH "/home/runner/work/EXP/EXP/_codeql_build_dir/src:/home/runner/work/EXP/EXP/_codeql_build_dir/expui:/home/runner/work/EXP/EXP/_codeql_build_dir/exputil:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/runner/work/EXP/EXP/_codeql_build_dir/extern/yaml-cpp:/usr/lib/x86_64-linux-gnu/hdf5/serial:"
         NEW_RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/liblogpot.so")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/runner/work/EXP/EXP/_codeql_build_dir/src/user/CMakeFiles/logpot.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libdisk.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libdisk.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libdisk.so"
         RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/user" TYPE SHARED_LIBRARY FILES "/home/runner/work/EXP/EXP/_codeql_build_dir/src/user/libdisk.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libdisk.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libdisk.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libdisk.so"
         OLD_RPATH "/home/runner/work/EXP/EXP/_codeql_build_dir/src:/home/runner/work/EXP/EXP/_codeql_build_dir/expui:/home/runner/work/EXP/EXP/_codeql_build_dir/exputil:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/runner/work/EXP/EXP/_codeql_build_dir/extern/yaml-cpp:/usr/lib/x86_64-linux-gnu/hdf5/serial:"
         NEW_RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libdisk.so")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/runner/work/EXP/EXP/_codeql_build_dir/src/user/CMakeFiles/disk.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libhalo.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libhalo.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libhalo.so"
         RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/user" TYPE SHARED_LIBRARY FILES "/home/runner/work/EXP/EXP/_codeql_build_dir/src/user/libhalo.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libhalo.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libhalo.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libhalo.so"
         OLD_RPATH "/home/runner/work/EXP/EXP/_codeql_build_dir/src:/home/runner/work/EXP/EXP/_codeql_build_dir/expui:/home/runner/work/EXP/EXP/_codeql_build_dir/exputil:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/runner/work/EXP/EXP/_codeql_build_dir/extern/yaml-cpp:/usr/lib/x86_64-linux-gnu/hdf5/serial:"
         NEW_RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libhalo.so")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/runner/work/EXP/EXP/_codeql_build_dir/src/user/CMakeFiles/halo.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libmndisk.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libmndisk.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libmndisk.so"
         RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/user" TYPE SHARED_LIBRARY FILES "/home/runner/work/EXP/EXP/_codeql_build_dir/src/user/libmndisk.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libmndisk.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libmndisk.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libmndisk.so"
         OLD_RPATH "/home/runner/work/EXP/EXP/_codeql_build_dir/src:/home/runner/work/EXP/EXP/_codeql_build_dir/expui:/home/runner/work/EXP/EXP/_codeql_build_dir/exputil:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/runner/work/EXP/EXP/_codeql_build_dir/extern/yaml-cpp:/usr/lib/x86_64-linux-gnu/hdf5/serial:"
         NEW_RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libmndisk.so")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/runner/work/EXP/EXP/_codeql_build_dir/src/user/CMakeFiles/mndisk.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libmwgala.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libmwgala.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libmwgala.so"
         RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/user" TYPE SHARED_LIBRARY FILES "/home/runner/work/EXP/EXP/_codeql_build_dir/src/user/libmwgala.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libmwgala.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libmwgala.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libmwgala.so"
         OLD_RPATH "/home/runner/work/EXP/EXP/_codeql_build_dir/src:/home/runner/work/EXP/EXP/_codeql_build_dir/expui:/home/runner/work/EXP/EXP/_codeql_build_dir/exputil:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/runner/work/EXP/EXP/_codeql_build_dir/extern/yaml-cpp:/usr/lib/x86_64-linux-gnu/hdf5/serial:"
         NEW_RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/user/libmwgala.so")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/runner/work/EXP/EXP/_codeql_build_dir/src/user/CMakeFiles/mwgala.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/home/runner/work/EXP/EXP/_codeql_build_dir/src/user/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
