find_package (Python3 COMPONENTS Interpreter Development)

# Build and configure the installation of the Python bindings
#
pybind11_add_module(pyEXP PyWrappers.cc CoefWrappers.cc
  BasisWrappers.cc FieldWrappers.cc ParticleReaderWrappers.cc MSSAWrappers.cc)
target_link_libraries(pyEXP PUBLIC expcoef exputil ${common_LINKLIB})
target_include_directories(pyEXP PUBLIC ${common_INCLUDE})
target_compile_options(pyEXP PUBLIC ${OpenMP_CXX_FLAGS})
get_target_property(cxxflags pyEXP COMPILE_OPTIONS)
message("The project has set the following flags: ${cxxflags}")

install(TARGETS pyEXP
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION "lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages"
  ARCHIVE DESTINATION "lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages")

