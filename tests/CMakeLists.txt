include(CTest)

# Is EXP configured for pyEXP? If yes, run pyEXP tests...
#
if(ENABLE_PYEXP)
  # Check that pyEXP can be loaded; does not execute anything
  add_test(NAME pyexpLoadTest
    COMMAND ${CMAKE_COMMAND} -E env
    PYTHONPATH=${CMAKE_BINARY_DIR}/pyEXP:$ENV{PYTHONPATH}
    ${PYTHON_EXECUTABLE} import.py
    WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}")
  
  # Make a spherical basis using pyEXP
  add_test(NAME pyexpSphBasisTest
    COMMAND ${CMAKE_COMMAND} -E env
    PYTHONPATH=${CMAKE_BINARY_DIR}/pyEXP:$ENV{PYTHONPATH}
    ${PYTHON_EXECUTABLE} sph_basis.py
    WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/Halo")

  # Make a cylindrical basis using pyEXP
  add_test(NAME pyexpCylBasisTest
    COMMAND ${CMAKE_COMMAND} -E env
    PYTHONPATH=${CMAKE_BINARY_DIR}/pyEXP:$ENV{PYTHONPATH}
    ${PYTHON_EXECUTABLE} cyl_basis.py
    WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/Disk")

  # Remove cache files
  add_test(NAME removeSphCache
    COMMAND ${CMAKE_COMMAND} -E remove .slgrid_sph_cache
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/Halo)

  add_test(NAME removeCylCache
    COMMAND ${CMAKE_COMMAND} -E remove .eof.cache.run0t test_adddisk_sl.0
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/Disk)

  set_tests_properties(removeSphCache PROPERTIES DEPENDS pyexpSphBasisTest
    REQUIRED_FILES ".slgrid_sph_cache")

  set_tests_properties(removeCylCache PROPERTIES DEPENDS pyexpCylBasisTest
    REQUIRED_FILES ".eof.cache.run0t;test_adddisk_sl.0")

  # Other tests for pyEXP go here ...
  
  # Set labels for pyEXP tests
  set_tests_properties(pyexpLoadTest pyexpSphBasisTest PROPERTIES LABELS "quick")
  set_tests_properties(pyexpCylBasisTest PROPERTIES LABELS "long")

endif()

# Is EXP configured for the N-body code? If yes, run N-body tests...
#
if(ENABLE_NBODY)
  # Check that exp executes; does not test any functionality
  add_test(NAME expExecuteTest
    COMMAND ${EXP_MPI_LAUNCH} ${CMAKE_BINARY_DIR}/src/exp -v)

  # Makes some spherical ICs using utils/ICs/gensph
  add_test(NAME makeICTest
    COMMAND ${EXP_MPI_LAUNCH} ${CMAKE_BINARY_DIR}/utils/ICs/gensph -N 1000 -i SLGridSph.model
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/Halo)

  # Runs those ICs using exp
  add_test(NAME expNbodyTest
    COMMAND ${EXP_MPI_LAUNCH} ${CMAKE_BINARY_DIR}/src/exp config.yml
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/Halo)

  set_tests_properties(expNbodyTest PROPERTIES DEPENDS makeICTest)

  # Check OUTLOG file for a sane 2T/W mean and stdv
  add_test(NAME expNbodyCheck2TW
    COMMAND ${CMAKE_COMMAND} -E env
    PYTHONPATH=${CMAKE_BINARY_DIR}/pyEXP:$ENV{PYTHONPATH}
    ${PYTHON_EXECUTABLE} check.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/Halo)

  set_tests_properties(expNbodyCheck2TW PROPERTIES DEPENDS makeNbodyTest)

  # This adds a coefficient read test using pyEXP only if
  # makeNbodyTest is run and pyEXP has been built
  if(ENABLE_PYEXP)
     # Read coefficient file with pyEXP
     add_test(NAME pyEXPCoefReadTest
              COMMAND ${CMAKE_COMMAND} -E env
              PYTHONPATH=${CMAKE_BINARY_DIR}/pyEXP:$ENV{PYTHONPATH}
              ${PYTHON_EXECUTABLE} readCoefs.py
              WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/Halo")
     set_tests_properties(pyEXPCoefReadTest PROPERTIES DEPENDS makeNbodyTest)
  endif()

  # A separate test to remove the generated files if they all exist;
  # perhaps there is a better way?
  add_test(NAME removeTempFiles
    COMMAND ${CMAKE_COMMAND} -E remove
    config.run0.yml current.processor.rates.run0 new.bods new.bods.0
    OUTLOG.run0 run0.levels SLGridSph.cache.run0 test.grid
    'outcoef.dark halo.run0'
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/Halo)

  # Remove the temporary files
  set_tests_properties(removeTempFiles PROPERTIES DEPENDS expNbodyCheck2TW
    REQUIRED_FILES "config.run0.yml;current.processor.rates.run0;new.bods;new.bods.0;OUTLOG.run0;run0.levels;SLGridSph.cache.run0;test.grid;"
    )

  # Makes some cube ICs using utils/ICs/cubeics
  add_test(NAME makeCubeICTest
    COMMAND ${EXP_MPI_LAUNCH} ${CMAKE_BINARY_DIR}/utils/ICs/cubeics -N 4000 -z -d 2,2,2
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/Cube)

  # Runs those ICs using exp
  add_test(NAME expCubeTest
    COMMAND ${EXP_MPI_LAUNCH} ${CMAKE_BINARY_DIR}/src/exp
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/Cube)

  set_tests_properties(expCubeTest PROPERTIES DEPENDS makeCubeICTest)

  # Check OUTLOG file for mean position
  add_test(NAME expCubeCheckPos
    COMMAND ${CMAKE_COMMAND} -E env
    PYTHONPATH=${CMAKE_BINARY_DIR}/pyEXP:$ENV{PYTHONPATH}
    ${PYTHON_EXECUTABLE} check.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/Cube)

  set_tests_properties(expCubeCheckPos PROPERTIES DEPENDS expCubeTest)

  # A separate test to remove the generated files if they all exist
  add_test(NAME removeCubeFiles
    COMMAND ${CMAKE_COMMAND} -E remove
    config.runS.yml current.processor.rates.runS cube.bods
    OUTLOG.runS runS.levels
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/Cube)

  # Remove the temporary files
  set_tests_properties(removeCubeFiles PROPERTIES DEPENDS expCubeCheckPos
    REQUIRED_FILES "config.runS.yml;current.processor.rates.runS;cube.bods;OUTLOG.runS;runS.levels;")

  # Set labels for pyEXP tests
  set_tests_properties(expExecuteTest PROPERTIES LABELS "quick")
  set_tests_properties(makeICTest expNbodyTest expNbodyCheck2TW
  pyEXPCoefReadTest removeTempFiles expCubeTest removeCubeFiles
  PROPERTIES LABELS "long")

endif()

# Tests that should work for any configuration.  Nothing so far.  This
# placeholder look for some file; nothing to do with EXP but might be
# useful later on...
#
add_test(NAME scriptTest
  COMMAND "${CMAKE_CURRENT_LIST_DIR}/script_test.sh"
  WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}")

set_tests_properties(scriptTest PROPERTIES LABELS "quick")
