# Is pyEXP built? If yes, run pyEXP tests.
#
if(ENABLE_PYEXP)
  # Check that pyEXP can be loaded; does not execute anything
  add_test(NAME pyexpLoadTest
    COMMAND python3 test.py
    WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}")

  # Append the pyEXP build directory to PYTHONPATH.  The property might
  # be changed depending on arch target.  Here, I assume that it's valid
  # for all arch targets.
  set_property(TEST pyexpLoadTest PROPERTY ENVIRONMENT_MODIFICATION
    "PYTHONPATH=path_list_prepend:${CMAKE_BINARY_DIR}/pyexp}")
endif()

# Is n-body EXP built? If yes, run n-body tests
#
if(ENABLE_NBODY)
  # Check that exp executes; does not test any functionality
  add_test(NAME expExecuteTest
    COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/src/exp -v)

  # Makes some spherical ICs using utils/ICs/gensph
  add_test(NAME makeICTest
    COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/utils/ICs/gensph -N 1000 -i SLGridSph.model
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/Halo)

  # Runs those ICs using exp
  add_test(NAME expNbodyTest
    COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/src/exp config.yml
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/Halo)

  # Check OUTLOG file for a sane 2T/W mean and stdv
  add_test(NAME expNbodyCheck2TW
    COMMAND COMMAND python3 check.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/Halo)

  # Remove the temporary files
  add_test(NAME removeTempFiles
    COMMAND ${CMAKE_COMMAND} -E remove
    config.run0.yml current.processor.rates.run0 new.bods new.bods.0
    OUTLOG.run0 run0.levels SLGridSph.cache.run0 test.grid
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/Halo)

  # Will removes these files if they all exist; perhaps there is a better way?
  set_tests_properties(removeTempFiles PROPERTIES DEPENDS expNbodyTest
    REQUIRED_FILES "config.run0.yml;current.processor.rates.run0;new.bods;new.bods.0;OUTLOG.run0;run0.levels;SLGridSph.cache.run0;test.grid;"
    )
endif()

# Look for some file; nothing to do with EXP but might be useful later
# on...
add_test(NAME scriptTest
  COMMAND "${CMAKE_CURRENT_LIST_DIR}/script_test.sh"
  WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}")
