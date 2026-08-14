include(CTest)

set(CTEST_OUTPUT_ON_FAILURE ON)

# HDF5 particle IC smoke checks use h5py for concise schema inspection.  Keep
# them optional so a Python installation without h5py does not disable CTest.
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import h5py, numpy"
    RESULT_VARIABLE EXP_H5PY_STATUS
    OUTPUT_QUIET ERROR_QUIET)
endif()

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
    COMMAND ${CMAKE_COMMAND} -E remove .eof.cache.run0t
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/Disk)

  set_tests_properties(removeSphCache PROPERTIES DEPENDS pyexpSphBasisTest
    REQUIRED_FILES ".slgrid_sph_cache")

  set_tests_properties(removeCylCache PROPERTIES DEPENDS pyexpCylBasisTest
    REQUIRED_FILES ".eof.cache.run0t") 

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

  if(ENABLE_UTILS)
    # Makes some spherical ICs using utils/ICs/gensph
    add_test(NAME makeICTest
      COMMAND ${EXP_MPI_LAUNCH} ${CMAKE_BINARY_DIR}/utils/ICs/gensph -N 10000 -i SLGridSph.model
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

    set_tests_properties(expNbodyCheck2TW PROPERTIES DEPENDS expNbodyTest LABELS "long")

    # This adds a coefficient read test using pyEXP only if
    # expNbodyTest is run and pyEXP has been built
    if(ENABLE_PYEXP)
       # Read coefficient file with pyEXP
       add_test(NAME pyEXPCoefReadTest
                COMMAND ${CMAKE_COMMAND} -E env
                PYTHONPATH=${CMAKE_BINARY_DIR}/pyEXP:$ENV{PYTHONPATH}
                ${PYTHON_EXECUTABLE} readCoefs.py
                WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/Halo")
       set_tests_properties(pyEXPCoefReadTest PROPERTIES DEPENDS expNbodyTest LABELS "long")

       add_test(NAME pyEXPCoefMatrixTest
                COMMAND ${CMAKE_COMMAND} -E env
                PYTHONPATH=${CMAKE_BINARY_DIR}/pyEXP:$ENV{PYTHONPATH}
                ${PYTHON_EXECUTABLE} changeCoefs.py
                WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/Halo")
       set_tests_properties(pyEXPCoefMatrixTest PROPERTIES DEPENDS expNbodyTest LABELS "long")

       add_test(NAME pyEXPCoefCreateTest
                COMMAND ${CMAKE_COMMAND} -E env
                PYTHONPATH=${CMAKE_BINARY_DIR}/pyEXP:$ENV{PYTHONPATH}
                ${PYTHON_EXECUTABLE} createCoefs.py
                WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/Halo")
       set_tests_properties(pyEXPCoefCreateTest PROPERTIES LABELS "quick")
    endif()

    # A separate test to remove the generated files if they all exist;
    # perhaps there is a better way?
    add_test(NAME removeTempFiles
      COMMAND ${CMAKE_COMMAND} -E remove
      config.run0.yml current.processor.rates.run0 new.bods
      OUTLOG.run0 run0.levels SLGridSph.cache.run0 test.grid
      outcoef.halo.run0 SLGridSph.cache.run0 coefcovar.halo.test_covar.h5
      WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/Halo)

    # Remove the temporary files
    set_tests_properties(removeTempFiles PROPERTIES DEPENDS expNbodyCheck2TW
      REQUIRED_FILES "config.run0.yml;current.processor.rates.run0;new.bods;run0.levels;SLGridSph.cache.run0;test.grid;"
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

    # These are basic sanity checks for writing and reading EXP HDF5 IC files
    if(Python3_Interpreter_FOUND AND EXP_H5PY_STATUS EQUAL 0 AND ENABLE_UTILS)

      # Create the banner / message test
      add_test(NAME "--- Description ---" 
        COMMAND ${CMAKE_COMMAND} -E echo "The following tests check for consistent HDF5 IC generation")

      # Check spherical body file generation
      add_test(NAME makeSphHDF5ICTest
        COMMAND ${EXP_MPI_LAUNCH} ${CMAKE_BINARY_DIR}/utils/ICs/gensph --hdf5 -N 16 -i SLGridSph.model --SEED 17 --zerovel -o hdf5-sph.bods
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/Halo)

      add_test(NAME checkSphHDF5ICTest
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_LIST_DIR}/check_hdf5_particles.py hdf5-sph.bods.h5 --count 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/Halo)
      set_tests_properties(checkSphHDF5ICTest PROPERTIES DEPENDS makeSphHDF5ICTest)

      add_test(NAME removeSphHDF5ICFiles
        COMMAND ${CMAKE_COMMAND} -E remove hdf5-sph.bods.h5
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/Halo)
      set_tests_properties(removeSphHDF5ICFiles PROPERTIES DEPENDS checkSphHDF5ICTest
        REQUIRED_FILES "hdf5-sph.bods.h5")
      set_tests_properties(makeSphHDF5ICTest checkSphHDF5ICTest removeSphHDF5ICFiles
        PROPERTIES LABELS "quick")

      # Check slab body file generation
      add_test(NAME makeSlabHDF5ICTest
        COMMAND ${CMAKE_BINARY_DIR}/utils/ICs/slabics --hdf5 --number 16 --seed 17 --cube --outfile hdf5-slab.bods
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/Slab)

      add_test(NAME checkSlabHDF5ICTest
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_LIST_DIR}/check_hdf5_particles.py hdf5-slab.bods.h5 --count 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/Slab)
      set_tests_properties(checkSlabHDF5ICTest PROPERTIES DEPENDS makeSlabHDF5ICTest)

      add_test(NAME removeSlabHDF5ICFiles
        COMMAND ${CMAKE_COMMAND} -E remove hdf5-slab.bods.h5
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/Slab)
      set_tests_properties(removeSlabHDF5ICFiles PROPERTIES DEPENDS checkSlabHDF5ICTest
        REQUIRED_FILES "hdf5-slab.bods.h5")
      set_tests_properties(makeSlabHDF5ICTest checkSlabHDF5ICTest removeSlabHDF5ICFiles
        PROPERTIES LABELS "quick")

      # Check cube body file generation
      add_test(NAME makeCubeHDF5ICTest
        COMMAND ${CMAKE_BINARY_DIR}/utils/ICs/cubeics --hdf5 --number 16 --seed 17 --cube --zerovel --file hdf5-cube.bods
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/Cube)

      add_test(NAME checkCubeHDF5ICTest
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_LIST_DIR}/check_hdf5_particles.py
        hdf5-cube.bods.h5 --count 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/Cube)
      set_tests_properties(checkCubeHDF5ICTest PROPERTIES DEPENDS makeCubeHDF5ICTest)

      add_test(NAME removeCubeHDF5ICFiles
        COMMAND ${CMAKE_COMMAND} -E remove hdf5-cube.bods.h5
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/Cube)
      set_tests_properties(removeCubeHDF5ICFiles PROPERTIES DEPENDS checkCubeHDF5ICTest
        REQUIRED_FILES "hdf5-cube.bods.h5")
      set_tests_properties(makeCubeHDF5ICTest checkCubeHDF5ICTest removeCubeHDF5ICFiles
        PROPERTIES LABELS "quick")
    endif()

    # Set labels for pyEXP tests
    set_tests_properties(expExecuteTest PROPERTIES LABELS "quick")
    set_tests_properties(makeICTest expNbodyTest expNbodyCheck2TW
    removeTempFiles makeCubeICTest expCubeTest removeCubeFiles
    PROPERTIES LABELS "long")
  endif()

endif()

# Tests that should work for any configuration.  Nothing so far.  This
# placeholder look for some file; nothing to do with EXP but might be
# useful later on...
#
add_test(NAME scriptTest
  COMMAND "${CMAKE_CURRENT_LIST_DIR}/script_test.sh"
  WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}")

set_tests_properties(scriptTest PROPERTIES LABELS "quick")
