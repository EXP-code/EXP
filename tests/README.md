# This directory contains unit tests for `exp` and `pyEXP`

## TL;DR

The tests are enabled by default.  You make disable the tests by including
`-DENABLE_TESTS=off` in CMake.

Run the tests from your build directory using the `make test` or
`ctest` command.

You may perform a specific test by name using the `ctest -R` option.
For example, the trivial execute test for `exp` can be performed as

```bash
ctest -R expExecuteTest
```

You can see the deatiled output from the test using the `ctest
--verbose` option.  For example,

```bash
ctest -V -R expExecuteTest
```

## Notes for developers

The `CMakeLists.txt` file is divided into three sections:
1. An N-body exp code section that is run if the `ENABLE_NBODY` flag is
   	set in CMake.  These tests, then, can run simulations and use
	stand-alone utilities from the `utils` directory structure.
2. A pyEXP section that is run if `ENABLE_PYEXP` is set in CMake.
3. A general section that is a placeholder for any tests that are not
   specific to a particular configuration. 

The `pyEXP` section needs more detailed unit tests that exercise each
of the major `pyEXP` classes.  TBD.

The `EXP` section currently includes a simple IC generation and N-body
simulation workflow, ending with a sanity check.  Additional unit
tests might check some of the additional stand-alone tools and perhaps
more of the specific N-body features such as the `Cylindrical` or
`Cube` bases.
