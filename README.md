# EXP: nbody EXPansion code

## Repo organization

| Files   |     |
| ---     | --- |
|README	  | This file |
| TODO	  | Wish list of features to add, things to fix, big and small |
| include | Include file for all common classes |
| src	  | Source for nbody code |
| exputil | The main EXP function and class library used by both pyEXP and EXP |
| coefs   | Source for standalone classes used by Python wrappers and other standalone utilities |
| pyEXP   | Source for Python wrappers |
| utils   | Older but still useful standalone C++ utilities |

## Version reporting

EXP automatically stashes its compile time, git branch, and git commit
hash when 'make' is invoked in the src directory.  You can see this
info using the -v flag, i.e. 'mpirun -np 1 exp -v' or 'exp -v'.  Note:
some MPI implementations require the MPI-aware executable to be run
using 'mpirun'.

## Compile hints

See README.build for a brief synposis.

If you are using DSMC, download and install the CHIANTI database from https://www.chiantidatabase.org

Add the path to where you install it by adding the following to your .bashrc file:

export CHIANTI_DATA=/path_to_install_location/CHIANTI/dbase

replacing the path_to_install_location as apropriate.

Also set the LD_LIBRARY_PATH by adding the following to your .bashrc file:

export LD_LIBRARY_PATH=${HOME}/lib:${LD_LIBRARY_PATH}

Log out and then back in in order for this to take effect.

