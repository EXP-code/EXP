# EXP: nbody EXPansion code

## Repo organization

| Files   | Description |
| ---     | ---         |
| README  | This file |
| TODO	  | Wish list of features to add, things to fix, big and small |
| include | Include file for all common classes |
| src	  | Source for nbody code |
| exputil | The main EXP function and class library used by both pyEXP and EXP |
| coefs   | Source for standalone classes used by Python wrappers and other standalone utilities |
| pyEXP   | Source for Python wrappers |
| utils   | Older but still useful standalone C++ utilities |

## Version reporting

EXP automatically stashes its compile time, git branch, and git commit
hash when `make` is invoked in the src directory.  You can see this
info using the -v flag, i.e. `mpirun -np 1 exp -v` or `exp -v`.  Note:
some MPI implementations require the MPI-aware executable to be run
using 'mpirun'.  Some recent HPC systems using `slurm` require the use
of `srun` instead of `mpirun` to better administer and schedule your
resource request.

## Compile hints

See README.build for a brief synposis.

A few quick additional notes. By default, both the n-body code and the
Python bindings, the pyEXP interface, will be compiled by default.
For those of you that only want pyEXP, add `-DENABLE_NBODY=OFF` to
your `cmake` invocation or toggle `ENABLE_NBODY` using `ccmake` or
your favorite gui configurator.

## Documentation

Currently, EXP is extensively documented using `doxygen`.  The
documentation is in `doc/html`.  You will need to enable the build
using the `cmake` flag `-DBUILD_DOCS=ON`.  A permanent online presense
is in the works.

## Companion repositories

We are developing two repositories of examples and tutorials:

| Repo   | Description |
| ---    | ---         |
| EXP-examples | Each subdirectory contains a full set of body files and configurations to run EXP with with model galaxy |
| pyEXP-examples | Tutorials and example workflows for a variety of envisioned use cases |

Both of these are available from the origin as EXP.

## pyEXP

Provides a collection of EXP tools for processing and analyzing
simulation data using BFE techniques and MSSA.

The main documentation is the many docstrings embedded in the code and
a set of examples provided in the auxiliary pyEXP-examples repository.
An online reference guide is in the works. We
recommend beginning with the example Python scripts and IPython
notebooks and adapting them to your own needs. You can explore
the available classes and member functions using the usual Python
`help` function.  The classes are organized into six submodules
that are described briefly below.  Run `help(pyEXP.xxxx)` for each
of the submodules below for more detailed usage info...

#### The pyEXP submodules

| Submodule | Description |
| ---       | ---         |
| read      | Read particle snapshots of various types.  Currently EXP, Gadget, Tipsy, and Bonzai types are supported. |
| basis     | Create and apply specific biorthogonal bases to generate coefficients from particle data and evaluate potential, density, and force fields. |
| coefs     | Classes for reading, passing, writing, converting, and querying coefficient sets. |
| field     | Create two- and three-dimension rectangular grids of fields for visualization. |
| mssa      | Tools to apply Multivariate Singular Spectrum Analysis (MSSA) to the coefficients computed using the 'basis' classes. |
| util      | Miscellaneous tools that support the others.  Currently this include centering algorithms.  While EXP has native methods for doing this, others will need to supply an estimated center. |

#### pyEXP example workflow

To provide some context, suppose you want to read some snapshots, make some coefficients, and then analyze them with MSSA. The class construction would go something like this:

1. Create a reader instance for your simulation, call it 'reader'.
2. Create a basis designed to represent the a particular particle type.  Star particles, perhaps, so let's call it 'disk'.
3. We then pass 'reader' to the createCoefficients member of 'disk' to get coefficients for your snapshots, called 'coefs'
4. We might then want to explore dynamical patterns in these coefficients by passing 'coefs' to 'expMSSA'.  'expMSSA' will return principal signals as an updated coefficient object, that we call 'newcoefs'
5. 'newcoefs' and 'disk' can then be passed to the FieldGenerator to provide density, potential, force fields, etc. for the each principal signal

This is only one example of many possible uses.  There are many
variants to this work flow, of course, and I expect that you will
invent some interesting ones.

