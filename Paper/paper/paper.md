---
title: 'EXP: a Python/C++ package for basis
	function expansion methods in galactic dynamics'
tags:
  - C++
  - Python
  - astronomy
  - dynamics
  - galactic dynamics
  - Milky Way
authors:
  - name: Michael S. Petersen
    orcid: 0000-0003-1517-3935
    affiliation: 1
  - name: Martin D. Weinberg
    orcid: 0000-0003-2660-2889
    affiliation: 2
affiliations:
 - name: University of Edinburgh, UK
   index: 1
 - name: University of Massachusetts/Amherst, USA
   index: 2
date: 01 June 2024
bibliography: paper.bib

---

# Summary

Galaxies are ensembles of dark and baryonic matter confined by their
mutual gravitational attraction. The dynamics of galactic evolution
have been studied with a combination of numerical simulations and
analytical methods [@Binney:2008] for decades.  Recently, data
describing the positions and motions of billions of stars from the
_Gaia satellite_ [@gaia; @gaia_DR2_disk] depict a Milky Way (our home
galaxy) much farther from equilibrium than we imagined and beyond the
range of many analytic models.  Simulations that allow collections of
bodies to evolve under their mutual gravity are capable of reproducing
such complexities but robust links to fundamental theoretical
explanations are still missing.

Basis Function Expansions (BFE) represent fields as a linear
combination of orthogonal functions. BFEs are particularly well-suited
for studies of perturbations from equilibrium, such as the evolution
of a galaxy.  For any galaxy simulation, a biorthogonal BFE can fully
represent the density, potential and forces by time series of
coefficients.  The coefficients have physical meaning: they represent
the gravitational potential energy in a given function.  The variation
the function coefficients in time encodes the dynamical evolution.
The representation of simulation data by BFE results in huge
compression of the information in the dynamical fields; for example,
1.5 TB of phase space data enumerating the positions and velocities of
millions of particles becomes 200 MB of coefficient data!

For optimal representation, the lowest-order basis function should be
similar to the mean or equilibrium profile of the galaxy. This allows
for use of the fewest number of terms in the expansion. For example,
the often-used basis set from @Hernquist:92 matches the @Hernquist:90
dark matter halo profile, yet this basis set is inappropriate for
representing the cosmologically-motivated @NFW profile.  The `EXP`
software package implements the adaptive empirical orthogonal function
(EOF; see \autoref{fig:examplecylinder}) basis strategy originally
described in @Weinberg:99 that matches any physical system close to an
equilibrium model. The package includes both a high performance N-body
simulation toolkit with computational effort scaling linearly with N
[@Petersen:22], and a user-friendly Python interface called `pyEXP`
that enables BFE and time-series analysis of any N-body simulation
dataset.

# Statement of need

The need for methodology that seamlessly connects theoretical
descriptions of dynamics, N-body simulations, and compact descriptions
of observed data gave rise to `EXP`. This package provides recent
developments from applied mathematics and numerical computation to
represent complete series of _Basis Function Expansions_ that describe
the variation of _any_ field in space.  In the context of galactic
dynamics, these fields may be density, potential, force, velocity
fields or any intrinsic field produced by simulations such as
chemistry data.  By combining the coefficient information through time
using multichannel singular spectral analysis (mSSA; @SSA), a
non-parametric spectral technique, `EXP` can deepen our understanding
by discovering the dynamics of galaxy evolution directly from
simulated, and by analogy, observed data.

`EXP` decomposes a galaxy into multiple bases for a variety of scales
and geometries and is thus able to represent arbitrarily complex
simulation with many components (e.g., disk, bulge, dark matter halo,
satellites).  `EXP` is able to efficiently summarize the degree and
nature of asymmetries through coefficient amplitudes tracked through
time and provide details at multiple scales. The
amplitudes themselves enable ex-post-facto dynamical discovery.  
`EXP` is a collection of object-oriented C++ libraries with an
associated modular N-body code and a suite of stand-alone analysis
applications.

`pyEXP` provides a full Python interface to the `EXP` libraries,
implemented with `pybind11` [@pybind11], which provides full
interoperability with major astronomical packages including Astropy
[@astropy] and `Gala` [@gala]. Example workflows based on previously
published work are available and distributed as accompanying [examples
and tutorials](https://github.com/EXP-code/pyEXP-examples).  The
examples and tutorials flatten the learning curve for
adopting BFE tools to generate and analyze the significance of
coefficients and discover dynamical relationships using time series
analysis such as mSSA.  We provide a [full online
manual](https://exp-docs.readthedocs.io) hosted by ReadTheDocs.

The software package brings published -- but difficult to implement --
applied-math technologies into the astronomical mainstream.  `EXP` and
the associated Python interface `pyEXP` accomplish this by providing
tools integrated with the Python ecosystem, and in particular are
well-suited for interactive Python [@iPython] use through (e.g.)
Jupyter notebooks [@jupyter]. `EXP` serves as the
scaffolding for new imaginative applications in galactic dynamics,
providing a common dynamical language for simulations and analytic
theory.

# Features and workflow

The core `EXP` library is built around methods to build the best basis
function expansion for an arbitrary data set in galactic dynamics.
The table below lists some of the available basis functions. All
computed bases and resulting coefficient data are stored in HDF5
[@hdf5] format.

| Name        | Description |
| ----------- | -------------        |
| sphereSL    | Sturm-Liouville basis function solutions to Poisson's equation for any arbitrary input spherical density |
| bessel      | Basis constructed from eigenfunctions of the spherical Laplacian |
| cylinder    | EOF solutions tabulated on the meridional plane for distributions with cylindrical geometries |
| flatdisk    | EOF basis solutions for the three-dimensional gravitational field of a razor-thin disk |
| cube        | Trigonometric basis solution for expansions in a cube with boundary criteria |
| field       | General-purpose EOF solution for scalar profiles |
| velocity    | EOF solution for velocity flow coefficients |

![Example cylinder basis functions, where the color encodes the amplitude of the function, for an exponential disk with a scalelength of 3 and a scaleheight of 0.3 in arbitrary units. We select three functions at low, medium, and higher order (corresponding to the number of nodes). The color scale has been normalised such that the largest amplitude is unity in each panel. \label{fig:examplecylinder}](examplefunctions.png)


## N-body simulation

Our design includes a wide choice of run-time summary diagnostics,
phase-space output formats, dynamically loadable user libraries, and
easy extensibility. Stand-alone routines include the EOF and mSSA 
methods described above, and the modular software architecture of 
EXP enables users to easily build and maintain extensions. The `EXP` 
code base is described in published papers [@Petersen:22; @Weinberg:23]
and has been used, enhanced, and rigorously tested for nearly two 
decades.


The design and implementation of the N-body tools allows for execution
on a wide variety of hardware, from personal laptops to high
performance computing centers, with communication between processes
handled by MPI [@mpi41] and GPU implementations in CUDA [@cuda]. Owing
to the linear scaling of computational effort with N and the novel
GPU implementation, the N-body
methods in `EXP` deliver performance in collisionless N-body simulations
previously only accessible with large dedicated CPU clusters.

The flexible N-body software design allows users to write their own
modules for on-the-fly execution during N-body integration.  Modules
enable powerful and intricate dynamical experiments in N-body
simulations, further reducing the gap between numerical simulations
and analytic dynamics. The package ships with several examples,
including imposed external potentials, as well as a basic example that
can be extended by users.

## Using pyEXP to represent simulations

`pyEXP` provides an interface to many of the classes in the `EXP` C++
library, allowing for both the generation of all bases listed in the
table above as well as coefficients for an input data set.  Each of
these tools are Python classes that accept `numpy` [@numpy] arrays for
immediate interoperability with `matplotlib` [@matplotlib] and
Astropy.  We include a verified set of stand-alone routines that read
phase-space files from many major cosmological tree codes and produce
BFE-based analyses.  The code suite includes adapters for reading and
writing phase space for many of the widely used cosmology codes, with
a base class for developing new ones.  There are multiple ways to use
the versatile and modular tools in `pyEXP`, and we anticipate
pipelines that we have not yet imagined.


## Using pyEXP to analyze time series

The `EXP` library includes multiple time series analysis tools,
documented in the manual. Here, we briefly highlight one technique
that we have already used in published work: mSSA [@Weinberg:21;
@Johnson:23].  Beginning with coefficient series from the previous
tools, mSSA summarizes signals _in time_ that describes dynamically
correlated responses and patterns.  Essentially, this is BFE in time
and space.  These temporal and spatial patterns allow users to better
identify dynamical mechanisms and enable intercomparisons and
filtering for features in simulation suites; e.g. computing the
fraction galaxies with grand design structure or hosting
bars. Random-matrix techniques for singular-value decomposition ensure
that analyses of large data sets is possible. All mSSA decompositions
are saved in HDF5 format for reuse.

# Acknowledgements

We acknowledge the ongoing support of the [B-BFE
collaboration](https://b-bfe.org).  We also acknowledge the support of
the Center for Computational Astrophysics (CCA).  The CCA is part of
the Flatiron Institute, funded by the Simons Foundation.  We thank
Robert Blackwell for invaluable help with HPC best practices.


# References

