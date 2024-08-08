#!/usr/bin/env python
# coding: utf-8

# import numpy as np
import pyEXP

# Make the disk basis config
#
disk_config = """
---
id: cylinder
parameters:
  acyl: 0.01       # The scale length of the exponential disk
  hcyl: 0.001      # The scale height of the exponential disk
  lmaxfid: 20      # The maximum spherical harmonic order for the input basis
  nmaxfid: 20      # The radial order for the input spherical basis
  mmax: 6          # The maximum azimuthal order for the cylindrical basis
  nmax: 8          # The maximum radial order of the cylindrical basis
  ncylnx: 128      # The number of grid points in mapped cylindrical radius
  ncylny: 64       # The number of grid points in mapped verical scale
  ncylodd: 3       # The number of anti-symmetric radial basis functions per azimuthal order m
  rnum: 32         # The number of radial integration knots in the inner product
  pnum: 0          # The number of azimuthal integration knots (pnum: 0, assume axisymmetric target density)
  tnum: 16         # The number of colatitute integration knots
  ashift: 0.5      # Target shift length in scale lengths to create more variance
  vflag: 16        # Verbosity flag: print diagnostics to stdout for vflag>0
  logr: false      # Log scaling in cylindrical radius
  density: true    # Compute the density functions
  eof_file: .eof.cache.run0t  # The cache file name
  ignore: true
...
"""

# Construct the basis instances
#
disk_basis = pyEXP.basis.Basis.factory(disk_config)


# Look into the cache file
#
node_cyl = disk_basis.cacheInfo('.eof.cache.run0t')
