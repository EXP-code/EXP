#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pyEXP

# Make the halo basis config
halo_config="""
---
id: sphereSL
parameters :
  numr:  2000       # Number of radial grid points
  rmin:  0.0001     # Minimum radius (make > 0 for a divergent cusp)
  rmax:  1.95       # Maximum radius
  Lmax:  4          # Maximum spherical harmonic order
  nmax:  10         # Maximum radial basis function order
  scale: 0.0667     # Characteristic scale for coordindate mapping
  modelname: SLGridSph.model   # The model file name
  cachename: .slgrid_sph_cache # The basis function cache file name
...
"""

# Construct the basis instances
#
halo_basis = pyEXP.basis.Basis.factory(halo_config)


# Look into the cache file
#
node_sph = halo_basis.cacheInfo('.slgrid_sph_cache')

# Orthogonality test
#
if (halo_basis.orthoTest()):
    exit(0)
else:
    exit(1)
