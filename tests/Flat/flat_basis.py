#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pyEXP

# Make the flatdisk basis config
#
disk_config = """
id           : flatdisk
parameters   :
  scale      : 0.01
  Mmax       : 6
  nmax       : 12
  nmaxfid    : 40
  numx       : 64
  numy       : 32 
  numr       : 2000
  acyltbl    : 0.6
  rcylmax    : 6.0
  logr       : false
  cachename  : .biorth_cache
"""

# Construct the basis instance
#
disk_basis = pyEXP.basis.Basis.factory(disk_config)

# Orthogonality test
#
if (disk_basis.orthoTest()):
    exit(0)
else:
    exit(1)
