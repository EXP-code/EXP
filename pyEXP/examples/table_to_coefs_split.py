import os
import time
import pyEXP
import numpy as np
import matplotlib.pyplot as plt

os.chdir('/data/Nbody/Better')

# Make the halo basis config
config="""
---
id: sphereSL
parameters :
  numr: 2000
  rmin: 0.0001
  rmax: 1.95
  Lmax: 4
  nmax: 10
  scale: 0.0667
  modelname: SLGridSph.model
...
"""

# Construct the basis instances
#
basis = pyEXP.basis.Basis.factory(config)

# Open the file and read the array.  The file was generated using
# psp2ascii.
#
bodyfile = 'comp.dark halo'
if not os.path.exists(bodyfile):
    print("You need to use psp2ascii to make an ascii input file, e.g.\n"
          "'psp2ascii -f OUT.run0.00010'")
    exit(1)

data = np.loadtxt(bodyfile, skiprows=1, usecols=(1, 2, 3, 4))
print(data.shape)


# Test partitioning
#

asize = data.shape[0]
npart = 10
bunch = int(asize/npart)

# Setup for coefficient accumulation
#
basis.initFromArray()

for i in range(npart+1):
    ibeg = bunch*i
    if ibeg<asize:
        iend = bunch*(i+1)
        if iend>asize: iend = asize
        if ibeg < iend:
            basis.addFromArray(data[ibeg:iend,0], data[ibeg:iend,1:4])
            print('beg={} end={} [{}]'.format(ibeg, iend, asize))

# Done, get the coefficient structure
#
coef = basis.makeFromArray(time=3.0)

# Make a Coef instance and add the coefficient structure for this time
#
coefs = pyEXP.coefs.Coefs.makecoefs(coef)
coefs.add(coef)

# Write it to HDF5
#
coefs.WriteH5Coefs("table_to_coefs_split");
print("Wrote coefficients")

