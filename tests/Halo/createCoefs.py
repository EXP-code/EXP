import os
import pyEXP
import numpy as np

# Make the halo basis config
config="""
---
id : sphereSL
parameters :
  numr: 4000
  rmin: 0.0001
  rmax: 1.95
  Lmax: 2
  nmax: 10
  rmapping : 0.0667
  self_consistent: true
  modelname: SLGridSph.model
  cachename: SLGridSph.cache.run0
...
"""

# Construct the basis instances
#
basis = pyEXP.basis.Basis.factory(config)

# Open the file and read the array.  The file was generated using
# psp2ascii.
#
bodyfile = 'new.bods'
if not os.path.exists(bodyfile):
    print('Body file <{}> does not exist'.format(bodyfile))
    exit(1)

data = np.loadtxt(bodyfile, skiprows=1, usecols=(1, 2, 3, 4))
print(data.shape)

# Call the basis to generate coefficients
#
coef = basis.createFromArray(data[:,0], data[:,1:4], time=3.0)

# Add the coefficient structure coefficient structure
coefs = pyEXP.coefs.SphCoefs(True)
coefs.add(coef)
print("Times:", coefs.Times())

exit(0)
