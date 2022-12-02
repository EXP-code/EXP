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

# Call the basis to generate coefficients
#
coef = basis.createFromArray(data[:,0], data[:,1:4], time=3.0)

# Print the data for a check
#
print("Time=", coef.time, " geometry=", coef.geometry)
print("Shape=", coef.data.shape)
print("Data=\n", coef.data)

# Make an HDF5 file
#
coefs = pyEXP.coefs.Coefs.makecoefs(coef)
coefs.add(coef)
coefs.WriteH5Coefs("table_to_coefs");
print("Wrote coefficients")

#
# DONE
