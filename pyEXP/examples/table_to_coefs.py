import os
import time
import pyEXP
import numpy as np
import matplotlib.pyplot as plt

os.chdir('/home/weinberg/Nbody/Better')

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
  rs: 0.0667
  modelname: SLGridSph.model
...
"""

# Construct the basis instances
#
basis = pyEXP.basis.Basis.factory(config)

# Open the file and read the array.  The file was generated using
# psp2ascii.
#
data = np.loadtxt('comp.dark halo', skiprows=1, usecols=(1, 2, 3, 4))

# Call the basis to generate coefficients
#
coef = basis.createFromArray(data[:,0], data[:,1:3], time=3.0)

# Print the data for a check
#
print("Time=", coef.time, " geometry=", coef.geometry)
print("Shape=", coef.data.shape)
print("Data=\n", coef.data)
