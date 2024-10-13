import os, time
import pyEXP
import random
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
  modelname: SLGridSph.model
  cachename: SLGridSph.cache.run0
...
"""

print("---- about to create basis")

# Construct the basis instances
#
basis = pyEXP.basis.Basis.factory(config)

print("---- created basis")

time.sleep(1)

# Create a coefficient structure
#
coefs = pyEXP.coefs.SphCoefs(True)

print("---- created coefficients")

time.sleep(1)

exit(0)

# Call the basis to generate coefficients
#
mass = []
xpos = []
ypos = []
zpos = []

print("---- creating list data")
for i in range(0, 100):
    mass.append(0.01)
    xpos.append(random.random()*2.0 - 1.0)
    ypos.append(random.random()*2.0 - 1.0)
    zpos.append(random.random()*2.0 - 1.0)

print("---- createFromArray usings lists")
time.sleep(1)
coef1 = basis.createFromArray(mass, [xpos, ypos, zpos], time=3.0)
time.sleep(1)
coefs.add(coef1)

print("Times:", coefs.Times())

exit(0)

print("---- creating array data from list data")

# Note: this overwrites mass and data variables.  But data is now
# passed by rference so pybind11 should have copies and not break the
# references though garbage collection on the Python side
#
mass = np.array(mass)
data = np.array([xpos, ypos, zpos])

print("---- createFromArray usings numpy arrays converted from lists")
coef2 = basis.createFromArray(mass, data, time=3.1)

coefs.add(coef2)

print("---- creating array data from list of numpy arrays")

mass = np.ones(100) * 1.0e-02
xpos = np.random.normal(0.0, 1.0, 100)
ypos = np.random.normal(0.0, 1.0, 100)
zpos = np.random.normal(0.0, 1.0, 100)

print("---- createFromArray using a list of numpy arrays")
coef3 = basis.createFromArray(mass, [xpos, ypos, zpos], time=3.2)

coefs.add(coef3)

data  = np.array([xpos, ypos, zpos])
print("---- position data shape is:", data.shape)

print("---- createFromArray usings pure numpy arrays")
coef4 = basis.createFromArray(mass, data, time=3.3)

coefs.add(coef4)

# Check the coefficient structure for the 4 added times
#
print("Times:", coefs.Times())

exit(0)
