import os
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

# Make some fake body data
#

# Call the basis to generate coefficients
#
mass = []
xpos = []
ypos = []
zpos = []

print("---- creating list data")
for i in range(0, 100):
    mass.append(0.001)
    xpos.append(random.random()*2.0 - 1.0)
    ypos.append(random.random()*2.0 - 1.0)
    zpos.append(random.random()*2.0 - 1.0)

print("---- createFromArray usings lists")
coef1 = basis.createFromArray(mass, [xpos, ypos, zpos], time=3.0)

# mass  = np.array(mass)
# data  = np.array([xpos, ypos, zpos])

# print("---- createFromArray usings numpy arrays converted from lists")
# coef2 = basis.createFromArray(mass, data, time=3.1)

# exit(0)

mass = np.ones(1000) * 1.0e6
xpos = np.random.normal(0.0, 1.0, 1000)
ypos = np.random.normal(0.0, 1.0, 1000)
zpos = np.random.normal(0.0, 1.0, 1000)

print("---- createFromArray using a list of numpy arrays")
coef3 = basis.createFromArray(mass, [xpos, ypos, zpos], time=3.2)

data  = np.array([xpos, ypos, zpos])

print("---- createFromArray usings pure numpy arrays")
coef4 = basis.createFromArray(mass, data, time=3.3)

# Add the coefficient structure coefficient structure
coefs = pyEXP.coefs.SphCoefs(True)
coefs.add(coef1)
# coefs.add(coef2)
coefs.add(coef3)
coefs.add(coef4)
print("Times:", coefs.Times())

exit(0)
