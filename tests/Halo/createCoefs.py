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

# Create a coefficient structure
#
coefs = pyEXP.coefs.SphCoefs(True)

print("---- created coefficients")

# Call the basis to generate coefficients
#
mass1 = []
xpos1 = []
ypos1 = []
zpos1 = []

print("---- creating list data")
for i in range(0, 100):
    mass1.append(0.001)
    xpos1.append(random.random()*2.0 - 1.0)
    ypos1.append(random.random()*2.0 - 1.0)
    zpos1.append(random.random()*2.0 - 1.0)

print("---- createFromArray usings lists")
coef1 = basis.createFromArray(mass1, [xpos1, ypos1, zpos1], time=3.0)

coefs.add(coef1)

print("---- creating array data from list data")

mass2 = np.array(mass1)
data2 = np.array([xpos1, ypos1, zpos1])

print("---- createFromArray usings numpy arrays converted from lists")
coef2 = basis.createFromArray(mass2, data2, time=3.1)

coefs.add(coef2)

print("---- creating array data from list of numpy arrays")

mass3 = np.ones(1000) * 1.0e6
xpos3 = np.random.normal(0.0, 1.0, 1000)
ypos3 = np.random.normal(0.0, 1.0, 1000)
zpos3 = np.random.normal(0.0, 1.0, 1000)

print("---- createFromArray using a list of numpy arrays")
coef3 = basis.createFromArray(mass3, [xpos3, ypos3, zpos3], time=3.2)

coefs.add(coef3)

data4  = np.array([xpos3, ypos3, zpos3])
print("---- position data shape is:", data4.shape)

print("---- createFromArray usings pure numpy arrays")
coef4 = basis.createFromArray(mass3, data4, time=3.3)

coefs.add(coef4)

# Check the coefficient structure for the 4 added times
#
print("Times:", coefs.Times())

exit(0)
