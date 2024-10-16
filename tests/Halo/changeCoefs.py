import os
import pyEXP
import numpy as np

# Read the coefs from the file created by the expNbodyTest
#
coefs = pyEXP.coefs.Coefs.factory('outcoef.halo.run0')

print("Got coefs for name=", coefs.getName())


times = coefs.Times()
data  = coefs.getAllCoefs()

print("The data sizes and coefficient orders are:")
print("M orders, N orders, Times:", data.shape)


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

# Let's zero all of the odd order data
#
for k in range(data.shape[0]):
    l, m = pyEXP.basis.SphericalSL.invI(k)
    if l%2!=0 or m%2!=0: data[k,:,:] *= 0.0

# Reset the coefficient data
#
for i in range(data.shape[2]):
    coefs.setMatrix(times[i], data[:,:,i])

# Check that odd coefficients are really zero
data1 = coefs.getAllCoefs()
minZero =  1.0e30
maxZero = -1.0e30
for k in range(data1.shape[0]):
    l, m = pyEXP.basis.SphericalSL.invI(k)
    if l%2!=0 or m%2!=0:
        for v in data1[k,:,:].flatten():
           minZero = min([minZero, abs(v)])
           maxZero = max([maxZero, abs(v)])

print('Zero test: min, max: {:13.7e}, {:13.7e}'.format(minZero, maxZero))

if minZero < -1.0e-18 or maxZero > 1.0e-18:
    exit(1)
else:
    exit(0)
