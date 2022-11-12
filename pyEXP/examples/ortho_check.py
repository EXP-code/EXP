# Using Sturm-Liouville

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pyEXP

# I'm using the model from the EXP example Better run here.  For
# others, you will want to change this a directory containing a
# spherical model for the basis or put the file in the working
# directory and omit the following line
#
os.chdir('/home/weinberg/Projects/EXP/examples/Better')

# Some parameters
#
rmin = 0.0001                   # Minimum radius
rmax = 1.99                     # Maximum radius
Lmax = 6                        # Maximum harmonic order
Nmax = 24                       # Maximum radial order
knot = 200                      # Number of quadrature knots

# 200 knots should be enough for most cases.  Try a smaller value (e.g. 20) to
# observe the effect on the orthogonality ...

# Construct the basis config for this model
#
bconfig = """
---
id: sphereSL
parameters :
  numr:  2000
  rmin:  {}
  rmax:  {}
  Lmax:  {}
  nmax:  {}
  scale: 0.0667
  modelname: SLGridSph.model
...
""".format(rmin, rmax, Lmax, Nmax)

# Construct the basis instance
#
basis = pyEXP.basis.Basis.factory(bconfig)

# Now compute the orthogonality matrices
#
ret   = basis.orthoCheck(knot)

# Plot the matrices as images with a greyscale color map
#
fig   = plt.figure()
ncol  = 4                       # Rows with 4 columns
nrow  = int(Lmax/ncol)

if ncol*nrow < Lmax: nrow += 1
ax = fig.subplots(nrow, ncol).flatten()

l = 0                           # Harmonic index counter

for i in range(0, nrow):
    for j in range(0, ncol):
        if l<=Lmax:
            ax[i*ncol+j].imshow(ret[l], interpolation='nearest', cmap=cm.Greys_r)
            ax[i*ncol+j].set_aspect('equal')
            ax[i*ncol+j].set_title('l={}'.format(l))
            l += 1
        else:
            # Remove unused frames
            fig.delaxes(ax[i*ncol+j])

plt.tight_layout()
plt.show()
