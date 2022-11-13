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

Mmax = 6

bconfig = """
---
id: cylinder
parameters:
  acyl: 0.01
  hcyl: 0.001
  lmax: 32
  mmax: 6
  nmax: 32
  ncylorder: 8
  ncylnx: 256
  ncylny: 128
  eof_file: .eof.cache.run0
...
"""

# Construct the basis instance
#
basis = pyEXP.basis.Basis.factory(bconfig)

# Now compute the orthogonality matrices
#
ret   = basis.orthoCheck()

# Plot the matrices as images with a greyscale color map
#
fig   = plt.figure()
ncol  = 4                       # Rows with 4 columns
nrow  = int(Mmax/ncol)

if ncol*nrow < Mmax: nrow += 1
ax = fig.subplots(nrow, ncol).flatten()

M = 0                           # Harmonic index counter

for i in range(0, nrow):
    for j in range(0, ncol):
        if M<=Mmax:
            ax[i*ncol+j].imshow(ret[M], interpolation='nearest', cmap=cm.Greys_r)
            ax[i*ncol+j].set_aspect('equal')
            ax[i*ncol+j].set_title('M={}'.format(M))
            M += 1
        else:
            # Remove unused frames
            fig.delaxes(ax[i*ncol+j])

plt.tight_layout()
plt.show()
