import os
import yaml
import pyEXP

# I'm using the EXP example Better run here; obviously another use
# will want to change this a directory containing a spherical model
# for the basis or put the file in the working directory and omit the
# following line
#
os.chdir('/home/weinberg/Projects/EXP/examples/Better')

# Construct the basis config for this model
#
bconfig = """
---
id: sphereSL
parameters :
  numr:  2000
  rmin:  0.0001
  rmax:  1.95
  Lmax:  4
  nmax:  10
  scale: 0.0667
  modelname: SLGridSph.model
...
"""

# Construct the basis instance
#
basis   = pyEXP.basis.Basis.factory(bconfig)

# Now compute the orthogonality matrices
#
ret     = basis.orthoCheck()

# Plot the matrices as images
#
Lmax = 4
fig  = plt.figure()
ncol = 3
nrow = int(Lmax/ncol)

if ncol*nrow < Lmax: nrow += 1
ax  = fig.subplots(nrow, ncol)

l = 0
for i in range(1, nrow+1):
    for j in range(1, ncol+1):
        if l<=Lmax:
            ax[i, j].imshow(a, interpolation='nearest', cmap=cm.Greys_r)
            ax[i, j].set_aspect('equal')
            l += 1

plt.show()
