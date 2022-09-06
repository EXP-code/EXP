import os
import yaml
import time
import pyEXP
import numpy as np
import matplotlib.pyplot as plt

os.chdir('/home/weinberg/Nbody/Better')

# Make the halo basis config
halo_config="""
---
id: SphereSL
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

# Make the disk basis config
#
disk_config = """
---
id: Cylinder
parameters:
  acyl: 0.01
  hcyl: 0.001
  lmax: 32
  mmax: 6
  nmax: 32
  ncylorder: 8
  ncylnx: 256
  ncylny: 128
  rnum: 200
  pnum: 0
  tnum: 80
  ashift: 0.5
  vflag: 0
  logr: false
  density: false
  eof_file: .eof.cache.run0
...
"""

# Construct the basis instances
#
halo_basis = pyEXP.basis.Basis.factory(halo_config)
disk_basis = pyEXP.basis.Basis.factory(disk_config)

print(halo_basis)

# Get the two basis grids
#
lrmin = -3.0
lrmax = 0.5
rnum  = 200
halo_grid = halo_basis.getBasis(lrmin, lrmax, rnum)

Rmin = 0.0
Rmax = 0.1
Rnum = 100
Zmin = -0.03
Zmax =  0.03
Znum = 40

disk_grid = disk_basis.getBasis(Rmin, Rmax, Rnum, Zmin, Zmax, Znum)

# Plot some halo basis functions
#
r = np.linspace(lrmin, lrmax, rnum)
r = np.power(10.0, r)

for l in range(3):
    for n in range(5):
        plt.semilogx(r, halo_grid[l][n], '-', label="n={}".format(n))
    plt.xlabel('r')
    plt.ylabel('potential')
    plt.title('l={}'.format(l))
    plt.legend()
    plt.show()
    
R = np.linspace(Rmin, Rmax, Rnum)
Z = np.linspace(Zmin, Zmax, Znum)

xv, yv = np.meshgrid(R, Z)

for m in range(3):
    for n in range(5):
        # Tranpose for contourf
        cx = plt.contourf(xv, yv, disk_grid[m][n].transpose())
        plt.xlabel('R')
        plt.ylabel('Z')
        plt.title('m, n={}, {}'.format(m, n))
        plt.colorbar(cx)
        plt.show()
    
