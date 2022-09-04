import os
import yaml
import time
import pyEXP
import numpy as np
import matplotlib.pyplot as plt

# I'm reading a previously computed native cylindrical coefficient set,
# saving it as HDF5, and executing some MSSA computations
#
os.chdir('/home/weinberg/Nbody/Better')

# Get the basis config
#
yaml_config = ""
with open('disk_basis.yaml') as f:
    config = yaml.load(f, Loader=yaml.FullLoader)
    yaml_config = yaml.dump(config)

# Construct the basis instance
#
basis = pyEXP.basis.Basis.factory(yaml_config)

# Make the coefficients by the factory method
#
coefs = pyEXP.coefs.Coefs.factory('outcoef.star.run0')
print('The coefficient time list is', coefs.Times())

# Try saving coefficients to an HDF5 file
#
# coefs.WriteH5Coefs('test_disk')
# print('Saved coefficients as "test_disk.h5"')

# Now try some slices for rendering
#
times = coefs.Times()[-2:]
pmin  = [-0.05, -0.05, 0.0]
pmax  = [ 0.05,  0.05, 0.0]
grid  = [    20,   20,   0]

print('Creating surfaces with times:', times)

fields = pyEXP.field.FieldGenerator(times, pmin, pmax, grid)

print('Created fields instance')

surfaces = fields.slices(basis, coefs)

print('Created surfaces')

print("We now have the following [time field] pairs")
final = 0.0
for v in surfaces:
    print('-'*40)
    for u in surfaces[v]:
        print("{:8.4f}  {}".format(v, u))
        final = v

# Print the potential image at the final time
# 
nx = surfaces[final]['p'].shape[0]
ny = surfaces[final]['p'].shape[1]

x = np.linspace(pmin[0], pmax[0], nx)
y = np.linspace(pmin[1], pmax[1], ny)
xv, yv = np.meshgrid(x, y)

cont1 = plt.contour(xv, yv, surfaces[final]['p'], colors='k')
plt.clabel(cont1, fontsize=9, inline=True)
cont2 = plt.contourf(xv, yv, surfaces[final]['p'])
plt.colorbar(cont2)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Potential at T={}'.format(final))
plt.show()

cont1 = plt.contour(xv, yv, surfaces[final]['fr'], colors='k')
plt.clabel(cont1, fontsize=9, inline=True)
cont2 = plt.contourf(xv, yv, surfaces[final]['fr'])
plt.colorbar(cont2)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Radial force at T={}'.format(final))
plt.show()

cont1 = plt.contour(xv, yv, surfaces[final]['ft'], colors='k')
plt.clabel(cont1, fontsize=9, inline=True)
cont2 = plt.contourf(xv, yv, surfaces[final]['ft'])
plt.colorbar(cont2)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Vertical force at T={}'.format(final))
plt.show()

# Okay, now try expMSSA
#
# Make some parameter flags as YAML.  The defaults should work fine
# for most people.  'chatty' turns on std::out diagnostics and
# 'output' is the prefix for written files.
#
flags ="""
---
# chatty: on
output: mytest2
groups: [[0, 1], [2, 3]]
...
"""

# Make a subkey sequence
#
keylst = coefs.makeKeys([2])
print("All m=2 keys=", keylst)

# Make some custom [m, n] pairs
keylst = [[2, 0], [2, 1], [2, 2], [2, 3], [2, 4], [2, 5], [2, 6], [2, 7]]
print("Custom keys=", keylst)

config = {"star disk": (coefs, keylst, [])}

startTime = time.time()
ssa = pyEXP.mssa.expMSSA(config, int(len(coefs.Times())/2), 20, flags)

ev = ssa.eigenvalues()
print('Computed eigenvalues in {:6.2f} seconds'.format(time.time() - startTime))

plt.plot(ev, 'o-')
plt.xlabel("Index")
plt.ylabel("EV")
plt.title("Eigenvalues by index")
plt.show()

times = coefs.Times()
pc = ssa.getPC()

rows, cols = pc.shape

for i in range(min(cols,4)):
    plt.plot(times[0:rows], pc[:,i], '-', label="{:d}".format(i))

plt.xlabel('Time')
plt.ylabel('PC')
plt.legend()
plt.title("Principal components (left-singular vectors)")
plt.show()

# Okay, now try a reconstruction
#
startTime = time.time()
print('Calling reconstruction')
ssa.reconstruct([0, 1])
print('Reconstruction took {:6.2f} seconds'.format(time.time() - startTime))

# newdata = ssa.getReconstructed(False) # <---replace with reconstructed
newdata = ssa.getReconstructed(True) # <---zero all but reconstructed
print('newdata is a', type(newdata))

surfaces = fields.slices(basis, newdata['star disk'])

print('Created surfaces from reconstruction')

print("We now have the following [time field] pairs")
final = 0.0
for v in surfaces:
    print('-'*40)
    for u in surfaces[v]:
        print("{:8.4f}  {}".format(v, u))
        final = v

# Print the potential image at the final time (I think there is a
# fencepost issue in this grid, no matter).
nx = surfaces[final]['p'].shape[0]
ny = surfaces[final]['p'].shape[1]

x = np.linspace(pmin[0], pmax[0], nx)
y = np.linspace(pmin[1], pmax[1], ny)
xv, yv = np.meshgrid(x, y)

cont1 = plt.contour(xv, yv, surfaces[final]['p'], colors='k')
plt.clabel(cont1, fontsize=9, inline=True)
cont2 = plt.contourf(xv, yv, surfaces[final]['p'])
plt.colorbar(cont2)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Potential at T={}'.format(final))
plt.show()

cont1 = plt.contour(xv, yv, surfaces[final]['fr'], colors='k')
plt.clabel(cont1, fontsize=9, inline=True)
cont2 = plt.contourf(xv, yv, surfaces[final]['fr'])
plt.colorbar(cont2)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Radial force at T={}'.format(final))
plt.show()

cont1 = plt.contour(xv, yv, surfaces[final]['ft'], colors='k')
plt.clabel(cont1, fontsize=9, inline=True)
cont2 = plt.contourf(xv, yv, surfaces[final]['ft'])
plt.colorbar(cont2)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Vertical force at T={}'.format(final))
plt.show()


# Try the kmeans analysis (not sure this is working correctly yet,
# although it used to work and nothing has changed)
#
print('Calling k-means')
ssa.kmeans()

# Test the PNG output
#
ssa.wcorrPNG()

