import os
import yaml
import pyEXP
import numpy as np
import matplotlib.pyplot as plt

# In this test, I assume that sample.py has already been run to
# generate a coefficient set.  This script points at that directory
# and does some additional analysis and plotting
#
os.chdir('/home/weinberg/Nbody/Better')

# Get the basis config
#
yaml_config = ""
with open('basis.yaml') as f:
    config = yaml.load(f, Loader=yaml.FullLoader)
    yaml_config = yaml.dump(config)

# Construct the basis instance
#
basis = pyEXP.basis.Basis.factory(yaml_config)

# Reread the coefs from the file
#
coefs = pyEXP.coefs.Coefs.factory('test_better.h5')

print("Got coefs for name=", coefs.getName())

# Now try some slices for rendering
#
times = coefs.Times()
pmin  = [-1.0, -1.0, 0.0]
pmax  = [ 1.0,  1.0, 0.0]
grid  = [  40,   40,   0]

fields = pyEXP.field.FieldGenerator(times, pmin, pmax, grid)

surfaces = fields.slices(basis, coefs)

print("We now have the following [time field] pairs")
final = 0.0
for v in surfaces:
    print('-'*40)
    for u in surfaces[v]:
        print("{:8.4f}  {}".format(v, u))
        final = v

# Print the potential image at the final time (I think there is a
# fencepost issue in this grid, no matter).
x = np.linspace(pmin[0], pmax[0], grid[0])
y = np.linspace(pmin[1], pmax[1], grid[1])
xv, yv = np.meshgrid(x, y)

cont1 = plt.contour(xv, yv, surfaces[final]['p'].transpose(), colors='k')
plt.clabel(cont1, fontsize=9, inline=True)
cont2 = plt.contourf(xv, yv, surfaces[final]['p'].transpose())
plt.colorbar(cont2)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Potential at T={}'.format(final))
plt.show()

# Okay, now try expMSSA
#

# Make a subkey sequence
#
keylst = coefs.makeKeys([1])
print("Keys=", keylst)

config = {"dark halo": (coefs, keylst, [])}

window = int(len(coefs.Times())/2)
npc = 10

print("Window={} PC number={}".format(window, npc))

ssa = pyEXP.mssa.expMSSA(config, window, npc)

ev = ssa.eigenvalues()

plt.plot(ev, 'o-')
plt.xlabel("Index")
plt.ylabel("EV")
plt.title("Eigenvalues by index")
plt.show()

times = coefs.Times();
pc    = ssa.getPC();

rows, cols = pc.shape

for i in range(cols):
    plt.plot(times[0:rows], pc[:,i], '-', label="{:d}".format(i))

plt.xlabel('Time')
plt.ylabel('PC')
plt.legend()
plt.title("Principal components (left-singular vectors)")
plt.show()

# Okay, now try a reconstruction
#
ssa.reconstruct()

newdata = ssa.getReconstructed()
print('newdata is a', type(newdata))

# Try the kmeans analysis (not sure this is working correctly yet,
# although it used to work and nothing has changed)
#
ssa.kmeans()

# Test the PNG output
#
ssa.wcorrPNG()

