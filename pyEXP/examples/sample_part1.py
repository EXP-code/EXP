import os
import yaml
import pyEXP

# I have the Better run here; obviously another use will want to
# change this a directory containing some snapshots of their own
#
os.chdir('/home/weinberg/Nbody/Better')

# Get the basis config
#
yaml_config = ""
with open('basis.yaml') as f:
    config = yaml.load(f, Loader=yaml.FullLoader)
    yaml_config = yaml.dump(config)

# Alternatively, you could construct this on the fly, e.g.
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
print('-'*60)
print('Read from file')
print('-'*60)
print(yaml_config)
print('-'*60)
print('Constructed from string')
print('-'*60)
print(bconfig)
print('-'*60)

# Construct the basis instance
#
basis   = pyEXP.basis.Basis.factory(yaml_config)

# Construct batches of files the particle reader.  One could use the
# parseStringList to create batches from a vector/list of files.  NB:
# a std::vector in C++ becomes a Python.list and vice versa
#
batches = pyEXP.read.ParticleReader.parseFileList('file.list', '')

# This will contain the coefficient container, need to start will a
# null instance to trigger construction
#
coefs   = None

for group in batches:

    print("file group is", group)

    # Make the reader for the desired type.  One could probably try to
    # do this by inspection but that's another project.
    #
    reader = pyEXP.read.ParticleReader.createReader('PSPout', group, 0, False);

    # Print the type list
    #
    print('The component names are:', reader.GetTypes())

    compname = 'dark halo'
    reader.SelectType(compname)
    print('Selected', compname)

    print('Call createFromReader at Time', reader.CurrentTime(), 'for', reader.CurrentNumber(), 'particles')

    coef = basis.createFromReader(reader)
    print("Created coef")

    # We need this stupid idiom here because None is not mapping to a
    # null pointer.  There is probably a way to do this.  Suggestions
    # anyone?
    #                          This is optional---+
    #                                             |
    if coefs is None:           #                 v
        coefs = pyEXP.coefs.Coefs.makecoefs(coef, compname)
    else:
        coefs.add(coef)

    print('Added coef')
    print('-'*60)

print('\nCompleted the file group list\n')

print('The coefficient time list is', coefs.Times())

# Now try some slices for rendering
#

times = coefs.Times()[0:3]
pmin  = [-1.0, -1.0, 0.0]
pmax  = [ 1.0,  1.0, 0.0]
grid  = [  40,   40,   0]

fields = pyEXP.field.FieldGenerator(times, pmin, pmax, grid)

surfaces = fields.slices(basis, coefs)

print("We now have the following [time field] pairs")
for v in surfaces:
    print('-'*40)
    for u in surfaces[v]:
        print("{:8.4f}  {}".format(v, u))

print("\nHere is the first one:")
for v in surfaces:
    for u in surfaces[v]:
        print('-'*40)
        print('----', u)
        print('-'*40)
        print(surfaces[v][u])
    break

# These could be make into images and so forth

# Okay, now try expMSSA
#
# Make some parameter flags as YAML.  The defaults should work fine
# for most people.  'chatty' turns on std::out diagnostics and
# 'output' is the prefix for written files.
#
flags ="""
---
# chatty: true
output: mytest
...
"""

config = {"dark halo": (coefs, [[1, 0, 0], [1, 0, 1], [1, 0, 2]], [])}
#          ^            ^      ^                                  ^
# tag -----+            |      |                                  |
# coefficient set-------+      |                                  |
# list of indices for MSSA-----+                                  |
# list of backround indices for reconstruction--------------------+
#
ssa = pyEXP.mssa.expMSSA(config, 5, 3, flags)

ev = ssa.eigenvalues()
print("The eigenvalues are:\n", ev)

pc = ssa.getPC();
print("The PC vectors are:\n", pc)

# Okay, now try a reconstruction
#
ssa.reconstruct()

newdata = ssa.getReconstructed(False)
print(type(newdata))
print(newdata)

# Try the kmeans analysis (not sure this is working yet)
#
ssa.kmeans()

# Test the PNG output
#
ssa.wcorrPNG()

# Make a subkey sequence
#
keylst = coefs.makeKeys([1])
print("Keys=", keylst)

# Try saving coefficients to an HDF5 file
#
coefs.WriteH5Coefs('test_better')

# Now try reading it in
#
coefs2 = pyEXP.coefs.Coefs.factory('test_better.h5')
print("Type is", coefs2.getGeometry())

# Now compare with the original
#
coefs2.CompareStanzas(coefs)

