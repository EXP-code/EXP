import os
import yaml
import pyEXP

os.chdir('/home/weinberg/Nbody/Better')

# Get the basis config
yaml_config = ""
with open('basis.yaml') as f:
    config = yaml.load(f, Loader=yaml.FullLoader)
    yaml_config = yaml.dump(config)

# Alternatively, you could construct this on the fly, e.g.
bconfig = """
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
basis   = pyEXP.basis.Basis.factory(yaml_config)

# Construct the particle reader
batches = pyEXP.read.ParticleReader.parseFileList('file.list', '')

# This will contain the coefficient container
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

    print('Call createCoefficients at Time', reader.CurrentTime(), 'for', reader.CurrentNumber(), 'particles')

    coef = basis.createCoefficients(reader)
    print("Created coef")

    # We need this stupid idiom here because None is not mapping to a
    # null pointer.  There is probably a way to fix this.
    if coefs is None:
        coefs = pyEXP.coefs.Coefs.makecoefs(coef)
    else:
        coefs.add(coef)

    print('Added coef')
    print('-'*60)

print('\nCompleted the file group list\n')

print('The coefficient time list is', coefs.Times())

# Now try some slices
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

# Okay, now try expMSSA

# Make some parameter flags as YAML
flags ="""
---
chatty: true
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
print("Eigenvalues:", ev)

pc = ssa.getPC();
print("PC:", pc)

# Okay, now try a reconstruction
ssa.reconstruct([0, 1])

newdata = ssa.getReconstructed(False)
print(type(newdata))
print(newdata)

# Try the kmeans
ssa.kmeans()

# Test the PNG output
ssa.wcorrPNG()

# Make a subkey sequence
keylst = coefs.makeKeys([1])
print("Keys=", keylst)

# Try saving coefficients to an HDF5 file
coefs.WriteH5Coefs('testC')

# Now try reading it in
coefs2 = pyEXP.coefs.Coefs.factory('testC.h5')
print("Type is", coefs2.getCoefType())

# Now compare with the original
coefs2.CompareStanzas(coefs)

