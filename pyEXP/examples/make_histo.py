import os, sys
import pyEXP
import numpy as np
import pickle

from os.path import exists

# Make the file list for the snapshot sequence
#
beg_seq = 0
end_seq = 6
file_list = []
for i in range(beg_seq, end_seq):
    file_list.append('SPL.run2Fd_2.{:05d}'.format(i))
#                     ^
#                     |
#   Change this depending on your phase-space type

# Construct batches of files the particle reader.  One could use the
# parseStringList to create batches from a vector/list of files.  NB:
# a std::vector in C++ becomes a Python.list and vice versa
#
batches = pyEXP.read.ParticleReader.parseStringList(file_list, '')
xy = {}
xz = {}
yz = {}

# For a unique map with fixed-point time as a key
#
def getTime(time):
    fixedD = 100000.0;
    return int(fixedD*time+0.5)/fixedD


times = []
lower = [-0.03, -0.03, -0.03]
upper = [ 0.03,  0.03,  0.03]
ngrid = [  80,   80,   80]

fg = pyEXP.field.FieldGenerator(times, lower, upper, ngrid)
gd = {}

for group in batches:

    okay = True
    for f in group:
        if not exists(f):
            okay = False
        
    if not okay: continue

    # Make the reader for the desired type.  One could probably try to
    # do this by inspection but that's another project.
    #
    reader = pyEXP.read.ParticleReader.createReader('PSPspl', group, 0, False);

    compname = 'star'
    reader.SelectType(compname)
    
    tim = getTime(reader.CurrentTime())
    gd[tim] = fg.histo(reader)

for v in gd: print(v)

file = open('imagePickle', 'wb')
db = {'image': gd, 'lower': lower, 'upper': upper, 'ngrid': ngrid}
pickle.dump(db, file)
file.close()
