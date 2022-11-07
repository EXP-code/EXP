import os, sys
from mpi4py import MPI

import pyEXP
import numpy as np
import pickle

from os.path import exists

#
# This script makes an image histogram in parallel using MPI.  This
# can easily be adapted for whatever snapshots you have and could be
# run on a cluster
#

# Run command is: "mpirun -np N python3 make_histo_MPI.py"
# where N is the number of processes.  For Slurm allocations, you can
# leave off "-np N" as usual.
#

# Make the file list for the snapshot sequence
#

# For making a unique fixed-point time for a dictionary key
#
def fixTime(time):
    fixedD = 100000.0;
    return int(fixedD*time+0.5)/fixedD

if __name__ == "__main__":

    # Parameters
    #
    beg_seq = 0
    end_seq = 633

    # Get basic information about the MPI communicator
    #
    world_comm = MPI.COMM_WORLD
    world_size = world_comm.Get_size()
    my_rank    = world_comm.Get_rank()

    # NB: we are not using standard MPI commands on the Python side,
    # but invoking mpi4py initializes MPI so it can be used correctly
    # by the C++ routines

    # Just for info and to convince yourself check that MPI is working
    #
    print("World size is {} and my rank is {}".format(world_size, my_rank))


    # Make a phase space file list
    #
    file_list = []
    for i in range(beg_seq, end_seq):
        file_list.append('SPL.run2Fd_3.{:05d}'.format(i))

    # Construct batches of files the particle reader.  One could use the
    # parseStringList to create batches from a vector/list of files.  NB:
    # a std::vector in C++ becomes a Python.list and vice versa
    #
    batches = pyEXP.read.ParticleReader.parseStringList(file_list, '')
    xy = {}
    xz = {}
    yz = {}

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
        
        tim = fixTime(reader.CurrentTime())
        gd[tim] = fg.histo(reader)

    if my_rank==0:
        keys = list(gd.keys())
        print('First time={}  last time={}'.format(keys[0], keys[-1]))
        file = open('imagePickle', 'wb')
        db = {'image': gd, 'lower': lower, 'upper': upper, 'ngrid': ngrid}
        pickle.dump(db, file)
        file.close()
