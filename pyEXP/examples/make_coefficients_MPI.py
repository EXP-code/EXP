import os
import time                     # Used for timing the coefficient construction
import pyEXP
from mpi4py import MPI

from os.path import exists

#
# This script makes an HDF5 coefficient set from some Gadget files in
# parallel using MPI.  This can easily be adapted for whatever
# snapshots you have and could be run on a cluster
#

# There is a companion Jupyter Notebook called 'test_read_gadget.ipynb'
# that reads in the coefficients made by this script to check that
# they seem sane.

# Usage:
#
# Run command is: "mpirun -np N python3 make_coefficients_MPI.py"
# where N is the number of processes.  For Slurm allocations, you can
# leave off "-np N" as usual.
#

if __name__ == "__main__":

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

    # Now switch the working directory where my Gadget simulation
    # lives.  Change this to your working directory.
    #
    os.chdir('/media/weinberg/Simulation data/Nbody/Sphere/RunG')

    # Make the spherical basis config.
    #
    halo_config = """
id          : sphereSL
parameters  :
  numr      : 1000
  rmin      : 0.000011
  rmax      : 1.99
  Lmax      : 6
  nmax      : 18
  rs        : 0.05
  modelname : mw_halo.model
"""

    # Construct the basis instance
    #
    halo_basis = pyEXP.basis.Basis.factory(halo_config)

    # Make the file list.  Here, I'm making the first 300 snaps.
    #
    file_list = []
    # for i in range(0, 1593): file_list.append('snapshot_{:04d}.hdf5'.format(i))
    for i in range(0, 300): file_list.append('snapshot_{:04d}.hdf5'.format(i))

    # Construct batches of files the particle reader.  One could use the
    # parseStringList to create batches from a vector/list of files.  NB:
    # a std::vector in C++ becomes a Python.list and vice versa
    #
    batches = pyEXP.read.ParticleReader.parseStringList(file_list, '')

    # This will contain the coefficient container, need to start will a
    # null instance to trigger construction
    #
    coefs = None

    for group in batches:
        okay = True
        for f in group:
            if not exists(f): okay = False
        
        # Skip a file that does not exist in the sequence without
        # going belly up.  I notice that Gadget-2 seems to not be
        # exactly sequential by 1 through a restart?
        #
        if not okay: continue

        if my_rank==0: print("file group is", group)

        # Make the reader for the desired type.  One could probably try to
        # do this by inspection but that's another project.
        #
        reader = pyEXP.read.ParticleReader.createReader('GadgetHDF5', group, 0, False);

        # Print the type list
        #
        if my_rank==0: print('The component names are:', reader.GetTypes())

        compname = 'Halo'
        reader.SelectType(compname)
        if my_rank==0: print('Selected', compname)

        startTime = time.time()
        coef = halo_basis.createCoefficients(reader)
        if my_rank==0:
            print('Created createCoefficients at Time', reader.CurrentTime(), 'for', reader.CurrentNumber(), 'particles in', time.time() - startTime, 'seconds')

        # We need this stupid idiom here because None is not mapping to a
        # null pointer.  There is probably a way to do this.  Suggestions
        # anyone?
        #                          This is optional---+
        #                                             |
        if coefs is None:           #                 v
            coefs = pyEXP.coefs.Coefs.makecoefs(coef, compname)
        else:
            coefs.add(coef)

        if my_rank==0:
            print('Added coef to container')
            print('-'*60)

    if my_rank==0:
        print('\nCompleted the file group list\n')
        print('The coefficient time list is', coefs.Times())
        print('Save the coefficients to a HDF5 file')

        # You can call the file something convenient.  The suffix 'h5'
        # will be appended. You only want the root process to write
        # the file.
        #
        coefs.WriteH5Coefs('RunG_halo') 
        
