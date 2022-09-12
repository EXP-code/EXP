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

    # Parameters
    #
    h5file  = 'RunG_halo_test'
    ctrfile = 'new.centers'
    beg_seq = 0
    end_seq = 400
    nskip   = 20

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

    # Make the file list for the snapshot sequence assumed have the
    # format: snapshot_XXXX.hdf5.  Change as necessary.
    #
    file_list = []
    for i in range(beg_seq, end_seq):
        file_list.append('snapshot_{:04d}.hdf5'.format(i))

    # Construct batches of files the particle reader.  One could use the
    # parseStringList to create batches from a vector/list of files.  NB:
    # a std::vector in C++ becomes a Python.list and vice versa
    #
    batches = pyEXP.read.ParticleReader.parseStringList(file_list, '')

    # This will contain the coefficient container, need to start will a
    # null instance to trigger construction
    #
    coefs = None

    centers = []

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
        if my_rank==0: print('We selected:', compname)

        # This computes an expansion center from a mean density
        # weighted position.  You could compute and cache the center
        # array . . . or supply it in a different way
        #
        startTime = time.time()
        center = pyEXP.util.getDensityCenter(reader, nskip)
        #                                            ^
        #                                            |
        # Choose every nskip particle for sample ----+
        # This is 10^6 samples for these snaps and nskip=20
        #
        if my_rank==0:
            print('Created center in {:4.2f} seconds'.
                  format(time.time() - startTime))
            print('Center is:', center)
            centers.append(center)

        # Now compute the coefficients using this center
        #
        startTime = time.time()
        coef = halo_basis.createCoefficients(reader, center)
        if my_rank==0:
            print('Created createCoefficients at Time {:5.3f} for {} particles '
                  'in {:4.2f} seconds'.
                  format(reader.CurrentTime(), reader.CurrentNumber(),
                         time.time() - startTime))

        # We need this stupid idiom here because None is not mapping to a
        # null pointer.  There is probably a way to do this.  Suggestions
        # anyone?
        #                          This is optional---+
        #                                             |
        if coefs is None:           #                 v
            coefs = pyEXP.coefs.Coefs.makecoefs(coef, compname)
        coefs.add(coef)

        if my_rank==0:
            print('Added coef to container')
            print('-'*60)

    if my_rank==0:
        print('\nCompleted the file group list\n')
        print('The coefficient time list is', coefs.Times())

        # You can call the file something convenient.  The suffix 'h5'
        # will be appended. You only want the root process to write
        # the file.
        #
        if exists(h5file + '.h5'):
            coefs.ExtendH5Coefs(h5file) # Update an existing HDF5
            print('Saved the coefficients to an existing HDF5 file')
        else:
            coefs.WriteH5Coefs(h5file) # Create a new HDF5
            print('Saved the coefficients to a new HDF5 file')

        # Save the center positions
        #
        with open(ctrfile, 'a') as f:
            times = coefs.Times()
            for i in range(len(times)):
                line = '{:13.6e} {:13.6e} {:13.6e} {:13.6e}\n'.format(times[i], centers[i][0], centers[i][1], centers[i][2])
                f.write(line)
                
