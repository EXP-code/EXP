#!/usr/bin/env python
# coding: utf-8

import pyEXP
import numpy as np
from mpi4py import MPI

def pprint(str="", end="\n", comm=MPI.COMM_WORLD):
    """Print for MPI parallel programs: Only rank 0 prints *str*."""
    if comm.rank == 0:
        print(str+end, end=' ')
        
# Get the basis config
#
disk_config = """
id           : cylinder
parameters   :
  acyl       : 1.0
  hcyl       : 0.1
  lmax       : 48
  mmax       : 10
  nmax       : 48
  ncylorder  : 32
  ncylodd    : 6
  ncylnx     : 256
  ncylny     : 128
  ncylr      : 2000
  rnum       : 200
  pnum       : 1
  tnum       : 80
  rcylmin    : 0.001
  rcylmax    : 20
  ashift     : 0
  logr       : true
  density    : true
  expcond    : true
  deproject  : true
  eof_file   : .eof.cache.file_new
  ignore     : true
  vflag      : 16
"""

# Initialize MPI
#
world_comm = MPI.COMM_WORLD
world_size = world_comm.Get_size()
my_rank    = world_comm.Get_rank()

# Just for info and to convince yourself check that MPI is working on
# your system
#
pprint("============================================================================")
pprint("Compute a Cylindrical basis on multiple nodes or on a laptop with MPI")
pprint("============================================================================")

world_comm.Barrier()
print("World size is {} and my rank is {}".format(world_size, my_rank))

world_comm.Barrier()
pprint("============================================================================")

# Begin calculation and start stopwatch
#
t_start = MPI.Wtime()

# Construct the basis instance
#
disk_basis = pyEXP.basis.Basis.factory(disk_config)

world_comm.Barrier()

pprint("============================================================================")
pprint("Calculation finished")

# Stop stopwatch
#
world_comm.Barrier()
t_diff = MPI.Wtime() - t_start

pprint("Computed basis in {}".format(t_diff))
pprint("============================================================================")
