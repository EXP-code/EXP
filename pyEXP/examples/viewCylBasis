#!/usr/bin/env python3

"""
Plot the cylindrical basis functions
"""

import os, sys, getopt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import gridspec
import pyEXP

def help(phrase: str) -> None:
   """Print some usage info"""
   print(phrase)

def main(prog, argv) -> int:
   """
   Plot the cylindrical basis functions
   """

   cfile = ''
   dir   = ''
   Rmin  = 0.0
   Rmax  = 0.1
   Zmax  = 0.03
   Rnum  = 100
   Znum  = 40

   phrase = prog + ': [-h] -c|--cache=file [-d|--dir cache_directory] [-R|--rmax val] [-Z|--zmax val]';

   try:
      opts, args = getopt.getopt(argv,"hc:d:R:Z:",["cache=","dir=","rmax=","zmax="])
   except getopt.GetoptError:
      help(phrase)
      sys.exit(2)

   for opt, arg in opts:
      if opt == '-h':
         help(phrase)
         sys.exit()
      elif opt in ("-d", "--dir"):
         dir = arg
      elif opt in ("-c", "--cache"):
         cfile = arg
      elif opt in ("-R", "--rmax"):
         Rmax = float(arg)
      elif opt in ("-Z", "--zmax"):
         Zmax = float(arg)

   # Check for directory
   if len(dir):
      os.chdir(dir)

   # Get the cache data
   #
   params = pyEXP.basis.Cylindrical.cacheInfo(cfile)

   Mmax = int(params['mmax'])

   bconfig = """
---
id: cylinder
parameters:
  acyl: {}
  hcyl: {}
  mmax: {}
  nmax: {}
  ncylorder: {}
  ncylnx: {}
  ncylny: {}
  eof_file: {}
...
""".format(params['ascl'], params['hscl'], params['mmax'], params['nmax'],
           params['norder'], params['numx'], params['numy'], cfile)
   
   # Construct the basis instance
   #
   basis = pyEXP.basis.Basis.factory(bconfig)

   Norder = int(params['norder'])

   # Plot the matrices as images with a greyscale color map
   #
   ncol  = 4                   # Rows with 4 columns
   nrow  = int(Norder/ncol)
    
   if ncol*nrow <= Norder: nrow += 1

   R = np.linspace(0.0, Rmax, Rnum)
   Z = np.linspace(-Zmax, Zmax, Znum)

   xv, yv = np.meshgrid(R, Z)

   grid = basis.getBasis(Rmin, Rmax, Rnum, -Zmax, Zmax, Znum)

   for m in range(Mmax+1):
      n = 0
      fig = plt.figure(figsize=(15, 15)) 
      gs = gridspec.GridSpec(nrow, ncol, width_ratios=[1, 1, 1, 1],
                             wspace=0.01, hspace=0.02) 

      for row in range(nrow):
         for col in range(ncol):
            if n<Norder:
               ax = plt.subplot(gs[n])
               plt.axis('off')
               cx = plt.contourf(xv, yv, grid[m][n].transpose())
               n += 1

      plt.annotate("m={}".format(m), (0.5, 0.95), xycoords='figure fraction')
      plt.show()

if __name__ == "__main__":
   main(sys.argv[0], sys.argv[1:])
