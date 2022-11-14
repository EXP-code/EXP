#!/usr/bin/env python3

"""
Plot the orthgonality matrices using Sturm-Liouville
"""

import os, sys, getopt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pyEXP

def help(phrase: str) -> None:
   """Print some usage info"""
   print(phrase)

def main(prog, argv):
    """
    Plot the orthgonality matrices using Sturm-Liouville
    """

    model = 'SLgridSph.model'
    dir   = ''
    rmin  = 1.0e-4
    rmax  = 2.0
    Lmax  = 6
    Nmax  = 18
    knot  = 200

    phrase = prog + ': [-h] -f|--model=modelfile -r|--rmin arg -R|--rmax arg [-l|--lmax arg] [-n/--nmax arg] [-k|--knots arg]';

    try:
        opts, args = getopt.getopt(argv,"hf:r:R:L:N:",["model=","rmin=","rmax=","lmax=","nmax="])
    except getopt.GetoptError:
        help(phrase)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            help(phrase)
            sys.exit()
        elif opt in ("-d", "--dir"):
            dir = arg
        elif opt in ("-f", "--model"):
            model = arg
        elif opt in ("-r", "--rmin"):
            rmin = float(arg)
        elif opt in ("-R", "--rmax"):
            rmin = float(arg)
        elif opt in ("-l", "--lmax"):
            Lmax = int(arg)
        elif opt in ("-n", "--nmax"):
            Nmax = int(arg)
        elif opt in ("-k", "--knot"):
            knot = int(arg)


    if len(dir): os.chdir(dir)

    # 200 knots should be enough for most cases.  Try a smaller value
    # (e.g. 20) to observe the effect on the orthogonality ...

    # Construct the basis config for this model
    #
    bconfig = """
---
id: sphereSL
parameters :
  numr:  2000
  rmin:  {}
  rmax:  {}
  Lmax:  {}
  nmax:  {}
  scale: 0.0667
  modelname: {}
...
""".format(rmin, rmax, Lmax, Nmax, model)

    # Construct the basis instance
    #
    basis = pyEXP.basis.Basis.factory(bconfig)

    # Now compute the orthogonality matrices
    #
    ret   = basis.orthoCheck(knot)

    # Plot the matrices as images with a greyscale color map
    #
    fig   = plt.figure()
    ncol  = 4                   # Rows with 4 columns
    nrow  = int(Lmax/ncol)

    if ncol*nrow < Lmax: nrow += 1
    ax = fig.subplots(nrow, ncol).flatten()

    l = 0                       # Harmonic index counter

    for i in range(0, nrow):
        for j in range(0, ncol):
            if l<=Lmax:
                ax[i*ncol+j].imshow(ret[l], interpolation='nearest', cmap=cm.Greys_r)
                ax[i*ncol+j].set_aspect('equal')
                ax[i*ncol+j].set_title('l={}'.format(l))
                l += 1
            else:
                # Remove unused frames
                fig.delaxes(ax[i*ncol+j])

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
   main(sys.argv[0], sys.argv[1:])
