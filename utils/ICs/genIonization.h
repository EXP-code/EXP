R"(#!/usr/bin/python3

import os
import re
import sys
import operator

# Try to prevent ChiantiPy from grabbing the mouse
#
import matplotlib as mpl
mpl.use('Agg')

from numpy import *
import ChiantiPy.core as ch
from optparse import OptionParser

def accumulate(iterable, func=operator.add):
    'Return running totals'
    # accumulate([1,2,3,4,5]) --> 1 3 6 10 15
    # accumulate([1,2,3,4,5], operator.mul) --> 1 2 6 24 120
    it = iter(iterable)
    total = next(it)
    yield total
    for element in it:
        total = func(total, element)
        yield total

def getIoneq(ofile, T=1.0e+04, numBeg=1, numEnd=2):
        """ Use ChiantiPy to get the ionization equilibrium"""
        file = open(ofile, 'w')
        file.write('{0:4d} {1:4d} {2:15.6e}\n'.format(numBeg, numEnd, T))
        for n in range(numBeg, numEnd+1):
                atom = ch.ioneq(n)
                atom.calculate([T, T])
                f = atom.Ioneq
                z = list(accumulate(f))
                z /= z[-1]      # Normalization, should not be needed
                file.write('{0:4d} {1:4d}\n'.format(n, len(z)))
                for v in f:
                    file.write('{0:15.6e}'.format(v[0]))
                file.write('\n')
                for v in z:
                    file.write('{0:15.6e}'.format(v[0]))
                file.write('\n')

def main():
        """
        Get the ionization equilibrium for some number of elements
        Print to a file
        """
        usage  = "usage: %prog [options] file"
        parser = OptionParser(usage)

        parser.add_option("-1", "--nBeg", default=1,
                      action="store", type="int", dest="nbeg",
                      help="beginning atomic number")
        parser.add_option("-2", "--nEnd", default=2,
                      action="store", type="int", dest="nend",
                      help="ending atomic number")
        parser.add_option("-T", "--temp", default=10000.0,
                      action="store", type="float", dest="temp",
                      help="temperature")
        parser.add_option("-o", "--output", default='ioneq.dat',
                      action="store", type="string", dest="ofile",
                      help="output file")

        (opt, args) = parser.parse_args()

        # Call the main routine
        getIoneq(opt.ofile, T=opt.temp, numBeg=opt.nbeg, numEnd=opt.nend)

if __name__ == "__main__":
        main()
)"
