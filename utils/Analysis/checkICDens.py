#!/usr/bin/python

# For Python 3 compatibility
#
from __future__ import absolute_import, division, print_function, unicode_literals

import sys, os, argparse
import astroval as C

# Argument parsing
#
parser = argparse.ArgumentParser(description='Evaluate global statistics for body file created by makeIonIC')
parser.add_argument('-b', '--body',    default='gas.bod', help='Name of body file')
parser.add_argument('-m', '--msun',    default=1.0, type=float, help='Number of M_sun per unit mass')
parser.add_argument('-l', '--lpc',     default=1.0, type=float, help='Number of parsecs per unit length')

args = parser.parse_args()

file = open(args.body)
file.readline()
mval = 1.0e32
mtot = 0.0
minD = [mval, mval, mval, mval, mval, mval]
maxD = [0.0,  0.0,  0.0,  0.0,  0.0,  0.0]
avgP = [0.0, 0.0, 0.0]
avgV = [0.0, 0.0, 0.0]
for line in file:
    v = line.split()
    m = float(v[0])
    mtot += m
    for i in range(0,6):
        minD[i] = min(minD[i], float(v[i+1]))
        maxD[i] = max(maxD[i], float(v[i+1]))
    for i in range(0,3):
        avgP[i] += m*float(v[i+1])
        avgV[i] += m*float(v[i+4])
        
mfac = args.msun*C.Msun/C.m_p/(args.lpc*C.pc)**3
dens = mtot/((maxD[0]-minD[0])*(maxD[1]-minD[1])*(maxD[2]-minD[2])) * mfac

for i in range(0,3):
    avgP[i] /= mtot
    avgV[i] /= mtot

print()
print("Value          : {:>13s}   {:>13s}   {:>13s}".format('x|u', 'y|v', 'z|w'))
print("-------------- : {:13s}   {:13s}   {:13s}".format('-'*13, '-'*13, '-'*13))
print("Minimum pos    : {:13.6e}   {:13.6e}   {:13.6e}".format(minD[0], minD[1], minD[2]))
print("Maximum pos    : {:13.6e}   {:13.6e}   {:13.6e}".format(maxD[0], maxD[1], maxD[2]))
print("Minimum vel    : {:13.6e}   {:13.6e}   {:13.6e}".format(minD[3], minD[4], minD[5]))
print("Maximum vel    : {:13.6e}   {:13.6e}   {:13.6e}".format(maxD[3], maxD[4], maxD[5]))
print("Average pos    : {:13.6e}   {:13.6e}   {:13.6e}".format(avgP[0], avgP[1], avgP[2]))
print("Average vel    : {:13.6e}   {:13.6e}   {:13.6e}".format(avgV[0], avgV[1], avgV[2]))
print("-------------- : {:13s}   {:13s}   {:13s}".format('-'*13, '-'*13, '-'*13))
print("Density (n/cc) : {:13.6e}".format(dens))
print("Total mass     : {:13.6e}".format(mtot))
print()

