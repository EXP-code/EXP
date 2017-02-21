#!/usr/bin/python

import sys, os, argparse
import astroval as C

parser = argparse.ArgumentParser(description='Evaluate global statistics for body file created by makeIonIC')
parser.add_argument('-b', '--body',    default='gas.bod', help='Name of body file')
parser.add_argument('-m', '--msun',    default=1.0, help='Number of M_sun per unit mass')
parser.add_argument('-l', '--lpc',     default=1.0, help='Name of parsecs per unit length')

args = parser.parse_args()

file = open(args.body)
file.readline()
maxD = [0.0, 0.0, 0.0, 0.0]
for line in file:
    v = line.split()
    maxD[0] += float(v[0])
    for i in range(1,4): maxD[i] = max(maxD[i], float(v[i]))
        
mfac = args.msun*C.Msun/C.m_p/(args.lpc*C.pc)**3
dens = maxD[0]/(maxD[1]*maxD[2]*maxD[3]) * mfac

print "Total mass:", maxD[0]
print "Maximum pos:", maxD[1:]
print "Density (n/cc) =", dens

