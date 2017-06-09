#!/usr/bin/python

# For Python 3 compatibility
from __future__ import absolute_import, division, print_function, unicode_literals

import sys, os, argparse
import math
import numpy as np
import astroval as C

parser = argparse.ArgumentParser(description='Evaluate global statistics for body and spec file created by makeIonIC\nwith specific hybrid-method values')
parser.add_argument('-b', '--body',    default='gas.bod', help='Name of body file')
parser.add_argument('-s', '--spec',    default='species.spec', help='Name of species map file')
parser.add_argument('-m', '--msun',    default=1.0, type=float, help='Number of M_sun per unit mass')
parser.add_argument('-l', '--lpc',     default=1.0, type=float, help='Number of parsecs per unit length')

args = parser.parse_args()

file = open(args.spec)
line = file.readline()
if line[0:6] != 'hybrid':
    print("This only makes sense for a hybrid IC file")
    exit(1)

line  = file.readline()
hpos  = int(line.split()[1])

file  = open(args.body)
line  = file.readline()
toks  = line.split()
numb  = int(toks[0])
inatr = int(toks[1])
dnatr = int(toks[2])

bpos  = 7 + inatr + hpos
maxD  = [0.0, 0.0, 0.0, 0.0]
masS  = {}
tots  = 0
engI  = {}
engE  = {}
spcS  = {}

atomic_masses = [0.000548579909, 1.00794, 4.002602]

for line in file:
    v = line.split()
    mass = float(v[0])
    maxD[0] += mass
    ityp = int(v[7])

    v2E  = 0.0
    masE = mass * atomic_masses[0]/atomic_masses[ityp]
    if ityp not in engI: engI[ityp] = 0.0
    if ityp not in engE: engE[ityp] = 0.0

    for i in range(3):
        vI = float(v[4+i])
        vE = float(v[-4+i])
        engI[ityp] += 0.5 * mass * vI*vI
        engE[ityp] += 0.5 * masE * vE*vE

    if ityp not in masS: masS[ityp] = (0, 0.0)
    masS[ityp] = (masS[ityp][0] + 1, masS[ityp][1] + mass)
    for i in range(1,4): maxD[i] = max(maxD[i], float(v[i]))
    sum = 0.0
    for i in range(0,3): sum += float(v[bpos+i])
    if math.fabs(sum - 1.0) > 1.0e-6: tots += 1
    if ityp not in spcS: spcS[ityp] = np.zeros(ityp+1)
    for i in range(0,ityp+1): spcS[ityp][i] += float(v[bpos+i])
        
mfac = args.msun*C.Msun/C.m_p/(args.lpc*C.pc)**3
dens = maxD[0]/(maxD[1]*maxD[2]*maxD[3]) * mfac

print('-'*60)
print("Total mass:", maxD[0])
print('-'*60)
print("Subspecies tally")
print('-'*60)
print("{:>4s}: {:>8s} {:>13s} {:>13s}".format("Sp", "Count", "Mass", "Mass frac"))
print("{:>4s}: {:>8s} {:>13s} {:>13s}".format("--", "-----", "----", "---------"))

for k in masS:
    print("{:-4d}: {:8d} {:13.6e} {:13.6e}".format(k, masS[k][0], masS[k][1], masS[k][1]/maxD[0]))

print('-'*60)
print('Subspecies states')
print('-'*60)
print("{:>4s}:".format("Sp"), end='')
zMax = max(spcS.keys())
for z in range(1,zMax+2): print(' {:>13s}'.format('C={}'.format(z)), end='')
print()
print('{:>4s}:'.format("--"), end='')
zMax = max(spcS.keys())
for z in range(1,zMax+2): print(' {:>13s}'.format('-----'), end='')
print()
for k in spcS:
    print('{:-4d}:'.format(k), end='')
    sum = 0.0
    for v in spcS[k]: sum += v
    for v in spcS[k]: print(' {:13.6e}'.format(v/sum), end='')
    print()
print('-'*60)
print('Subspecies energies')
print('-'*60)
print('{:>4s}: {:>13s} {:>13s}'.format("Sp", "Ion", "Elec"))
print('{:>4s}: {:>13s} {:>13s}'.format("--", "--------", "--------"))
sumI = 0.0
sumE = 0.0
for k in engI:
    print('{:>4d}: {:>13.6e} {:>13.6e}'.format(k, engI[k], engE[k]))
    sumI += engI[k]
    sumE += engE[k]
print('{:>4s}: {:>13.6e} {:>13.6e}'.format("Tot", sumI, sumE))
print('-'*60)
tmpl0 = "{:<18s} = [{:9.6f} {:9.6f} {:9.6f}]"
tmpl1 = "{:<18s} = {}"
tmpl2 = "{:<18s} = {}/{}"
print(tmpl0.format("Maximum pos", maxD[1], maxD[2], maxD[3]))
print(tmpl1.format("Density (n/cc)", dens))
print(tmpl2.format("Out of bounds", tots, numb))
print('-'*60)
