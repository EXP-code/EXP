#!/usr/bin/python

import sys
import math
import numpy as np

fname = "gas.bod"
argc  = len(sys.argv)
potf  = 7
T     = 1.0e+05
L     = 1.0
M     = 0.1

for i in range(len(sys.argv)):
    if sys.argv[i] == '-f' or sys.argv[i] == '--file':
        fname = sys.argv[i+1]
        ptof  = 8
    if sys.argv[i] == '-T' or sys.argv[i] == '--time':
        T = float(sys.argv[i+1])
    if sys.argv[i] == '-M' or sys.argv[i] == '--mass':
        M = float(sys.argv[i+1])
    if sys.argv[i] == '-L' or sys.argv[i] == '--length':
        L = float(sys.argv[i+1])

Tunit = 31556925.19 * T
Lunit = 3.08568025e+18 * L
Munit = 1.98892e+33 * M
Dunit = Munit/Lunit**3/1.674927211e-24

file = open("species.spec")
line = file.readline()
if line[0:6] != 'hybrid':
    print "This only makes sense for a hybrid IC file"
    exit(1)

line  = file.readline()
hpos  = int(line.split()[1])

file  = open(fname)
line  = file.readline()
toks  = line.split()
numb  = int(toks[0])
inatr = int(toks[1])
dnatr = int(toks[2])

bpos  = potf + inatr + hpos
maxD  = [0.0, 0.0, 0.0, 0.0]
masS  = {}
tots  = 0

spcS  = {}

for line in file:
    v = line.split()
    mass = float(v[0])
    maxD[0] += mass
    ityp = int(v[potf])
    if ityp not in masS: masS[ityp] = (0, 0.0)
    masS[ityp] = (masS[ityp][0] + 1, masS[ityp][1] + mass)
    for i in range(1,4): maxD[i] = max(maxD[i], float(v[i]))
    sum = 0.0
    for i in range(0,3): sum += float(v[bpos+i])
    if math.fabs(sum - 1.0) > 1.0e-6: tots += 1
    if ityp not in spcS: spcS[ityp] = np.zeros(ityp+1)
    for i in range(0,ityp+1): spcS[ityp][i] += float(v[bpos+i])
        
dens = maxD[0]/(maxD[1]*maxD[2]*maxD[3]) * Dunit
print '-'*60
print "Total mass:", maxD[0]
print '-'*60
print "Subspecies tally"
print '-'*60
print "{:>4s}: {:>8s} {:>13s} {:>13s}".format("Sp", "Count", "Mass", "Mass frac")
print "{:>4s}: {:>8s} {:>13s} {:>13s}".format("--", "-----", "----", "---------")

for k in masS:
    print "{:-4d}: {:8d} {:13.6e} {:13.6e}".format(k, masS[k][0], masS[k][1], masS[k][1]/maxD[0])

print '-'*60
print "Subspecies states"
print '-'*60
print "{:>4s}:".format("Sp"),
zMax = max(spcS.keys())
for z in range(1,zMax+2): print " {:>13s}".format("C={}".format(z)),
print ""
print "{:>4s}:".format("--"),
zMax = max(spcS.keys())
for z in range(1,zMax+2): print " {:>13s}".format("-----"),
print ""
for k in spcS:
    print "{:-4d}:".format(k),
    sum = 0.0
    for v in spcS[k]: sum += v
    for v in spcS[k]: print " {:13.6e}".format(v/sum),
    print ""
print '-'*60
print "Maximum pos:", maxD[1:]
print "Density (n/cc) =", dens
print "Out of bounds =", tots, "/", numb
print '-'*60
