#!/usr/bin/python

file = open("gas.bod")
file.readline()
maxD = [0.0, 0.0, 0.0, 0.0]
for line in file:
    v = line.split()
    maxD[0] += float(v[0])
    for i in range(1,4): maxD[i] = max(maxD[i], float(v[i]))
        
dens = maxD[0]/(maxD[1]*maxD[2]*maxD[3]) * 4.041742
print "Total mass:", maxD[0]
print "Maximum pos:", maxD[1:]
print "Density (n/cc) =", dens

