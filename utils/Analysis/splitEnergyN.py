#!/usr/bin/python

# -*- Python -*-
# -*- coding: utf-8 -*-

import sys, os, argparse
import numpy as np
import matplotlib.pyplot as plt
import getSpecies as gs

parser = argparse.ArgumentParser(description='Read DSMC species file and plot energies per particle for each species')
parser.add_argument('-t', '--tscale', default=1000.0,  type=float, help='System time units in years')
parser.add_argument('-T', '--Tmax', default=1000000.0, type=float, help='Maximum time in years')
parser.add_argument('tags', nargs='*', help='Files to process')

args = parser.parse_args()

labs = args.tags
        
fields = [ ['Eelc(1)', 'Eion(1)'], \
           ['Eelc(2)', 'Eion(2)'] ]
nfield = [ ['Nelc(1)', 'Nion(1)'], \
           ['Nelc(2)', 'Nion(2)'] ]
labels = [ 'Hydrogen', 'Helium' ]

atomic_masses = [0.000548579909, 1.00794, 4.002602]

d = {}
for v in labs:
	d[v] = gs.readDB(v)

f, ax = plt.subplots(2, 1, sharex='col')
for x in ax:
    box = x.get_position()
    newBox = [box.x0+0.05*box.width, box.y0, 0.85*box.width, box.height]
    x.set_position(newBox)

for i in range(2):
    for j in range(2):
        for v in labs:
            f = fields[i][j]
            b = nfield[i][j]
            if f in d[v]:
                indx = np.searchsorted(d[v]['Time'], args.Tmax/args.tscale)
                fv = d[v][f][0:indx]
                bf = d[v][b][0:indx]
                if f.find('elc')>=0:
                    bf /= atomic_masses[0]
                else:
                    if f.find('1')>=0:
                        bf /= atomic_masses[1]
                    if f.find('2')>=0:
                        bf /= atomic_masses[2]
                ax[i].plot(d[v]['Time'][0:indx]*args.tscale, fv/bf, '-', label=v+':'+f)
    if i>0: ax[i].set_xlabel('Time')
    ax[i].set_ylabel('Energy')
    ax[i].set_title(labels[i])
    ax[i].legend(prop={'size':8}, bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.0)
plt.show()
