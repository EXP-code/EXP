#!/usr/bin/python

# -*- Python -*-
# -*- coding: utf-8 -*-

import sys, os, argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Read DSMC species file and plot temperatures for each species')
parser.add_argument('-t', '--tscale', default=1000.0, type=float,   help='System time units in years')
parser.add_argument('-T', '--Tmax',   default=1.0e32, type=float,   help='Maximum time in years')
parser.add_argument('-w', '--lwidth', default=2.0,    type=float,   help='linewidth for curves')
parser.add_argument('-l', '--log',    action="store_true",          help='log vertical scale')
parser.add_argument('tags',           nargs='*',                    help='Files to process')

args = parser.parse_args()

labs = args.tags
        
lw   = args.lwidth

def readDB(tag):
        data = {}
        if not os.path.exists(tag+'.species'):
                print "Path <{}> does not exist".format(tag+'.species')
                exit(-1)
        file = open(tag+'.species', 'r')
        labs = file.readline()[1:].split() # Labels
        line = file.readline()             # Indices
        line = file.readline()             # Separator
        for v in labs: data[v] = []
        for line in file:       # Parse the remaining data
                for t in list(enumerate([float(v) for v in line.split()])):
                        data[labs[t[0]]].append(t[1])
        return data

fields = [ ['Temp_i', 'Temp_e'],
           ['(1,1)', '(1,2)', '(2,1)', '(2,2)', '(2,3)'] ]

d = {}
for v in labs:
	d[v] = readDB(v)

f, ax = plt.subplots(2, 1, sharex='col')
plots  = [0, 1]

for x in ax:
        box = x.get_position()
        newBox = [box.x0, box.y0, 0.9*box.width, box.height]
        x.set_position(newBox)

for v in labs:
        indx = np.searchsorted(d[v]['Time'], args.Tmax/args.tscale)
        for n in plots:
                for F in fields[n]:
                        if F in d[v]:
                                x = np.array(d[v]['Time'][0:indx])*args.tscale
                                y = np.array(d[v][F][0:indx])
                                if args.log:
                                        ax[n].semilogy(x, y, '-', linewidth=lw, label=v+':'+F)
                                else:
                                        ax[n].plot(x, y, '-', linewidth=lw, label=v+':'+F)

for n in plots:
        ax[n].legend(prop={'size':8}, bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.0)

ax[0].set_ylabel('Temperature')
ax[1].set_xlabel('Time')
ax[1].set_ylabel('Species density')

plt.show()
