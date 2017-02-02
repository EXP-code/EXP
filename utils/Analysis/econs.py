#!/usr/bin/python

# -*- Python -*-
# -*- coding: utf-8 -*-

import sys, os, argparse
import numpy as np
import matplotlib.pyplot as plt
import string

parser = argparse.ArgumentParser(description='Read ION_Coll file and plot energy disgnostics')
parser.add_argument('-t', '--tscale', default=1000.0,    type=float,   help='System time units in years')
parser.add_argument('-T', '--Tmax',   default=1000000.0, type=float,   help='Maximum time in years')
parser.add_argument('-l', '--log',    default=False, action='store_true', help='Logarithmic vertical scale')
parser.add_argument('-a', '--aux',    default=False, action='store_true', help='Sum energy fields')
parser.add_argument('tags',           nargs='*',                       help='Files to process')

args = parser.parse_args()

labs = args.tags
        
# Translation table to convert vertical bars and comments to spaces
#
trans = string.maketrans("#|", "  ")

d = {}
lead = 2

field = [ ('Time', 0),
          ('Temp', 1),
          ('Disp', 2),
          ('Etotl', -1),
          ('Efrac', -2),
          ('EdspE', -3),
          ('misE', -4),
          ('clrE', -5),
          ('delE', -6),
          ('delI', -7),
          ('PotI', -8),
          ('EkeE', -9),
          ('EkeI', -10),
          ('ElosC', -11),
          ('Elost', -12) ]

f, ax = plt.subplots(1, 1)

box = ax.get_position()
newBox = [box.x0, box.y0, 0.9*box.width, box.height]
ax.set_position(newBox)

toplot = ['Etotl', 'ElosC', 'Elost', 'EkeE', 'EkeI', 'delI', 'delE']

aux = ['ElosC', 'EkeE', 'EkeI', 'delI', 'delE']

for v in labs:
    # Read and parse the file
    #
    file  = open(v + ".ION_coll")
    for line in file:
        if line.find('Time')>=0:    # Get the labels
            next = True
            labels = line.translate(trans).split()
            nlabs  = len(labels)
            tindx  = labels.index('Elost')
            tail   = nlabs - tindx
            if 'Disp' in labels: lead = 3
            ncol   = (tindx-lead)/5
            d[v] = {}
            for e in field: d[v][e[0]] = []
        if line.find('#')<0:        # Read the data lines
            toks = line.translate(trans).split()
            allZ = True             # Skip lines with zeros only
            for i in range(lead,len(toks)):
                if float(toks[i])>0.0: 
                    allZ = False
                    break
            if not allZ:            
                # A non-zero line . . .  Make sure field counts are the
                # same (i.e. guard against the occasional badly written
                # output file
                if len(toks) == len(labels):
                    for e in field:
                        if e[0] in labels:
                            d[v][e[0]].append(float(toks[e[1]]))

    indx = np.searchsorted(d[v]['Time'], args.Tmax/args.tscale)
    x = np.array(d[v]['Time'][0:indx])*args.tscale
            
    for f in toplot:
        if f in d[v]:
            y = np.abs(np.array(d[v][f][0:indx]))
            if args.log:
                ax.semilogy(x, y, '-', label=v+':'+f)
            else:
                ax.plot(x, y, '-', label=v+':'+f)

    if args.aux:
        y = np.copy(x) * 0.0
        for f in aux: y += np.array(d[v][f][0:indx])
        if args.log:
            ax.semilogy(x, y, '-o', label=v+':Esum')
        else:
            ax.plot(x, y, '-o', label=v+':Esum')

ax.legend(prop={'size':8}, bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.0)
ax.set_xlabel('Time')
ax.set_ylabel('Energy')
plt.show()
