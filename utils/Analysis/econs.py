#!/usr/bin/python

# -*- Python -*-
# -*- coding: utf-8 -*-

import sys, os, argparse
import numpy as np
import matplotlib.pyplot as plt
import string

parser = argparse.ArgumentParser(description='Read ION_Coll file and plot energy disgnostics')
parser.add_argument('-t', '--tscale', default=1000.0,    type=float,   help='System time units in years')
parser.add_argument('-T', '--Tmax',    default=1.0e32, type=float,         help='Maximum time in years')
parser.add_argument('-l', '--log',     default=False, action='store_true', help='Logarithmic vertical scale')
parser.add_argument('-a', '--aux',     default=False, action='store_true', help='Sum energy fields')
parser.add_argument('-k', '--ke',      default=False, action='store_true', help='Total kinetic energy')
parser.add_argument('-d', '--delta',   default=False, action='store_true', help='Plot fraction of deferred energy to total')
parser.add_argument('-c', '--compare', default=False, action='store_true', help='Total energy minus kinetic energy')
parser.add_argument('-b', '--both',    default=False, action='store_true', help='Plot KE and Total E separately')
parser.add_argument('-w', '--lwidth',  default=1.0, type=float,            help='linewidth for curves')
parser.add_argument('tags',            nargs='*',                          help='Files to process')

args = parser.parse_args()

labs = args.tags
        
lw   = args.lwidth

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
          ('KlosC', -11),
          ('Klost', -12),
          ('ElosC', -13),
          ('Elost', -14) ]

f, ax = plt.subplots(1, 1)

box = ax.get_position()
newBox = [box.x0, box.y0, 0.9*box.width, box.height]
ax.set_position(newBox)

toplot = ['Etotl', 'ElosC', 'KlosC', 'Elost', 'Klost', 'EkeE', 'EkeI', 'delI', 'delE']

kesum = ['EkeE', 'EkeI']

delE = ['delE', 'delI', 'clrE']

# aux  = {'ElosC':1, 'KlosC':1, 'EkeE':1, 'EkeI':1, 'delI':-1, 'delE':-1}
aux  = {'ElosC':1, 'EkeE':1, 'EkeI':1, 'delI':-1, 'delE':-1}

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
            
    if args.both:
        y = np.copy(x) * 0.0
        yt = np.array(d[v]['Etotl'][0:indx])
        for f in kesum: y += np.array(d[v][f][0:indx])
        if args.log:
            ax.semilogy(x, y,  '-', linewidth=lw, label=v+':KE')
            ax.semilogy(x, yt, '-', linewidth=lw, label=v+':Total E')
        else:
            ax.plot(x, y,  '-', linewidth=lw, label=v+':KE')
            ax.plot(x, yt, '-', linewidth=lw, label=v+':Total E')

    elif args.compare:
        y = np.copy(x) * 0.0
        yt = np.array(d[v]['Etotl'][0:indx])
        for f in kesum: y += np.array(d[v][f][0:indx])
        if args.log:
            ax.semilogy(x, (yt - y)/yt, '-', linewidth=lw, label=v+':Delta E/E')
        else:
            ax.plot(x, (yt - y)/yt, '-', linewidth=lw, label=v+':Delta E/E')

    elif args.delta:
        y = np.copy(x) * 0.0
        denom = np.array(d[v]['Etotl'][0:indx])
        for f in delE:
            y = np.array(d[v][f][0:indx])/denom
            if args.log:
                ax.semilogy(x, y, '-', linewidth=lw, label=v+':'+f)
            else:
                ax.plot(x, y, '-', linewidth=lw, label=v+':'+f)
    else:

        for f in toplot:
            if f in d[v]:
                y = np.abs(np.array(d[v][f][0:indx]))
                if args.log:
                    ax.semilogy(x, y, '-', linewidth=lw, label=v+':'+f)
                else:
                    ax.plot(x, y, '-', linewidth=lw, label=v+':'+f)

        if args.ke:
            y = np.copy(x) * 0.0
            for f in kesum: y += np.array(d[v][f][0:indx])
            if args.log:
                ax.semilogy(x, y, 'o', linewidth=lw, label=v+':KEsum')
            else:
                ax.plot(x, y, 'o', linewidth=lw, label=v+':KEsum')

        if args.aux:
            y = np.copy(x) * 0.0
            for f in aux.keys(): y += np.array(d[v][f][0:indx])*aux[f]
            if args.log:
                ax.semilogy(x, y, '-o', linewidth=lw, label=v+':Esum')
            else:
                ax.plot(x, y, '-o', linewidth=lw, label=v+':Esum')

leg = ax.legend(prop={'size':8}, bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.0)
# set the linewidth of each legend object
for legobj in leg.legendHandles:
    legobj.set_linewidth(2.0)
ax.set_xlabel('Time')
ax.set_ylabel('Energy')
plt.show()
