#!/usr/bin/python

# -*- Python -*-
# -*- coding: utf-8 -*-

import sys, os, argparse
import numpy as np
import matplotlib.pyplot as plt
import getSpecies as gs

parser = argparse.ArgumentParser(description='Read DSMC species file and plot temperatures')
parser.add_argument('-t', '--tscale', default=1000.0,  help='System time units in years')
parser.add_argument('-T', '--Tmax', default=1.0e32, help='Maximum time in years')
parser.add_argument('tags', nargs='*', help='Files to process')

args = parser.parse_args()

labs = args.tags

fields = ['Telc(1)', 'Telc(2)']
fields = ['Telc(1)', 'Telc(2)', 'Tion(1)', 'Tion(2)']

d = {}
for v in labs:
	d[v] = gs.readDB(v)

f, ax = plt.subplots(1, 1)

box = ax.get_position()
newBox = [box.x0, box.y0, 0.9*box.width, box.height]
ax.set_position(newBox)

cnt = 0
for v in labs:
    for f in fields:
        if f in d[v]:
                indx = np.searchsorted(d[v]['Time'], args.Tmax/args.tscale)
                ax.plot(d[v]['Time'][0:indx]*args.tscale, d[v][f][0:indx], '-', label=v+':'+f)

ax.set_xlabel('Time')
ax.set_ylabel('Temperature')
ax.legend(prop={'size':8}, bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.0)
plt.show()
