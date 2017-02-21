#!/usr/bin/python

# -*- Python -*-
# -*- coding: utf-8 -*-

import sys, os, argparse
import numpy as np
import matplotlib.pyplot as plt
import getSpecies as gs

parser = argparse.ArgumentParser(description='Read DSMC species file and plot temperatures for each species')
parser.add_argument('-t', '--tscale', default=1000.0,    type=float,   help='System time units in years')
parser.add_argument('-T', '--Tmax',   default=1.0e32, type=float,   help='Maximum time in years')
parser.add_argument('tags',           nargs='*',                       help='Files to process')

args = parser.parse_args()

labs = args.tags
        
fields = [ ['Telc(1)', 'Telc(2)'], ['Tion(1)', 'Tion(2)'] ]
titles = ['Electrons', 'Ions']

d = {}
for v in labs:
	d[v] = gs.readDB(v)

f, ax = plt.subplots(2, 1, sharex='col')
for x in ax:
        box = x.get_position()
        newBox = [box.x0, box.y0, 0.9*box.width, box.height]
        x.set_position(newBox)

cnt = 0
for v in labs:
        for n in range(2):
                for F in fields[n]:
                        if F in d[v]:
                                indx = np.searchsorted(d[v]['Time'], args.Tmax/args.tscale)
                                ax[n].plot(d[v]['Time'][0:indx]*args.tscale, d[v][F][0:indx], '-', label=v+':'+F)

                ax[n].legend(prop={'size':8}, bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.0)
                if n>0: ax[n].set_xlabel('Time')
                ax[n].set_ylabel('Temperature')
                ax[n].set_title(titles[n])
plt.show()
