#!/usr/bin/python

# -*- Python -*-
# -*- coding: utf-8 -*-

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import getSpecies as gs

def is_number(s):
        try:
                float(s)
                return True
        except ValueError:
                return False

tmax   = 1000000.0
if len(sys.argv)==1:
        print "Usage: ", sys.argv[0], " tag1 tag2 . . . [Tmax]"
        os.exit(-1)

labs = []
if is_number(sys.argv[-1]):
        tmax = float(sys.argv[-1])
        labs = sys.argv[1:-1]
else:
        labs = sys.argv[1:]
        
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
                indx = np.searchsorted(d[v]['Time'], tmax)
                ax.plot(d[v]['Time'][0:indx], d[v][f][0:indx], '-', label=v+':'+f)

ax.set_xlabel('Time')
ax.set_ylabel('Temperature')
ax.legend(prop={'size':8}, bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.0)
plt.show()
