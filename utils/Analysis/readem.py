#! /usr/bin/python -i

"""
This is a test
"""

import os, sys
import numpy as np
import matplotlib.pyplot as plt

data   = {}
fields = []
fname  = ""

def readem(*args):
    global fname, data, fields
    if len(args)==1:
        fname  = args[0] + '.ION_coll'

    for line in open(fname, 'r'):
        if line.find('Time') >= 0:
            line = line[1:].replace('|', '')
            fields = line.split()
            for v in fields:
                data[v] = []
        elif line[0] != '#':
            line = line[1:].replace('|', '')
            vals = [float(v) for v in line.split()]
            for i in range(len(fields)):
                data[fields[i]].append(vals[i])

    for f in data:
        data[f] = np.array(data[f])

def slab():
    global fields
    for i in range(len(fields)):
        print '{:5d} {}'.format(i+1, fields[i])
    
def plotem(*fields, **kwargs):
    global data
    for f in fields:
        if f in data:
            plt.plot(data['Time'], data[f], '-', label=f)
        else:
            if f not in ['deltaKE', 'KE']:
                print "No field <{}> in data".format(f)
        
    if 'deltaKE' in kwargs:
        y = data['EkeI'] + data['EkeE']
        plt.plot(data['Time'], y-y[0], '-', label='Delta KE')

    if 'KE' in kwargs:
        y = data['EkeI'] + data['EkeE']
        plt.plot(data['Time'], y, '-', label='KE')

    plt.xlabel('Time')
    plt.legend().draggable()
    plt.show()

readem(sys.argv[1])
