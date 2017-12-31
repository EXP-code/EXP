#!/usr/bin/python -i

# For Python 3 compatibility
from __future__ import absolute_import, division, print_function, unicode_literals

import re, sys
import numpy as np
import numpy as np

import matplotlib
matplotlib.use('GTKAgg')
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter

data = {}
labs = []
rtag = sys.argv[-1]
Tscale = 1e3

def readem(Tag=''):
    global data, labs, rtag

    if Tag != '': rtag = Tag

    file = open('{}.ION_coll'.format(rtag), 'r')
    data = {}
    labs = []
    for line in file:
        if line[0] == '#':
            if 'Time' in line:
                line = re.sub('[#| ]+', ' ', line)
                toks = line.split()
                for tok in toks:
                    labs.append(tok)
                    data[tok] = []
        else:
            line = re.sub('[| ]+', ' ', line)
            toks = line.split()
            for v in zip(labs, toks): data[v[0]].append(float(v[1]))

def plotem(l, Log=False, Smooth=False, Tag='', Intvl=10):
    global rtag

    if Tag != '': rtag = Tag
    readem()

    t = np.array(data['Time'])*Tscale
    for v in l:
        y = np.array(data[v])
        if Log:
            plt.semilogy(t, y, '-', label=v)
        else:
            plt.plot(t, y, '-', label=v)

        if Smooth:
            if Log: y = np.log(y)
            # print("tsize=", t.shape[0], " ysize=", y.shape[0])
            itp = interp1d(t, y, kind='linear')
            wsize = 11
            dsize = int(y.shape[0]/Intvl)
            if dsize > wsize:
                wsize = dsize
                if 2*int(wsize/2) == wsize: wsize += 1
            porder = 3
            # print("wsize=", wsize)
            yy = savgol_filter(itp(t), wsize, porder)
            if Log:
                plt.semilogy(t, np.exp(yy), '-', linewidth=2)
            else:
                plt.plot(t, yy, '-', linewidth=2)

    plt.xlabel('Time (years)')
    plt.ylabel(v)
    plt.grid()
    plt.legend().draggable()
    plt.show()
    
def plotLoss(Log=True, Smooth=True, Tag='', Intvl=10):
    plotem(['Elost'], Log=Log, Smooth=Smooth, Tag=Tag, Intvl=Intvl)
    
def plotKE(Log=False, Smooth=True, Tag='', Intvl=10):
    plotem(['EkeI', 'EkeE'], Log=Log, Smooth=Smooth, Tag=Tag, Intvl=Intvl)
    
