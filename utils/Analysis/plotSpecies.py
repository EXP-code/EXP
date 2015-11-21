#!/usr/bin/python -i

"""Module with predefined plotting and targets for UserTreeDSMC::CollideIon output

The main functions begin with make* function.  After make* (or a
readDB() which is called by a make*), you may use showlabs() to see
the name tags for all available fields.

   setTag(tags)            : set run tag or list of run tags (default: "run")

   readDB()                : read OUTLOG files and build database

   showLabs()              : show the fields available for plotting

   pview(xl, labs)         : plot field(s) in "labs" against field "xl"

   pview(xl, labs, True)   : plot sum of list "labs" against field "xl"

   xl may be label or a 3-tuple (label, min_value, max_value)

"""

import matplotlib.pyplot as plt
import numpy as np
import os, sys, bisect

tags     = ["run"]
flab     = []

H_spc    = ['(1,1)', '(1,2)']
He_spc   = ['(2,1)', '(2,2)', '(2,3)']
temps    = ['Tion(1)','Tion(2)','Telc(1)','Telc(2)']
energies = ['Eion(1)','Eion(2)','Eelc(1)','Eelc(2)']
E_cons   = ['Cons_E', 'Cons_G', 'Ions_E', 'Elec_E', 'Totl_E']
E_sum    = ['Ions_E', 'Elec_E']

def setTag(x):
    """ Set desired *.species file tag for plotting"""
    global tags
    
    if isinstance(x, str):
        tags = [x]
    elif isinstance(x, list):
        tags = x
    else:
        print "Parameter must be a string or a list of strings"
        return None

def readDB():
    global tags, flab

    db   = {}
    flab = []

    for tag in tags:
        file = open(tag + ".species")
        data = []
        line = file.readline()  # Label line
        labs = [v for v in line[1:].split()]
        for line in file:
            if line.find('#')<0:
                data.append([float(v) for v in line.split()])

        a = np.array(data).transpose()

        db[tag] = {}

        for i in range(a.shape[0]):
            db[tag][labs[i]] = a[i]

        for v in labs:
            if v not in flab: flab.append(v)

    return db


def showLabs():
    icnt = 0
    for v in flab:
        icnt += 1
        print "{:10s}".format(v),
        if icnt % 6 == 0: print


def pview(xl, x, do_sum=False):
    if isinstance(x, str):
        return makeP(xl, x)
    elif isinstance(x, list):
        if do_sum:
            return makeS(xl, x)
        else:
            return makeM(xl, x)
    else:
        return None


def makeS(xL, labs):

    db = readDB()
    xl = []

    bound = False
    if isinstance(xL, tuple):
        xl   = xL[0]
        minv = xL[1]
        maxv = xL[2]
        bound = True
    else:
        xl = xL

    for lab in labs+[xl]:
        if lab not in flab:
            print "No such field, available data is:"
            showLabs()
            return

    # Set the figure
    fig = plt.figure()
    ax  = plt.subplot(111)

    # Shrink current axis's width by 25%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])



    for t in tags:
        sum = np.zeros(db[t][labs[0]].shape)
        for lab in labs: sum+= db[t][lab]
        if bound:
            imin = bisect.bisect_left (db[t][xl], minv)
            imax = bisect.bisect_right(db[t][xl], maxv)
            ax.plot(db[t][xl][imin:imax], db[t][lab][imin:imax],
                     '-', label=t+":sum")
        else:
            ax.plot(db[t][xl], db[t][lab], '-', label=t+":sum")

    # Put a legend to right of plot
    legend = ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1.0),
                       fancybox=True, shadow=True, ncol=1)

    for label in legend.get_lines():
        label.set_linewidth(2.0)  # the legend line width

    plt.xlabel(xl)
    plt.ylabel("Sum")
    plt.show()

def makeM(xL, labs):

    db = readDB()
    xl = []
    
    bound = False
    if isinstance(xL, tuple):
        xl   = xL[0]
        minv = xL[1]
        maxv = xL[2]
        bound = True
    else:
        xl = xL

    for lab in labs+[xl]:
        if lab not in flab:
            print "No such field, available data is:"
            showLabs()
            return

    # Set the figure
    fig = plt.figure()
    ax  = plt.subplot(111)

    # Shrink current axis's width by 25%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0,
                     box.width * 0.75, box.height])

    for t in tags:
        for lab in labs:
            l = t + ":" + lab
            if bound:
                imin = bisect.bisect_left (db[t][xl], minv)
                imax = bisect.bisect_right(db[t][xl], maxv)
                ax.plot(db[t][xl][imin:imax], db[t][lab][imin:imax],
                         '-', label=l)
            else:
                ax.plot(db[t][xl], db[t][lab], '-', label=l)

    # Put a legend to right of plot
    legend = ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1.0),
                       fancybox=True, shadow=True, ncol=1)

    for label in legend.get_lines():
        label.set_linewidth(2.0)  # the legend line width

    ax.set_xlabel(xl)
    ax.set_ylabel(lab)
    plt.show()

def makeP(xl, lab):
    db = readDB()

    for l in [xl, lab]:
        if l not in flab:
            print "No such field, available data is:"
            showLabs()
            return
    
    for t in tags:
        plt.plot(db[t][xl], db[t][lab], '-', label=t)
        
    plt.legend()
    plt.xlabel(xl)
    plt.ylabel(lab)
    plt.show()
