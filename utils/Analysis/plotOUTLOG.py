#!/usr/bin/python -i

"""Module with predefined plotting and targets for EXP output

The main functions begin with make* function.  After make* (or a
readDB() which is called by a make*), you may use showlabs() to see
the name tags for all available fields.

   setTags(list)           : set list of run tags (default: "run")

   readDB()                : read OUTLOG files and build database

   showLabs()              : show the fields available for plotting

   makeM(xl, labs)         : plot fields in list "labs" against field "xl"

   makeP(xl, lab)          : plot field "lab" against field "xl"

   timeStep()              : plot time-step CPU time against time
                             "file" assuming temperature "temp"

   energies()              : plot energy equilibrium values

   pos()                   : plot center-of-mass position values against time

   vel()                   : plot center-of-mass velocity values against time

   angmom()                : plot total angular momentum values against time

"""


import matplotlib.pyplot as plt
import string
import numpy as np
import os, sys

tags = ["run"]
flab = []

def setTags(intags):
    """ Set desired OUTLOG file tags for plotting"""
    tags = intags

def readDB():
    global flab

    db   = {}
    flab = []

    for tag in tags:
        file = open("OUTLOG." + tag)
        data = []
        line = file.readline()  # Label line
        line = file.readline()  # Separator
        line = file.readline()  # Labels
        labs = [v.strip() for v in line[1:].split('|')]
        line = file.readline()  # Separator
        line = file.readline()  # Indices
        line = file.readline()  # Separator
        for line in file:
            data.append([float(v) for v in line.split('|')])

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


def makeM(xl, labs):
    db = readDB()

    for lab in labs+[xl]:
        if lab not in flab:
            print "No such field, available data is:"
            showLabs()
            return

    for t in tags:
        for lab in labs:
            l = t + ":" + lab
            plt.plot(db[t][xl], db[t][lab], '-', label=l)
    plt.legend()
    plt.xlabel(xl)
    plt.ylabel(lab)
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

def timeStep():
    makeP('Time', 'Clock')

def equil(comp=''):
    s = ''
    if len(c)>0: s = comp + ' '
    makeP('Time', s+' 2T/VC')

def energies(comp=''):
    s = ''
    if len(comp)>0: s = comp + ' '
    makeM('Time', [s+'KE', s+'PE', s+'VC', s+'E'])

def pos(comp=''):
    s = ''
    if len(comp)>0: s = comp + ' '
    makeM('Time', [s+'R(x)', s+'R(y)', s+'R(z)'])
        
def vel(comp=''):
    s = ''
    if len(comp)>0: s = comp + ' '
    makeM('Time', [s+'V(x)', s+'V(y)', s+'V(z)'])
        
def angmom(comp=''):
    s = ''
    if len(comp)>0: s = comp + ' '
    makeM('Time', [s+'L(x)', s+'L(y)', s+'L(z)'])
        
