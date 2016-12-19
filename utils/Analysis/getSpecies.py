#!/usr/bin/python

"""Module with predefined plotting and targets for UserTreeDSMC::CollideIon output

The main functions begin with make* function.  After make* (or a
readDB() which is called by a make*), you may use showlabs() to see
the name tags for all available fields.

   readDB(tag)             : read *.species files and build database

   listAll()               : show the fields available for plotting

"""

# For Python 3 compatibility
from __future__ import absolute_import, division, print_function, unicode_literals

import matplotlib.pyplot as plt
import numpy as np
import os, sys, re, getopt

def readDB(tag):
    global species, db, labs

    # <head> leading variables and <stanza> per species variables
    head = 2
    stanza = 16

    # Clear DB
    db = {}

    # Labels
    labs = []

    # Check file
    filename = tag + '.species'
    if not os.path.exists(filename):
        raise NameError('Path <{}> does not exist'.format(filename))

    # Open file
    file = open(tag + '.species')

    # Get all the labels
    line = file.readline()
    labs = re.findall('([a-zA-Z0-9()\,\_\[\]]+)', line)
            
    # Skip the separator line
    line = file.readline()

    # Number of fields for checking integrity of data line
    #
    nlab = len(labs)

    # Initialize db from labels and species info
    db = {}
    for v in labs: db[v] = []

    # Process the data from the file
    #
    for line in file:
        toks = re.findall('([\+\-]*(?:inf|INF|nan|NAN|[\+\-0-9.eE]+))', line)
        nvec = len(toks)
        # Check number of fields with expected
        #
        if nvec == nlab:
            for i in range(nvec):
                try:
                    db[labs[i]].append(float(toks[i]))
                except:
                    print("Trouble reading float value??")
                    print("Toks[{}]={}".format(len(toks), toks))
                    print("Line=", line)
                    print("labels[{}]={}".format(len(labs), labs))
                    print("Attempting to read column {} out of {} expect {}".format(i, len(toks), len(labs)))

    # Convert lists to numpy arrays
    #
    for k in db:
        db[k] = np.array(db[k])

    return db
