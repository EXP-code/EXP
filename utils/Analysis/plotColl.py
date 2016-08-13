#!/usr/bin/python -i

"""Module with predefined plotting and targets for UserTreeDSMC::CollideIon output

The main functions begin with make* function.  After make* (or a
readDB() which is called by a make*), you may use showlabs() to see
the name tags for all available fields.

   readDB(tag)             : read OUTLOG files and build database

   listAll()               : show the fields available for plotting

"""


# For Python 3 compatibility
from __future__ import absolute_import, division, print_function, unicode_literals

import matplotlib.pyplot as plt
import numpy as np
import os, sys, bisect, re

species  = []
db = {}

def readDB(tag):
    global species, db

    head = 2
    stanza = 16

    file = open(tag + '.ION_coll')
    # Look for species line
    while True:
        line = file.readline()
        species = []
        if re.search('Species==>', line) is not None:
            spc = re.findall('[ |]*(\([0-9]+[ ]*,[ ]*[0-9]+\))[ |]*', line)
            for v in spc: species.append(v.replace(" ", ""))
            break
                
    # Skip the separator line
    line = file.readline()

    # Get all the labels
    line = file.readline()
    labs = re.findall('([a-zA-Z0-9()]+)', line)
            
    # Skip the column count line
    line = file.readline()

    # Skip the separator line
    line = file.readline()

    # Initialize db
    nspc = len(species)
    for i in range(head): db[labs[i]] = []
    for j in range(nspc):
        indx = head + stanza*j
        db[species[j]] = {}
        for k in range(stanza):
            db[species[j]][labs[indx+k]] = []
    indx = head + stanza*nspc
    for i in range(indx,len(labs)): db[labs[i]] = []

    # Process the data
    for line in file:
        toks = re.findall('([+-]*(?:inf|INF|nan|NAN|[+\-0-9.eE]+))', line)
        for i in range(head): db[labs[i]].append(float(toks[i]))
        for j in range(nspc):
            indx = head + stanza*j
            for k in range(stanza):
                db[species[j]][labs[indx+k]].append(float(toks[indx+k]))
        indx = head + stanza*nspc
        for i in range(indx,len(toks)):
            try:
                db[labs[i]].append(float(toks[i]))
            except:
                print("Column={} Variable={} token={}".format(i, labs[i],toks[i]))
                print("Toks=", toks)
                print("Line=", line)
                print("Unexpected error: {}".format(sys.exc_info()[0]))
                raise

    # Convert to numpy arrays
    for k in db:
        if k in species:
            for j in db[k]:
                db[k][j] = np.array(db[k][j])
        else:
            db[k] = np.array(db[k])

def showAll():
    for k in db:
        if k in species:
            print("Species {}:".format(k), end="")
            cnt = 0
            for j in db[k]:
                if cnt % 6 == 0: print("\n     ", end="")
                print("{}".format(j), end=" ")
                cnt += 1
            print()
        else:
            print("Value: {}".format(k))

def listAll():
    for k in db:
        if k in species:
            print("{:10} [species]".format(k))
        else:
            print("{:10} [value]".format(k))

def listGlobal():
    for k in db:
        if k not in species:
            print("* {}".format(k))

def listSpecies():
    if k in species:
        print("* {}".format(k))

def showSpecies(v):
    if v not in species:
        print("No species <{}> found".format(v))
    else:
        cnt = 0
        print("     ", end="")
        for j in db[v]:
            cnt += 1
            print("{}".format(j), end=" ")
            if cnt % 6 == 0: print("\n     ", end="")
        print()
        

def main(argv):
        """ Parse the command line and call the parsing and plotting routine """
        #
        # Last argument should be filename and must exist
        #
        if len(argv) <=1:
                print("Usage: {} runtag".format(argv[0]))
                exit(1)

        readDB(argv[1])
        listAll()

if __name__ == "__main__":
        main(sys.argv)
