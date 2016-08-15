#!/usr/bin/python -i

"""Module with predefined plotting and targets for UserTreeDSMC::CollideIon output

The main functions begin with make* function.  After make* (or a
readDB() which is called by a make*), you may use showlabs() to see
the name tags for all available fields.

   readDB(tag)             : read OUTLOG files and build database

   listAll()               : show the fields available for plotting

   plotEnergy()            : plot global energy quantities

"""

# For Python 3 compatibility
from __future__ import absolute_import, division, print_function, unicode_literals

import matplotlib.pyplot as plt
import numpy as np
import os, sys, re, getopt

species = []
labs    = []
db      = {}

def readDB(tag):
    global species, db, labs

    # <head> leading variables and <stanza> per species variables
    head = 2
    stanza = 16

    # Clear DB
    db = {}

    # Open file
    file = open(tag + '.ION_coll')

    # Look for species line, discarding all the user info
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

    # Initialize db from labels and species info
    nspc = len(species)
    for i in range(head): db[labs[i]] = []
    for j in range(nspc):
        indx = head + stanza*j
        db[species[j]] = {}
        for k in range(stanza):
            db[species[j]][labs[indx+k]] = []
    indx = head + stanza*nspc
    for i in range(indx,len(labs)): db[labs[i]] = []

    # Process the data from the file
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
                print("Trouble reading float value??\nColumn={} Variable={} token={}".format(i, labs[i],toks[i]))
                print("Toks=", toks)
                print("Line=", line)
                print("Unexpected error: {}".format(sys.exc_info()[0]))
                raise

    # Convert lists to numpy arrays
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
        

def plotEnergy(Log=False, lw=2, xtag='Time', scale=1000.0, maxT=1.0e+20, rtag=''):
    """ Plot critical energy conservation quantities """

    if len(rtag)>0: readDB(rtag)

    x = np.copy(db[xtag])
    k = len(x)
    if xtag=='Time':
        x *= scale
        for i in range(len(x)):
            if maxT <= x[i]:
                k = i
                break
    if Log:
        plt.semilogy(x[0:k], db['Etotl'][0:k], '-o', linewidth=lw, label='E total')
        plt.semilogy(x[0:k], db['ElosC'][0:k], '-', linewidth=lw, label='E lost')
        plt.semilogy(x[0:k], db['Elost'][0:k], '-', linewidth=lw, label='dE lost')
        plt.semilogy(x[0:k], db['PotI'][0:k],  '-', linewidth=lw, label='Ion pot')
        plt.semilogy(x[0:k], db['EkeI'][0:k] + db['EkeE'][0:k], '-', linewidth=lw, label='Tot KE')
        plt.semilogy(x[0:k], db['EkeI'][0:k], '-', linewidth=lw, label='Ion KE')
        plt.semilogy(x[0:k], db['EkeE'][0:k], '-', linewidth=lw, label='Elec KE')
        if 'delC' in labs:
            plt.semilogy(x[0:k], db['delC'][0:k], '-', linewidth=lw, label='E excess')
        else:
            plt.semilogy(x[0:k], db['delI'][0:k] + db['delE'][0:k], '-', linewidth=lw, label='E excess')
    else:
        plt.plot(x[0:k], db['Etotl'][0:k], '-o', linewidth=lw, label='E total')
        plt.plot(x[0:k], db['ElosC'][0:k], '-', linewidth=lw, label='E lost')
        plt.plot(x[0:k], db['Elost'][0:k], '-', linewidth=lw, label='dE lost')
        plt.plot(x[0:k], db['PotI'][0:k],  '-', linewidth=lw, label='Ion pot')
        plt.plot(x[0:k], db['EkeI'][0:k] + db['EkeE'][0:k], '-', linewidth=lw, label='Tot KE')
        plt.plot(x[0:k], db['EkeI'][0:k], '-', linewidth=lw, label='Ion KE')
        plt.plot(x[0:k], db['EkeE'][0:k], '-', linewidth=lw, label='Elec KE')
        if 'delC' in labs:
            plt.plot(x[0:k], db['delC'][0:k], '-', linewidth=lw, label='E excess')
        else:
            plt.plot(x[0:k], db['delE'][0:k] + db['delI'][0:k], '-', linewidth=lw, label='E excess')
            plt.plot(x[0:k], db['delI'][0:k], '-', linewidth=lw, label='Ion excess')
            plt.plot(x[0:k], db['delE'][0:k], '-', linewidth=lw, label='Elec excess')

    if xtag=='Time':
        plt.xlabel('Time')
    else:
        plt.xlabel('Temperature')

    plt.ylabel('Energy')
    plt.legend().draggable()
    plt.show()

def main(argv):
    """ Parse the command line and call the parsing and plotting routine """

    helpstring = \
    ' [-t <timescale> | --timescale=<timescale>]' + \
    ' [-T <max time> | --maxT=<max time>]' + \
    ' [--time] [--temp]'
    ' <runtag>'

    energy   = False
    logscale = False
    lw       = 2.0
    xtag     = 'Time'
    tscale   = 1000.0
    maxT     = 1.0e+20

    try:
        opts, args = getopt.getopt(argv, "helw:t:T:",
                                   ["help", "energy", "logscale", "linewidth=", "time", "temp", "timescale=", "maxT="])
    except getopt.GetoptError:
        print(sys.argv[0], helpstring)
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], helpstring)
            sys.exit()
        elif opt in ("-e", "--energy"):
            energy = True
        elif opt in ("-l", "--logscale"):
            logscale = True
        elif opt in ("-w", "--linewidth"):
            lw = float(arg)
        elif opt in ("--time"):
            xtag = 'Time'
        elif opt in ("--temp"):
            xtag = 'Temp'
        elif opt in ("-t", "--timescale"):
            tscale = float(arg)
        elif opt in ("-T", "--maxT"):
            maxT = float(arg)

    #
    # Last argument should be filename and must exist
    #
    if len(args)<=0:
        print("Usage: {} runtag".format(argv[0]))
        exit(1)

    readDB(args[-1])
    if energy:
        plotEnergy(Log=logscale, lw=lw, scale=tscale, maxT=maxT, xtag=xtag)
    else:
        listAll()

if __name__ == "__main__":
        main(sys.argv[1:])
