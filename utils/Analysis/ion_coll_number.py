#!/usr/bin/python

# -*- Python -*-
# -*- coding: utf-8 -*-

"""Program to display the particle collision counts for each type
using diagnostic output from the CollideIon class in the
UserTreeDSMC module.

There are two simple routines here.  The main routine that parses the
input command line and a plotting/parsing routine.

Examples:

	$ python ion_coll_number -d 2 run2

Plots the collision count types for each ion type.

"""

import sys, getopt
import copy
import string
import numpy as np
import matplotlib.pyplot as pl
import scipy.interpolate as ip


def plot_data(filename, msz, line, dot, elastic, useTime, stride, trange):
    """Parse and plot the *.ION_coll output files generated by
    CollideIon

    Parameters:

    filename (string): is the input datafile name

    line (bool): if True, connect the dots

    dot (bool): if True, markers set to dots

    elastic (bool): if True, plot elastic interaction counts

    stride (int): plot every stride points
    """

    # Marker type
    #
    mk = ''
    if line: mk = '-'
    if dot:  mk += '.'
    else:    mk += '*'

    # Translation table to convert vertical bars and comments to spaces
    #
    trans = string.maketrans("#|", "  ")

    # Initialize data and header containers
    #
    tabl  = {}
    time  = []
    temp  = []
    etot  = []
    ncol  = 9
    head  = 2
    tail  = 2
    data  = {}

    # Species
    #
    spec  = ['H', 'H+', 'He', 'He+', 'He++']
    for v in spec: data[v] = {}

    # Read and parse the file
    #
    file  = open(filename)
    for line in file:
        if line.find('Time')>=0:    # Get the labels
            next = True
            labels = line.translate(trans).split()
            if line.find("W(") >= 0: 
                ncol = 13
            if line.find("N(nn)") >= 0: 
                ncol = 16
            if line.find("EratC") >= 0: 
                tail = 12
        if line.find('[1]')>=0:     # Get the column indices
            toks = line.translate(trans).split()
            for i in range(2, len(toks)-tail):
                j = int(toks[i][1:-1]) - 1
                tabl[labels[j]] = i
                idx = (i-head) / ncol
                data[spec[idx]][labels[j]] = []
        if line.find('#')<0:        # Read the data lines
            toks = line.translate(trans).split()
            allZ = True             # Skip lines with zeros only
            for i in range(head, len(toks)):
                if float(toks[i])>0.0: 
                    allZ = False
                    break
            if not allZ:            
                # A non-zero line . . .  Make sure field counts are the
                # same (i.e. guard against the occasional badly written
                # output file
                if len(toks) == len(labels):
                    time.append(float(toks[0]))
                    temp.append(float(toks[1]))
                    etot.append(float(toks[-1]))
                    for i in range(head, len(toks)-tail):
                        idx = (i-head) / ncol
                        data[spec[idx]][labels[i]].append(float(toks[i]))
                else:
                    print "toks=", len(toks), " labels=", len(labels)

    # Fields to plot
    #
    nkeys = 4
    ncols = 2
    if elastic:
        ekeys = ['N(ne)', 'N(ie)', 'N(ce)', 'N(ci)', 'N(ff)', 'N(rr)']
        elabs = ['neut-elec', 'ion-elec', 'collide', 'ionize', 'free-free', 'recomb']
        nkeys = 6
        ncols = 3
    else:
        ekeys = ['N(ce)', 'N(ci)', 'N(ff)', 'N(rr)']
        elabs = ['collide', 'ionize', 'free-free', 'recomb']

    if useTime: xaxis = time
    else:       xaxis = temp

    icnt = 0
    for k in range(0, nkeys):
        icnt += 1
        pl.subplot(ncols, 2, icnt)
        if useTime: pl.xlabel('Time')
        else:       pl.xlabel('Temperature')
        pl.ylabel('Counts')
        for v in spec:
            x = xaxis
            y = data[v][ekeys[k]]
            if stride>1 and stride<len(xaxis)/2:
                x = []
                y = []
                for i in range(0,len(xaxis),stride):
                    if time[i] >= trange[0] and time[i] <= trange[1]:
                        x.append(xaxis[i])
                        y.append(data[v][ekeys[k]][i])
            pl.plot(x, y, mk, label=v, markersize=msz)
        pl.title(elabs[k])
        if icnt==nkeys:
            leg = pl.legend(loc='best',borderpad=0,labelspacing=0)
            leg.get_title().set_fontsize('6')
            pl.setp(pl.gca().get_legend().get_texts(), fontsize='12')

    pl.get_current_fig_manager().full_screen_toggle()
    pl.show()

    # Fields to plot
    #
    ekeys = ['W(ce)', 'W(ci)', 'W(ff)', 'W(rr)']
    elabs = ['collide', 'ionize', 'free-free', 'recomb']

    icnt = 0
    for k in range(0, 4):
        icnt += 1
        pl.subplot(2, 2, icnt)
        if useTime: pl.xlabel('Time')
        else:       pl.xlabel('Temperature')
        pl.ylabel('Weights')
        # pl.ylim((0, 30))
        for v in spec:
            x = xaxis
            y = data[v][ekeys[k]]
            if stride>1 and stride<len(xaxis)/2:
                x = []
                y = []
                for i in range(0,len(xaxis),stride):
                    if time[i] >= trange[0] and time[i] <= trange[1]:
                        x.append(xaxis[i])
                        y.append(data[v][ekeys[k]][i])
            pl.plot(x, y, mk, label=v, markersize=msz)
        pl.title(elabs[k])
        if icnt==4:
            leg = pl.legend(loc='best',borderpad=0,labelspacing=0)
            leg.get_title().set_fontsize('6')
            pl.setp(pl.gca().get_legend().get_texts(), fontsize='12')

    pl.get_current_fig_manager().full_screen_toggle()
    pl.show()

    # Fields to plot
    #
    ekeys = ['W(ce)', 'W(ci)', 'W(ff)', 'W(rr)']
    elabs = ['collide', 'ionize', 'free-free', 'recomb']

    for v in spec:
        pl.subplot(1, 1, 1)
        if useTime: pl.xlabel('Time')
        else:       pl.xlabel('Temperature')
        pl.ylabel('Weights')
        cnt = 0
        for k in range(4):
            x = []
            y = []
            for i in range(0,len(xaxis),stride):
                if time[i] >= trange[0] and time[i] <= trange[1]:
                    if data[v][ekeys[k]][j]>0.0:
                        x.append(xaxis[j])
                        y.append(data[v][ekeys[k]][j])
            if len(x)>0:
                pl.semilogy(x, y, mk, label=elabs[k], markersize=msz)
                cnt += 1
        if cnt==0: continue
        pl.title(v)
        leg = pl.legend(loc='best',borderpad=0,labelspacing=0)
        leg.get_title().set_fontsize('6')
        pl.setp(pl.gca().get_legend().get_texts(), fontsize='12')
        pl.get_current_fig_manager().full_screen_toggle()
        pl.show()


    # Fields to plot
    #
    if elastic:
        ekeys = ['N(nn)', 'N(ne)', 'N(ie)', 'N(ce)', 'N(ci)', 'N(ff)', 'N(rr)']
        elabs = ['neutral', 'neut-elec', 'ion-elec', 'collide', 'ionize', 'free-free', 'recomb']
    else:
        ekeys = ['N(ce)', 'N(ci)', 'N(ff)', 'N(rr)']
        elabs = ['collide', 'ionize', 'free-free', 'recomb']

    for v in spec:
        pl.subplot(1, 1, 1)
        if useTime: pl.xlabel('Time')
        else:       pl.xlabel('Temperature')
        pl.ylabel('Counts')
        cnt = 0
        for k in range(len(ekeys)):
            x = xaxis
            y = data[v][ekeys[k]]
            if stride>1 and stride<len(xaxis)/2:
                x = []
                y = []
                for i in range(0,len(xaxis),stride):
                    if time[i] >= trange[0] and time[i] <= trange[1]:
                        x.append(xaxis[i])
                        y.append(data[v][ekeys[k]][i])
            pl.plot(x, y, mk, label=elabs[k], markersize=msz)
        pl.title(v)
        leg = pl.legend(loc='best',borderpad=0,labelspacing=0)
        leg.get_title().set_fontsize('6')
        pl.setp(pl.gca().get_legend().get_texts(), fontsize='12')
        pl.get_current_fig_manager().full_screen_toggle()
        pl.show()

    # Fields to plot
    #
    ekeys = ['E(ce)', 'E(ci)', 'E(ff)', 'E(rr)']
    elabs = ['collide', 'ionize', 'free-free', 'recomb']

    for v in spec:
        pl.subplot(1, 1, 1)
        if useTime: pl.xlabel('Time')
        else:       pl.xlabel('Temperature')
        pl.ylabel('Energy')
        cnt = 0
        for k in range(4):
            x = xaxis
            y = data[v][ekeys[k]]
            if stride>1 and stride<len(xaxis)/2:
                x = []
                y = []
                for i in range(0,len(xaxis),stride):
                    if time[i] >= trange[0] and time[i] <= trange[1]:
                        x.append(xaxis[i])
                        y.append(data[v][ekeys[k]][i])
            pl.plot(x, y, mk, label=elabs[k], markersize=msz)
        pl.title(v)
        leg = pl.legend(loc='best',borderpad=0,labelspacing=0)
        leg.get_title().set_fontsize('6')
        pl.setp(pl.gca().get_legend().get_texts(), fontsize='12')
        pl.get_current_fig_manager().full_screen_toggle()
        pl.show()


def main(argv):
    """ Parse the command line and call the parsing and plotting routine """

    elastic = False
    line = False
    dot  = False
    useTime = False
    stride = 1
    msz = 4
    trange = [0.0, 1.0e20]

    useString = '[-p | --point | -l | --line | -t | --time | -e | --elastic | -m <size> | --msize=<size> | -s <stride> | --stride=<stride> | --tmin=<time> | --tmax=<time>] <runtag>'

    try:
        opts, args = getopt.getopt(argv,"hm:plets:", ["msize=", "point", "line", "elastic", "time", "stride=", "tmin=", "tmax="])
    except getopt.GetoptError:
        print 'Syntax Error'
        print sys.argv[0], usageString
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print sys.argv[0], usateString
            sys.exit()
        elif opt in ("-l", "--line"):
            line = True
        elif opt in ("-p", "--point"):
            dot = True
        elif opt in ("-e", "--elastic"):
            elastic = True
        elif opt in ("-m", "--msize"):
            msz = int(arg)
        elif opt in ("-t", "--time"):
            useTime = True
        elif opt in ("-s", "--stride"):
            stride = int(arg)
        elif opt in ("--tmin"):
            trange[0] = float(arg)
        elif opt in ("--tmax"):
            trange[1] = float(arg)

    suffix = ".ION_coll"
    if len(args)>0:
        filename = args[0] + suffix;
    else:
        filename = "run" + suffix;

    plot_data(filename, msz, line, dot, elastic, useTime, stride, trange)

if __name__ == "__main__":
   main(sys.argv[1:])
