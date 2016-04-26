#!/usr/bin/python -i

"""Plot cumulative value of a PSP attribute column

"""

import matplotlib.pyplot as plt
import string
import numpy as np
import os, sys, getopt, glob, re, subprocess

def getFiles(runtag):
    """ Get list of all PSP files with given tag """
    tags = sorted(glob.glob("OUT.{}.*".format(runtag)))
    return tags

def getTime(file):
    p = subprocess.Popen(['pspinfo', file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    num = '([0-9.+\-e]+)'
    tim = re.compile('Time=' + num)
    res = tim.match(out)
    if res is not None:
        Time = float(res.group(1))
        return Time
    else:
        print "Error parsing: {}".format(file)
        return -1.0

def plotData(runtag, cname, col):
    t = []
    k = []
    e = []
    files = getFiles(runtag)
    for file in files:
        time  = getTime(file)
        os.system("psp2ascii {}".format(file))
        ifil = open("comp.{}".format(cname))
        line = ifil.readline()
        toks = line.split()
        numb = int(toks[0])     # Number of particles
        iatr = int(toks[1])     # Number of integer attributes
        datr = int(toks[2])     # Number of float attributes
        eion = 0.0              # Ion KE
        econ = 0.0              # Deferred E
        #
        # Parse the loop
        #
        for i in range(numb):
            line = ifil.readline()
            toks = line.split()
            # Compute the KE
            m = float(toks[0])
            v = [float(toks[4]), float(toks[5]), float(toks[6])]
            ke = 0.0
            for u in v: ke += 0.5*m*u*u
            eion += ke
            # Compute the deferred energy loss
            econ += float(toks[8+iatr+col])
        # Append data
        t.append(time)
        k.append(eion)
        e.append(econ)
        # Diagnostic printing
        print "name={}, T={}, KE={}, EC={}".format(runtag, time, eion, econ)

    print '-'*72
    print "Time"
    print '-'*72
    print t
    print '-'*72
    print "Ion KE"
    print '-'*72
    print k
    print '-'*72
    print "E cons"
    print '-'*72
    print e
    print '-'*72
    #
    # Make the plot
    #
    if len(t)>1:
        plt.plot(t, e, '-')
        plt.xlabel("Time")
        plt.ylabel("Energy")
        plt.show()

def main(argv):
    """ Parse the command line and call the parsing and plotting routine """

    info = '[-r <runtag> | --runtag=<degree> | -i <column> | --col <column>] <runtag>'

    runtag = "run"
    cname  = "gas"
    icol   = 7

    try:
        opts, args = getopt.getopt(argv,"r:c:i:h", ["runtag=", "cname=", "col=", "help"])
    except getopt.GetoptError:
        print sys.argv[0], info
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print sys.argv[0], info
            sys.exit()
        elif opt in ("-r", "--runtag"):
            runtag = arg
        elif opt in ("-c", "--cname"):
            cname = arg
        elif opt in ("-i", "--col"):
            icol = int(arg)

    plotData(runtag, cname, icol)

if __name__ == "__main__":
   main(sys.argv[1:])
