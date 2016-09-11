#!/usr/bin/python

# -*- Python -*-
# -*- coding: utf-8 -*-

"""Program to display the ionization state for atomic type
using diagnostic output from the CollideIon class in the
UserTreeDSMC module.

Example:

	$ plot_species.py run2

"""

import sys, getopt
import numpy as np
import matplotlib.pyplot as plt
import sys
import plotColl as pC

def scanSeq(a):
        """ 
        Attempt to put sequence in order.  
        Returns False if no changes are made
        """
        x = a[0,0]
        for i in range(1, a.shape[1]):
                if a[0,i] < x:
                        # Find location
                        k = 0
                        for j in range(i):
                                if a[0,j] > a[0,i]:
                                        k = j
                                        break
                        # print "a={}, b={}, c={}".format(x, a[0,i], a[0,k])
                        return True, (k, i)
                x = a[0,i]
        return False, (0, a.shape[1])

def scanLoc(t, a):
        """ 
        Find index k, such that a[0,k]>=t
        Returns False if it does not exist.
        """
        k = 0
        for i in range(1, a.shape[1]):
                if t<=a[0,i] and t>a[0,i-1]:
                        return True, i
                x = a[0,i]
        return False, a.shape[1]

def plot_data(argv):
        """
        Parse and plot the *.species file
        """

        #
        # Check for point type and log time
        #
        fmt  = '-'
        logT = False
        Temp = False
        Tscl = 1
        Tbeg = 0.0
        Tend = -1.0

        argc = len(argv)

        for i in range(1,argc):
                if argv[i] == '-p' or argv[i] == '--points':
                        fmt = '-o'
                if argv[i] == '-l' or argv[i] == '--log':
                        logT = True
                if argv[i] == '-T' or argv[i] == '--temp':
                        Temp = True
                if argv[i] == '-t' or argv[i] == '--timescale':
                        Tscl = float(argv[i+1])
                if argv[i] == '-b' or argv[i] == '--Tbeg':
                        Tbeg = float(argv[i+1])
                if argv[i] == '-e' or argv[i] == '--Tend':
                        Tend = float(argv[i+1])
                if argv[i] == '-h' or argv[i] == '--help':
                        print "Usage: {} [-p|--points] [-l|--log] [-t scale|--timescale scale] [-b|--Tbeg] [-e|--Tend] [-h|--help] runtag".format(argv[0])
                        return
        #
        # Parse data file
        #
        try:
                file = open(argv[-1] + ".species")
        except:
                print "Error opening file <{}>".format(argv[-1] + ".species")
                return
        data = []
        nvec = 0
        for line in file:
                if line.find('#') < 0:
                        pvec = [float(v) for v in line.split()]
                        lvec = len(pvec)
                        # Enforce equal size vectors; data line integrity
                        #
            if nvec==0: nvec = lvec
                        if nvec == lvec: data.append(pvec)
        a = np.array(data).transpose()
        #
        # Search for appended sequence(s)
        #
        print "-------- Trim statistics --------"
        print "Initial shape:", a.shape

        status = True
        segcnt = 0
        while status:
                status, segment = scanSeq(a)
                if status:
                        b = a[:,np.s_[0:segment[0]:]]
                        c = a[:,np.s_[segment[1]::]]
                        a = np.concatenate((b, c), axis=1)
                        segcnt += 1
                        print "Segment [{:3d}]: ({}, {})".format(segcnt,
                                                                 segment[0],
                                                                 segment[1])
        if Tbeg > 0.0:
                status, indx = scanLoc(Tbeg/Tscl, a)
        if status: a = a[:,np.s_[indx::]]

        if Tend > Tbeg:
                status, indx = scanLoc(Tend/Tscl, a)
                if status: a = a[:,np.s_[:indx:]]


        print "  Final shape:", a.shape
        print "---------------------------------"

        pC.readDB(argv[-1])
        pc_time = pC.db['Time']
        pc_frac = pC.db['Efrac']
        pc_len  = len(pc_time)
        if pc_len != len(pc_frac) or pc_len==0: pc_len = 0

        #
        # Species plot
        #
        x = a[0] * Tscl
        if Temp: x = a[1]
        if logT:
                plt.semilogx(x, a[2], fmt, label='H')
                plt.semilogx(x, a[3], fmt, label='H+')
                plt.semilogx(x, a[4], fmt, label='He')
                plt.semilogx(x, a[5], fmt, label='He+')
                plt.semilogx(x, a[6], fmt, label='He++')
                if pc_len:
                        plt.semilogx(pc_time*Tscl, pc_frac, fmt, label='e')
        else:
                plt.plot(x, a[2], fmt, label='H')
                plt.plot(x, a[3], fmt, label='H+')
                plt.plot(x, a[4], fmt, label='He')
                plt.plot(x, a[5], fmt, label='He+')
                plt.plot(x, a[6], fmt, label='He++')
                if pc_len:
                        plt.plot(pc_time*Tscl, pc_frac, fmt, label='e')
        plt.legend().draggable()
        if Temp: plt.xlabel('Temperature')
        else:
                if Tscl == 1:
                        plt.xlabel('Time')
                else:
                        plt.xlabel('Time (years)')
        plt.ylabel('Species')
        plt.show()
        #
        # Species plot
        #
        if logT:
                plt.loglog(x, a[2], fmt, label='H')
                plt.loglog(x, a[3], fmt, label='H+')
                plt.loglog(x, a[4], fmt, label='He')
                plt.loglog(x, a[5], fmt, label='He+')
                plt.loglog(x, a[6], fmt, label='He++')
                if pc_len:
                        plt.loglog(pc_time*Tscl, pc_frac, fmt, label='e')
        else:
                plt.semilogy(x, a[2], fmt, label='H')
                plt.semilogy(x, a[3], fmt, label='H+')
                plt.semilogy(x, a[4], fmt, label='He')
                plt.semilogy(x, a[5], fmt, label='He+')
                plt.semilogy(x, a[6], fmt, label='He++')
                if pc_len:
                        plt.semilogy(pc_time*Tscl, pc_frac, fmt, label='e')
        #
        plt.legend().draggable()
        if Temp: plt.xlabel('Temperature')
        else:
                if Tscl == 1:
                        plt.xlabel('Time')
                else:
                        plt.xlabel('Time (years)')

        plt.ylabel('Species')
        plt.show()
        #
        # Temperature plot
        #
        if logT:
                plt.semilogx(x, a[1 ], fmt, label='Temp')
                plt.semilogx(x, a[10], fmt, label='Temp_e')
                plt.semilogx(x, a[16], fmt, label='Temp(1)_i')
                plt.semilogx(x, a[22], fmt, label='Temp(1)_e')
                plt.semilogx(x, a[28], fmt, label='Temp(2)_i')
                plt.semilogx(x, a[34], fmt, label='Temp(2)_e')
        else:
                plt.plot(x, a[1 ], fmt, label='Temp')
                plt.plot(x, a[10], fmt, label='Temp_e')
                plt.plot(x, a[16], fmt, label='Temp(1)_i')
                plt.plot(x, a[22], fmt, label='Temp(1)_e')
                plt.plot(x, a[28], fmt, label='Temp(2)_i')
                plt.plot(x, a[34], fmt, label='Temp(2)_e')
        #
        plt.legend().draggable()
        if Temp: plt.xlabel('Temperature')
        else:
                if Tscl == 1:
                        plt.xlabel('Time')
                else:
                        plt.xlabel('Time (years)')

        plt.ylabel('Species temperature')
        plt.show()


def main(argv):
        """ Parse the command line and call the parsing and plotting routine """

        #
        # Last argument should be filename and must exist
        #
        if len(argv) <=1:
                print "Usage: {} runtag".format(argv[0])
                exit(1)

        plot_data(argv)

if __name__ == "__main__":
        main(sys.argv)
