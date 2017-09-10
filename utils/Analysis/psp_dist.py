#!/usr/bin/python

# -*- coding: utf-8 -*-

"""Program to compute the energy distributions based on the DSMC_log
file data

Examples:

	$ python psp_dist.py -f electron run2

Plots the energy distributions for the named field (in the case above,
electrons) for the run with tag "run2".  Field value is "electron" by
default.  Other fields are "ion" and "interact" for the electron-ion
interaction kinetic energy.

Only for Trace method, so far.
"""

import os, re, sys, copy, getopt
import flufl.enum as enum
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import psp_io

class FitType(enum.Enum):
    Analytic = 1
    AmpOnly  = 2
    TempAmp  = 3
    Fixed    = 4

class PlotType(enum.Enum):
    Linear = 1
    Log    = 2
    Both   = 3

Munit = 1.9891e33
Lunit = 3.08568e18
Tunit = 3.15569e10
Vunit = Lunit/Tunit

m_H   = 1.008
m_He  = 4.002602
X_H   = 0.76
X_He  = 0.24
Mu_i  = 1.0/(X_H/m_H + X_He/m_He)
Mu_e  = 0.000548579909

amu   = 1.660539e-24
eV    = 1.60217653e-12

b0    = 0.0
a0    = 1.0

def func1(x, a):
    """Fit energy distribution for amplitude only"""
    global b0
    if a<=0.0: return 1e30
    return a * np.sqrt(x) * np.exp(-b0 * x) / b0**1.5

def func2(x, a, b):
    """Fit energy distribution for amplitude and temperature"""
    if a<=0.0: return 1e30
    if b<=0.0: return 1e30
    return a * np.sqrt(x) * np.exp(-b * x) / b**1.5

def plot_data(runtag, field, defaultT, start, stop, ebeg, efin, dE, fit, pTyp):
    """Parse and plot the OUT psp output files
 
    Parameters:

    runtag (string): is the input datafile name

    field (string): is either "electron", "ion" or "interact"

    defaultT (float): is the default temperature (probably not needed
    in most cases)

    start (int): initial psp index

    stop(int): final pspindex

    ebeg(float): initial energy

    efin(float): final energy

    delta(float): energy interval

    fit(FitType): how to chose comparison curve

    pTyp(PlotType): type of plots

    """
    global b0

    slopeFac = 11604.50560112828
    slope = slopeFac/defaultT

    #
    # Set up bins
    #
    nbin = int( (efin - ebeg)/dE )
    efin = ebeg + dE*nbin
    yy = np.zeros(nbin)
    oab = 0

    #
    # Loop through phase space
    #
    
    efac = 0.5*Vunit*Vunit*amu/eV
    if field=='ion':
        efac *= Mu_i
    elif field=='interact':
        efac *= Mu_i*Mu_e/(Mu_i + Mu_e)
    else:
        efac *= Mu_e
        
    for n in range(start, stop+1):
        filename = 'OUT.{}.{:05d}'.format(runtag,n)
        if not os.path.isfile(filename): continue

        O = psp_io.Input(filename, comp='gas')

        dv = np.zeros(3)
        for i in range(O.mass.size):
            if field=='ion':
                dv[0] = O.xvel[i]
                dv[1] = O.yvel[i]
                dv[2] = O.zvel[i]
            elif field=='interact':
                dv[0] = O.xvel[i] - O.d12[i]
                dv[1] = O.yvel[i] - O.d13[i]
                dv[2] = O.zvel[i] - O.d14[i]
            else:
                dv[0] = O.d12[i]
                dv[1] = O.d13[i]
                dv[2] = O.d14[i]

            EE = efac * np.dot(dv, dv)
            if EE<=ebeg or EE>=efin:
                oab += 1
            else:
                indx = int( (EE - ebeg)/dE)
                yy[indx] += 1

    # Compute bin centers
    xx = np.zeros(nbin)
    for i in range(nbin): xx[i] = ebeg + dE*(0.5+i)

    # Normalize
    norm = (sum(yy) + oab)*dE
    for i in range(nbin): yy[i] /= norm

    # Fit for temperature
    fc = 11604.5/defaultT # Initial guess for exponent
    tt = []
    if fit == FitType.Fixed:
        # Temp
        bT = defaultT
        # Inverse temp
        b0 = fc
        # 
        tt = np.zeros(nbin)
        for i in range(nbin): tt[i] = func1(xx[i], a0)
    elif fit == FitType.Analytic:
        # Temp
        bT = defaultT
        # Inverse temp
        b0 = fc
        # 
        tt = np.zeros(nbin)
        nrm = 2.0*b0**3/np.sqrt(np.pi)
        for i in range(nbin): tt[i] = func1(xx[i], nrm)
    elif fit == FitType.TempAmp:
        p0 = [sum(yy)/len(yy)*fc**1.5,fc] # Amplitude
        popt, pcov = curve_fit(func2, xx, yy, p0, sigma=np.sqrt(yy)+0.1)
        # Temp
        bT = slopeFac / popt[1]
        for v in xx: tt.append(func2(v, popt[0], popt[1]))
        print 'Amplitude={}  Temperature={}'.format(popt[0], bT)
    elif fit == FitType.AmpOnly:
        b0 = fc
        p0 = [sum(yy)/len(yy)*fc**1.5] # Amplitude
        popt, pcov = curve_fit(func1, xx, yy, p0, sigma=np.sqrt(yy)+1)
        # Temp
        bT = defaultT
        for v in xx: tt.append(func1(v, popt[0]))
        print 'Amplitude={}  Temperature={}'.format(popt[0], bT)
    else:
        print "This is impossible"
        exit

    # Make curves
    lt = []
    for v in tt: lt.append(np.log(v+1))
    #
    ly = []
    for v in yy: ly.append(np.log(v+1))
    
    if pTyp == PlotType.Both:
        fig, axes = plt.subplots(nrows=2, ncols=1)
    else:
        fig, axes = plt.subplots(nrows=1, ncols=1)
    
    if pTyp.value & PlotType.Log.value:

        if pTyp == PlotType.Both:  ax = axes[0]
        else:                      ax = axes

        ax.plot(xx, ly, '-o')
        ax.plot(xx, lt, '-')
    
        if type==FitType.TempAmp:
            ax.set_title("{}: T(fit)={}".format(field,bT))
        else:
            ax.set_title("{}: T={}".format(field,bT))

        if pTyp == PlotType.Log:
            ax.set_xlabel("Energy (eV)")
        else:
            ax.tick_params(axis='x', labelbottom='off')
        ax.set_ylabel("Log(counts)")

    if pTyp.value & PlotType.Linear.value:

        if pTyp == PlotType.Both: ay = axes[1]
        else:                     ay = axes

        ay.plot(xx, yy, '-o')
        ay.plot(xx, tt, '-')
        
        ay.set_xlabel("Energy (eV)")
        ay.set_ylabel("Counts")
        
    fig.tight_layout()
    plt.show()


def main(argv):
    """ Parse the command line and call the parsing and plotting routine """

    global a0

    field = "electron"
    start = 0
    stop  = 99999
    ebeg  = 0.05
    efin  = 60.0
    delta = 0.05
    defT  = 100000.0
    fit   = FitType.Analytic
    plt   = PlotType.Both

    options = '[-f <field> | --field=<field> | -n <start> | --start=<start> | -N <stop> | --stop=<stop> | -e <low eV> | --low=<low eV> | -E <high eV> | --high=<high eV> | -d <delta eV> | --delta=<delta eV> | -T <default T> | --temp=<default T> | -t | --type | -F | --fixed | -p <type> | --plot=<type>] <runtag>'

    try:
        opts, args = getopt.getopt(argv,"hf:n:N:e:E:d:T:t:Fp:", ["help","field=","start=","stop=","low=", "high=", "delta=", "temp=", "type=", "fixed=", "plot="])
    except getopt.GetoptError:
        print sys.argv[0], 
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print sys.argv[0], options
            sys.exit()
        elif opt in ("-f", "--field"):
            field = arg
        elif opt in ("-n", "--start"):
            start = int(arg)
        elif opt in ("-N", "--stop"):
            stop = int(arg)
        elif opt in ("-e", "--low"):
            ebeg = float(arg)
        elif opt in ("-E", "--high"):
            efin = float(arg)
        elif opt in ("-d", "--delta"):
            delta = float(arg)
        elif opt in ("-T", "--temp"):
            defT = float(arg)
        elif opt in ("-t", "--type"):
            if   FitType.Analytic.name == arg: fit = FitType.Analytic
            elif FitType.AmpOnly.name  == arg: fit = FitType.AmpOnly
            elif FitType.TempAmp.name  == arg: fit = FitType.TempAmp
            else:
                print "No such fit type: ", arg
                print "Valid types are:"
                for v in FitType: print v.name
                exit
        elif opt in ("-F", "--fixed"):
            a0 = float(arg)
            fit = FitType.Fixed
        elif opt in ("-p", "--plot"):
            if   PlotType.Linear.name == arg: plt = PlotType.Linear
            elif PlotType.Log.name    == arg: plt = PlotType.Log
            elif PlotType.Both.name   == arg: plt = PlotType.Type
            else:
                print "No such plot type: ", arg
                print "Valid types are:"
                for v in PlotType: print v.name
                exit

    if len(args)>0:
        runtag = args[0]
    else:
        runtag = "run"

    plot_data(runtag, field, defT, start, stop, ebeg, efin, delta, fit, plt)

if __name__ == "__main__":
   main(sys.argv[1:])
                