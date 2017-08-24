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

import re, sys, copy, getopt
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import psp_io


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

def plot_data(runtag, field, defaultT, start, stop, ebeg, efin, dE, fixed):
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

    """
    global b0

    slopeFac = 11604.50560112828
    slope = slopeFac/defaultT

    #
    # Set up bins
    #
    nbin = int( (efin - ebeg)/dE )
    efin = efin + dE*nbin
    yy = np.zeros(nbin)

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
        O = psp_io.Input(filename, comp='gas')

        for i in range(O.mass.size):
            EE = 0.0
            if field=='ion':
                EE = O.xvel[i]*O.xvel[i] + O.yvel[i]*O.yvel[i] + O.zvel[i]*O.zvel[i]
            elif field=='interact':
                EE = (O.xvel[i] - O.d12)*(O.xvel[i] - O.d12) + (O.yvel[i] - O.d13)*(O.yvel[i] - O.d13) + (O.zvel[i] - O.d14)*(O.zvel[i] - O.d14)
            else:
                EE = O.d12[i]*O.d12[i] + O.d13[i]*O.d13[i] + O.d14[i]*O.d14[i]

            EE *= efac
            indx = int( (EE - ebeg)/dE)
            if indx>=0 and indx<nbin: yy[indx] += 1

    # Compute bin centers
    xx = np.zeros(nbin)
    for i in range(nbin): xx[i] = ebeg + dE*(0.5+i)

    # Fit for temperature
    fc = 11604.5/defaultT # Initial guess for exponent
    tt = []
    if not fixed:
        p0 = [sum(yy)/len(yy)*fc**1.5,fc] # Amplitude
        popt, pcov = curve_fit(func2, xx, yy, p0, sigma=np.sqrt(yy)+0.1)
        # Temp
        bT = slopeFac / popt[1]
        for v in xx: tt.append(func2(v, popt[0], popt[1]))
    else:
        b0 = fc
        p0 = [sum(yy)/len(yy)*fc**1.5] # Amplitude
        popt, pcov = curve_fit(func1, xx, yy, p0, sigma=np.sqrt(yy)+1)
        # Temp
        bT = defaultT
        for v in xx: tt.append(func1(v, popt[0]))

    # Make curves
    lt = []
    for v in tt: lt.append(np.log(v+1))
    #
    ly = []
    for v in yy: ly.append(np.log(v+1))
    
    fig, axes = plt.subplots(nrows=2, ncols=1)
    
    ax = axes[0]
    ax.plot(xx, ly, '-o')
    ax.plot(xx, lt, '-')
    
    if fixed:
        ax.set_title("{}: T={}".format(field,bT))
    else:
        ax.set_title("{}: T(fit)={}".format(field,bT))
    ax.set_ylabel("Log(counts)")
    ax.tick_params(axis='x', labelbottom='off')

    ay = axes[1]
    ay.plot(xx, yy, '-o')
    ay.plot(xx, tt, '-')
        
    ay.set_xlabel("Energy (eV)")
    ay.set_ylabel("Counts")
        
    fig.tight_layout()
    plt.show()


def main(argv):
    """ Parse the command line and call the parsing and plotting routine """

    field = "electron"
    start = 0
    stop  = 1000000
    ebeg  = 0.05
    efin  = 60.0
    delta = 0.05
    defT  = 100000.0
    fixed = False

    options = '[-f <field> | --field=<field> | -n <start> | --start=<start> | -N <stop> | --stop=<stop> | -e <low eV> | --low=<low eV> | -E <high eV> | --high=<high eV> | -d <delta eV> | --delta=<delta eV> | -T <default T> | --temp=<default T> | -F | --fixed] <runtag>'

    try:
        opts, args = getopt.getopt(argv,"hf:n:N:e:E:d:T:F", ["help","field=","start=","stop=","low=", "high=", "delta=", "temp=", "fixed"])
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
        elif opt in ("-F", "--fixed"):
            fixed = True

    if len(args)>0:
        runtag = args[0]
    else:
        runtag = "run"

    plot_data(runtag, field, defT, start, stop, ebeg, efin, delta, fixed)

if __name__ == "__main__":
   main(sys.argv[1:])
                
