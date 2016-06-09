#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import os.path
import psp_io
import sys
from scipy.optimize import curve_fit

def func(x, a, b):
    """Fit energy distribution for amplitude and temperature"""
    if a<=0.0: return 1e30
    if b<=0.0: return 1e30
    return a * np.sqrt(x) * np.exp(-b * x) / b**1.5

#
# Last argument should be filename and must exist
#

help_string = "Usage: {} [-k n | --key n] [--beg n] [--end n] [--stride n] [-e n | --elec n] [-p | --plot] [-l | --log] [-h | --help] runtag".format(sys.argv[0])

argc = len(sys.argv)

if argc<=1:
        print help_string
        exit(1)

#
# Check for point type and log time
#
key  = 0
epos = 10
nbeg = 0
nend = -1
strd = 1
tstP = False
logP = False

for i in range(1,argc):
    if sys.argv[i] == '-k' or sys.argv[i] == '--key':
        key = int(sys.argv[i+1])
    if sys.argv[i] == '--beg':
        nbeg = int(sys.argv[i+1])
    if sys.argv[i] == '--end':
        nend = int(sys.argv[i+1])
    if sys.argv[i] == '--stride':
        strd = int(sys.argv[i+1])
    if sys.argv[i] == '-e' or sys.argv[i] == '--elec':
        epos = int(sys.argv[i+1])
    if sys.argv[i] == '-p' or sys.argv[i] == '--plot':
        tstP = True
    if sys.argv[i] == '-l' or sys.argv[i] == '--log':
        logP = True
    if sys.argv[i] == '-h' or sys.argv[i] == '--help':
        print help_string
        exit(1)
#
# Parse data files
#

filep = 'OUT.{}.{:05d}'
count = nbeg

atomic_masses = [0.000548579909, 1.00794, 4.002602];

# atomic mass unit
amu = 1.660011e-24

# Parsec (in cm)
pc  = 3.08567758e18

# Solar mass (in g)
msun = 1.9891e33

# Seconds per year
year = 365.242*24.0*3600.0

# electron volt in (cgs)
eV =  1.60217653e-12

Munit = msun*0.1
Lunit = pc
Tunit = year*1e5
Vunit = Lunit/Tunit

e2ev = 0.5*amu*Vunit*Vunit/eV

slopeFac = 11604.50560112828

ionsT = {}
elecT = {}

while os.path.isfile(filep.format(sys.argv[-1], count)):

    psp = psp_io.Input(filep.format(sys.argv[-1], count), comp='gas')
    s   = getattr(psp, 'i{}'.format(key))
    ex  = [getattr(psp, 'd{}'.format(epos+0)), getattr(psp, 'd{}'.format(epos+1)), getattr(psp, 'd{}'.format(epos+2))]
    Ei  = []
    Ee  = []
    Nf  = []
    Ms  = {}
    Mm  = {}
    for i in range(len(s)):
        Ei.append(e2ev*atomic_masses[s[i]]*(psp.xvel[i]**2 + psp.yvel[i]**2 + psp.zvel[i]**2))
        Ee.append(e2ev*atomic_masses[0]*(ex[0][i]**2 + ex[1][i]**2  + ex[2][i]**2))
        Nf.append(psp.mass[i]/atomic_masses[s[i]])
        if s[i] not in Ms:
            Ms[s[i]] = 0.0
            Mm[s[i]] = Nf[-1]
        Ms[s[i]] += Nf[-1]


    time = psp.ctime
    for ss in Ms:
        if ss not in ionsT: ionsT[ss] = [[],[]]
        if ss not in elecT: elecT[ss] = [[],[]]

    # Bin distribution
    for ss in Ms:
        spE = []
        spI = []
        for i in range(len(s)):
            if ss==s[i]:
                spI.append(Ei[i])
                spE.append(Ee[i])
        
        minE = min(spI)
        maxE = max(spI)
        delE = (maxE - minE)/100.0
        xx = np.arange(minE, maxE, delE)
        yy = np.zeros(len(xx))
        for v in spI:
            indx = int( (v - minE)/delE )
            if indx >= 0 and indx < len(yy): yy[indx] += 1
            
        popt, pcov = curve_fit(func, xx, yy)

        ionsT[ss][0].append(time)
        ionsT[ss][1].append(slopeFac/popt[1])

        if tstP:
            tt = []
            for j in range(len(xx)): tt.append(func(xx[j], popt[0], popt[1]))
            if logP:
                plt.semilogy(xx, yy, '-', label="data")
                plt.semilogy(xx, tt, '-', label="fit")
            else:
                plt.plot(xx, yy, '-', label="data")
                plt.plot(xx, tt, '-', label="fit")
            plt.xlabel("Energy (eV)")
            plt.ylabel("Counts")
            plt.title("Ion (Z={}): t={} T={}".format(ss, time, slopeFac/popt[1]))
            plt.legend()
            plt.show()
            
        minE = min(spE)
        maxE = max(spE)
        delE = (maxE - minE)/100.0
        xx = np.arange(minE, maxE, delE)
        yy = np.zeros(len(xx))
        for v in spE:
            indx = int( (v - minE)/delE )
            if indx >= 0 and indx < len(yy): yy[indx] += 1
            
        popt, pcov = curve_fit(func, xx, yy)

        elecT[ss][0].append(time)
        elecT[ss][1].append(slopeFac/popt[1])
        
        if tstP:
            tt = []
            for j in range(len(xx)): tt.append(func(xx[j], popt[0], popt[1]))
            if logP:
                plt.semilogy(xx, yy, '-', label="data")
                plt.semilogy(xx, tt, '-', label="fit")
            else:
                plt.plot(xx, yy, '-', label="data")
                plt.plot(xx, tt, '-', label="fit")
            plt.xlabel("Energy (eV)")
            plt.ylabel("Counts")
            plt.title("Electron (Z={}): t={} T={}".format(ss, time, slopeFac/popt[1]))
            plt.legend()
            plt.show()


    count += strd
    if nend>=0 and count>=nend: break

for s in ionsT:
    plt.plot(ionsT[s][0], ionsT[s][1], '-', label="Ion (Z={})".format(s))
    plt.plot(elecT[s][0], elecT[s][1], '-', label="Electron (Z={})".format(s))

plt.xlabel("Time")
plt.ylabel("Temperature")
plt.legend().draggable()
plt.show()
