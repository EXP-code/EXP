#!/usr/bin/python

import re, sys, copy
import numpy as np
import matplotlib.pyplot as plt

#
# Patterns
#
num = '([0-9.+\-e]+)'
beg = '^[\[ ]+'
end = '[\] ]+'
sep = '[, ]+'
bar = '[\]| =]+'
pattern = beg + num + sep + num + end + bar + num

#
# Regex compilations
#
prog = re.compile(pattern)
begn = re.compile('.*Electron interaction energy.*')
ctim = re.compile('^[ ]+' + num + '[ ]+' + 'current time')
clev = re.compile('^[ ]+' + num + '[ ]+' + 'current level')
temp = re.compile('^[ ]+' + num + '[ ]+' + 'mass-weighted temperature')

file = open(sys.argv[1] + ".DSMC_log")

look = False

#
# Plotting data
#
xb = []
xe = []
xx = []
yy = []

#
# Default time
#
time = 0.0

slopeFac = 7736.337067418854
defaultT = float(sys.argv[2])
slope = slopeFac/defaultT
level = 0

until = 0.0
if len(sys.argv)==4:
    until = float(sys.argv[3])

for line in file:
    #
    # Look for current time
    #
    result = ctim.match(line)
    if result is not None:
        time = float(result.group(1))
    #
    # Look for current temp
    #
    result = temp.match(line)
    if result is not None:
        ttemp = float(result.group(1))
        if ttemp>0.0:
            defaultT = ttemp
            slope = slopeFac/defaultT
    #
    # Look for current level
    #
    result = clev.match(line)
    if result is not None:
        level = int(result.group(1))
    #
    # If in a desired stanza, look for data
    #
    if look:
        result = prog.match(line)
        if result is not None:
            if len(result.groups()) == 3:
                xb.append(float(result.group(1)))
                xe.append(float(result.group(2)))
                xx.append(0.5*(xb[-1]+xe[-1]))
                yy.append(float(result.group(3)) + 0.1)
        elif len(xb)>0:         # End of data: make the plot
            ly = []
            for v in yy: ly.append(np.log(v))
            lz = []
            mid = len(ly)/2
            for v in xx: lz.append(ly[mid]-slope*(v-xx[mid]))

            zz = []
            for v in lz: zz.append(np.exp(v))

            fig, axes = plt.subplots(nrows=2, ncols=1)

            ax = axes[0]
            ax.plot(xx, ly, '-o')
            ax.plot(xx, lz, '-')

            lz = []
            for v in xx: lz.append(ly[mid]-slope*1.2*(v-xx[mid]))

            ax.plot(xx, lz, '-')

            lz = []
            for v in xx: lz.append(ly[mid]-slope*0.8*(v-xx[mid]))

            ax.plot(xx, lz, '-')
            ax.set_title("interaction: t={} T={}[+/-]{}".format(time,defaultT,0.2*defaultT))
            ax.set_ylabel("Log(counts)")
            ax.tick_params(axis='x', labelbottom='off')

            ay = axes[1]
            ay.plot(xx, yy, '-o')
            ay.plot(xx, zz, '-')
            ay.set_xlabel("Energy (eV)")
            ay.set_ylabel("Counts")

            fig.tight_layout()
            plt.show()
            look = False
    else:                       # Look for beginning of stanza
        result = begn.match(line)
        if result is not None and time >= until and level==0:
            xb = []             # Zero out data
            xe = []
            xx = []
            yy = []
            look = True
