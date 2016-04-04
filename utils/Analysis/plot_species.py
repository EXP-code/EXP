#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import sys

argc = len(sys.argv)

#
# Last argument should be filename and must exist
#
if argc<=1:
        print "Usage: {} file".format(sys.argv[0])
        exit(1)

#
# Check for point type and log time
#
fmt  = '-'
logT = False
Temp = False
for i in range(1,argc-1):
        if sys.argv[i] == '-p' or sys.argv[i] == '--points':
                fmt = '-o'
        if sys.argv[i] == '-l' or sys.argv[i] == '--log':
                logT = True
        if sys.argv[i] == '-T' or sys.argv[i] == '--temp':
                Temp = True
#
# Parse data file
#
file = open(sys.argv[-1])
data = []
for line in file:
	if line.find('#') < 0:
		data.append([float(v) for v in line.split()])
a = np.array(data).transpose()
#
# Species plot
#
x = a[0]
if Temp: x = a[1]
if logT:
        plt.semilogx(x, a[2], fmt, label='H')
        plt.semilogx(x, a[3], fmt, label='H+')
        plt.semilogx(x, a[4], fmt, label='He')
        plt.semilogx(x, a[5], fmt, label='He+')
        plt.semilogx(x, a[6], fmt, label='He++')
else:
        plt.plot(x, a[2], fmt, label='H')
        plt.plot(x, a[3], fmt, label='H+')
        plt.plot(x, a[4], fmt, label='He')
        plt.plot(x, a[5], fmt, label='He+')
        plt.plot(x, a[6], fmt, label='He++')
plt.legend().draggable()
if Temp: plt.xlabel('Temperature')
else:    plt.xlabel('Time')
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
else:
        plt.semilogy(x, a[2], fmt, label='H')
        plt.semilogy(x, a[3], fmt, label='H+')
        plt.semilogy(x, a[4], fmt, label='He')
        plt.semilogy(x, a[5], fmt, label='He+')
        plt.semilogy(x, a[6], fmt, label='He++')
#
plt.legend().draggable()
if Temp: plt.xlabel('Temperature')
else:    plt.xlabel('Time')
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
else:    plt.xlabel('Time')
plt.ylabel('Species temperature')
plt.show()
