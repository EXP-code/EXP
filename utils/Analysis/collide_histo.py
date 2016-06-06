#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import sys
import psp_io

#
# Last argument should be filename and must exist
#

help_string = "Usage: {} [-s n | --species n] [-c n | -coll n]  [-k n | -key n] [-d | -diff] [-o n | --offset n] runtag".format(sys.argv[0])

argc = len(sys.argv)

if argc<=1:
        print help_string
        exit(1)

#
# Check for point type and log time
#
spc  = 2
key  = 0
col  = 1
diff = False
offs = 1

for i in range(1,argc):
        if sys.argv[i] == '-k' or sys.argv[i] == '--key':
                key = int(sys.argv[i+1])
        if sys.argv[i] == '-s' or sys.argv[i] == '--species':
                spc = int(sys.argv[i+1])
        if sys.argv[i] == '-c' or sys.argv[i] == '--coll':
                col = int(sys.argv[i+1])
        if sys.argv[i] == '-d' or sys.argv[i] == '--diff':
                diff = True
        if sys.argv[i] == '-o' or sys.argv[i] == '--offset':
                offs = int(sys.argv[i+1])
        if sys.argv[i] == '-h' or sys.argv[i] == '--help':
                print help_string
                exit(1)
#
# Parse data file
#
if diff:
        hi   = int(sys.argv[-1])
        lo   = max(hi - offs, 0)
        psp0 = psp_io.Input('OUT.{}.{:05d}'.format(sys.argv[-2], lo), comp='gas')
        psp1 = psp_io.Input('OUT.{}.{:05d}'.format(sys.argv[-2], hi), comp='gas')
        s   = getattr(psp0, 'i{}'.format(key))
        c0  = getattr(psp0, 'i{}'.format(col))
        c1  = getattr(psp1, 'i{}'.format(col))
else:
        nd   = int(sys.argv[-1])
        psp1 = psp_io.Input('OUT.{}.{:05d}'.format(sys.argv[-2], nd), comp='gas')
        s    = getattr(psp1, 'i{}'.format(key))
        c1   = getattr(psp1, 'i{}'.format(col))

data = []
for i in range(len(s)):
    if s[i] == spc:
            if diff:
                    data.append(float(c1[i] - c0[i]))
            else:
                    data.append(float(c1[i]))

a = np.array(data);

minC = min(a)
maxC = max(a)           

num  = int(np.sqrt(len(data)))
dC   = (maxC - minC)/(num - 1)

hist = np.zeros(num)

for v in a:
    indx = int( (v - minC)/dC )
    hist[indx] += 1

for v in range(num):
    print "{:13.4g} {:13.4g}  :  {}".format(minC + dC*v, minC + dC*(v+1), hist[v])

n, bins, patches = plt.hist(a, num, normed=1, alpha=0.75)

plt.xlabel("Collisions")
plt.ylabel("Frequency")
plt.grid(True)
plt.show()

    

    
    
