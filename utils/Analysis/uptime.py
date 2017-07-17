#!/usr/bin/python

import os, sys, re
import subprocess

pattern = re.compile("[0-9-]+")

def parse(seq):
    toks = seq.split(',')
    n = []
    for t in toks:
        if pattern.match(t):
            be = t.split('-')
            if len(be)==1:
                n.append(int(be[0]))
            if len(be)==2:
                for i in range(int(be[0]), int(be[1])+1):
                    n.append(i)
    return n

sep = '---------+--------'

print '{:7s}  |  {:5s}'.format('Host', '1 min')
print sep

n = parse(sys.argv[1])
for i in n:
    p = subprocess.Popen(['ssh', 'node{:03d}'.format(i), 'uptime'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    start = out.find('load average:')
    if start >=0:
        want = out[start+len('load average:'):].split(',')
        # print 'node{:03d}:  {:5.2f}  {:5.2f}  {:5.2f}'.format(i, float(want[0]), float(want[1]), float(want[2]))
        print 'node{:03d}  |  {:5.2f}'.format(i, float(want[0]))

print sep
