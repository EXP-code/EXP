import sys, os
import numpy as np
import matplotlib.pyplot as plt
import getSpecies as gs

def is_number(s):
        try:
                float(s)
                return True
        except ValueError:
                return False

tmax   = 1000000.0
if len(sys.argv)==1:
        print "Usage: ", sys.argv[0], " tag1 tag2 . . . [Tmax]"
        os.exit(-1)

labs = []
if is_number(sys.argv[-1]):
        tmax = float(sys.argv[-1])
        labs = sys.argv[1:-1]
else:
        labs = sys.argv[1:]
        
fields = [ ['Telc(1)', 'Telc(2)'], ['Tion(1)', 'Tion(2)'] ]
titles = ['Electrons', 'Ions']

d = {}
for v in labs:
	d[v] = gs.readDB(v)

f, ax = plt.subplots(2, 1, sharex='col')
for x in ax:
        box = x.get_position()
        newBox = [box.x0, box.y0, 0.9*box.width, box.height]
        x.set_position(newBox)

cnt = 0
for v in labs:
        for n in range(2):
                for F in fields[n]:
                        if F in d[v]:
                                indx = np.searchsorted(d[v]['Time'], tmax)
                                ax[n].plot(d[v]['Time'][0:indx], d[v][F][0:indx], '-', label=v+':'+F)

                ax[n].legend(prop={'size':8}, bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.0)
                if n>0: ax[n].set_xlabel('Time')
                ax[n].set_ylabel('Temperature')
                ax[n].set_title(titles[n])
plt.show()
