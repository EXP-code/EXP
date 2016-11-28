import sys
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
        
fields = [ ['Ions_E', 'Eion(1)', 'Eion(2)'], ['Elec_E', 'Eelc(1)', 'Eelc(2)'], 'Totl_E', 'Cons_G', 'Cons_E']

d = {}
for v in labs:
	d[v] = gs.readDB(v)

totl = len(fields)
cols = 2
rows = totl/cols
if cols * rows < totl: rows += 1
f, ax = plt.subplots(rows, cols, sharex='col')

for x in ax:
        for y in x:
                box = y.get_position()
                newBox = [box.x0, box.y0, 0.8*box.width, box.height]
                y.set_position(newBox)

ax[rows-1, cols-1].axis('off')

def plotme(d, v, f, ax, l):
        indx = np.searchsorted(d[v]['Time'], tmax)
        ax.plot(d[v]['Time'][0:indx], d[v][f][0:indx], '-', label=l)
        if f == 'Totl_E':
                tt = d[v]['Ions_E'][0:indx] + d[v]['Elec_E'][0:indx]
                ax.plot(d[v]['Time'][0:indx], tt, '-', label=v+':comb')

cnt = 0
for f in fields:
        row = cnt/cols
        col = cnt - row*cols
        for v in labs:
                if type(f) is list:
                        for ff in f:
                                plotme(d, v, ff, ax[row, col], v+':'+ff)
                else:
                        plotme(d, v, f, ax[row, col], v)
                                
        if type(f) is list:
                ax[row, col].set_title(f[0])
        else:
                ax[row, col].set_title(f)
        ax[row, col].legend(prop={'size':8}, bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.0)
        if col==0:
                ax[row, col].set_ylabel('Energy')
                if row+1==rows: ax[row, col].set_xlabel('Time')
        cnt += 1

plt.show()
