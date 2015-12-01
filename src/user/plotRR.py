import sys
import matplotlib.pyplot as plt
import numpy as np

for i in range(1,len(sys.argv)):
        name = sys.argv[i]
        file = open(name)
        data = []
        for line in file:
                data.append([float(v) for v in line.split()])
        a = np.array(data).transpose()
        lx = np.log(a[0])/np.log(10.0)
        ly = np.log(a[2])/np.log(10.0)
        plt.plot(lx, ly, '-',label=name)

plt.xlabel("E(eV)")
plt.ylabel("Mb")
plt.legend()
plt.show()
