import numpy as np

# Read the output log file
data = np.loadtxt("OUTLOG.run0", skiprows=6, delimiter="|")

# Column 16 is -2T/VC.  The mean should be 1
mean = np.mean(data[:,16])
stdv = np.std (data[:,16])

# If the values are within 6 sigma of 1, assume that the simulation worked
if np.abs(mean - 1.0) > 6.0*stdv:
    exit(1)
else:
    exit(0)
