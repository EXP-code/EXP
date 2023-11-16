import numpy as np

# Read the output log file
data = np.loadtxt("OUTLOG.run0", skiprows=6, delimiter="|")

# Column 16 is -2T/VC.  The mean should be 1 with a std dev < 0.03
mean = np.mean(data[:,16])
stdv = np.std (data[:,16])

# If the values are within 3 sigma, assume that the simulation worked
if np.abs(mean - 1.0) > 3.0*stdv:
    exit(1)
else:
    exit(0)
