import numpy as np

# Read the output log file
data = np.loadtxt("OUTLOG.run0", skiprows=6, delimiter="|")

# Column 16 is -2T/VC.  The mean should be 1 with some small number of std dev
mean = np.mean(data[:,16])
stdv = np.std (data[:,16])

# print("Mean = " + str(mean) + ", std dev = " + str(stdv))
# print("Check = {}".format(np.abs(mean-1.0) - 5.0*stdv))

# If the values are within 3 sigma, assume that the simulation worked
if np.abs(mean - 1.0) > 5.0*stdv:
    exit(1)
else:
    exit(0)
