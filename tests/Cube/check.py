import numpy as np

# Read the output log file
data = np.loadtxt("OUTLOG.runS", skiprows=6, delimiter="|")

# Columns 4, 5, 6 is mean position
x = np.mean(data[:,3])
y = np.mean(data[:,4])
z = np.mean(data[:,5])

# If the values are close to 0.5, assume it worked
if np.abs(x - 0.5) > 0.15 or np.abs(y - 0.5) > 0.15 or np.abs(z - 0.5) > 0.15:
    exit(1)
else:
    exit(0)
