# Open the output log file
file = open("OUTLOG.run0")

n = 0                           # Count lines
mean = 0.0                      # Accumulate 2T/VC values

# Open the output log file
#
while (line := file.readline()) != "":
    if n >= 6:                  # Skip the header stuff
        v = [float(x) for x in line.split('|')]
        mean += v[16]           # This is the 2T/VC column
    n = n + 1                   # Count lines

if n>6: mean /= n-6             # Sanity check

# Check closeness to 1.0
#
if (mean-1.0)*(mean-1.0) > 0.003:
    exit(1)
else:
    exit(0)
